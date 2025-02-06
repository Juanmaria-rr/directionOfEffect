#### 08.01.2025
#### ALL PHASES
from array import ArrayType
from functions import (
    relative_success,
    spreadSheetFormatter,
    discrepancifier,
    temporary_directionOfEffect,
)
from stoppedTrials import terminated_td
from DoEAssessment import directionOfEffect
from membraneTargets import target_membrane
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
from datetime import datetime
from datetime import date
from pyspark.sql.types import StructType, StructField, StringType, IntegerType
from pyspark.sql.types import (
    StructType,
    StructField,
    DoubleType,
    DecimalType,
    StringType,
    FloatType,
)
import pandas as pd

spark = SparkSession.builder.getOrCreate()

path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"

target = spark.read.parquet(f"{path}targets/")

diseases = spark.read.parquet(f"{path}diseases/")

evidences = spark.read.parquet(f"{path}evidence")

credible = spark.read.parquet(f"{path}credibleSet")

### index with new fix" "gs://ot-team/irene/gentropy/study_index_2412_fixed"
index = spark.read.parquet(f"gs://ot-team/irene/gentropy/study_index_2412_fixed")

new = spark.read.parquet(f"{path}colocalisation/coloc")

variantIndex = spark.read.parquet(f"{path}variantIndex")

biosample = spark.read.parquet(f"{path}biosample")


#### Fixing scXQTL as XQTLs:
## code provided by @ireneisdoomed
pd.DataFrame.iteritems = pd.DataFrame.items

raw_studies_metadata_schema: StructType = StructType(
    [
        StructField("study_id", StringType(), True),
        StructField("dataset_id", StringType(), True),
        StructField("study_label", StringType(), True),
        StructField("sample_group", StringType(), True),
        StructField("tissue_id", StringType(), True),
        StructField("tissue_label", StringType(), True),
        StructField("condition_label", StringType(), True),
        StructField("sample_size", IntegerType(), True),
        StructField("quant_method", StringType(), True),
        StructField("pmid", StringType(), True),
        StructField("study_type", StringType(), True),
    ]
)
raw_studies_metadata_path = "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/fe3c4b4ed911b3a184271a6aadcd8c8769a66aba/data_tables/dataset_metadata.tsv"

study_table = spark.createDataFrame(
    pd.read_csv(raw_studies_metadata_path, sep="\t"),
    schema=raw_studies_metadata_schema,
)

# index = spark.read.parquet("gs://open-targets-pre-data-releases/24.12-uo_test-3/output/genetics/parquet/study_index")

study_index_w_correct_type = (
    study_table.select(
        F.concat_ws(
            "_",
            F.col("study_label"),
            F.col("quant_method"),
            F.col("sample_group"),
        ).alias("extracted_column"),
        "study_type",
    )
    .join(
        index
        # Get eQTL Catalogue studies
        .filter(F.col("studyType") != "gwas").filter(
            ~F.col("studyId").startswith("UKB_PPP")
        )
        # Remove measured trait
        .withColumn(
            "extracted_column",
            F.regexp_replace(F.col("studyId"), r"(_ENS.*|_ILMN.*|_X.*|_[0-9]+:.*)", ""),
        ).withColumn(
            "extracted_column",
            # After the previous cleanup, there are multiple traits from the same publication starting with the gene symbol that need to be removed (e.g. `Sun_2018_aptamer_plasma_ANXA2.4961.17.1..1`)
            F.when(
                F.col("extracted_column").startswith("Sun_2018_aptamer_plasma"),
                F.lit("Sun_2018_aptamer_plasma"),
            ).otherwise(F.col("extracted_column")),
        ),
        on="extracted_column",
        how="right",
    )
    .persist()
)

fixed = (
    study_index_w_correct_type.withColumn(
        "toFix",
        F.when(
            (F.col("study_type") != "single-cell")
            & (F.col("studyType").startswith("sc")),
            F.lit(True),
        ).otherwise(F.lit(False)),
    )
    # Remove the substring "sc" from the studyType column
    .withColumn(
        "newStudyType",
        F.when(
            F.col("toFix"), F.regexp_replace(F.col("studyType"), r"sc", "")
        ).otherwise(F.col("studyType")),
    ).drop("toFix", "extracted_column", "study_type")
).persist()
all_studies = index.join(
    fixed.selectExpr("studyId", "newStudyType"), on="studyId", how="left"
).persist()
fixedIndex = all_studies.withColumn(
    "studyType",
    F.when(F.col("newStudyType").isNotNull(), F.col("newStudyType")).otherwise(
        F.col("studyType")
    ),
).drop("newStudyType")
#### fixed

newColoc = (
    new.join(
        credible.selectExpr(  #### studyLocusId from credible set to uncover the codified variants on left side
            "studyLocusId as leftStudyLocusId",
            "StudyId as leftStudyId",
            "variantId as leftVariantId",
            "studyType as credibleLeftStudyType",
        ),
        on="leftStudyLocusId",
        how="left",
    )
    .join(
        credible.selectExpr(  #### studyLocusId from credible set to uncover the codified variants on right side
            "studyLocusId as rightStudyLocusId",
            "studyId as rightStudyId",
            "variantId as rightVariantId",
            "studyType as credibleRightStudyType",
        ),
        on="rightStudyLocusId",
        how="left",
    )
    .join(
        fixedIndex.selectExpr(  ### bring modulated target on right side (QTL study)
            "studyId as rightStudyId",
            "geneId",
            "projectId",
            "studyType as indexStudyType",
            "condition",
            "biosampleId",
        ),
        on="rightStudyId",
        how="left",
    )
    .persist()
)
# remove columns without content (only null values on them)
df = evidences.filter((F.col("datasourceId") == "gwas_credible_sets"))

# Use an aggregation to determine non-null columns
non_null_counts = df.select(
    *[F.sum(F.col(col).isNotNull().cast("int")).alias(col) for col in df.columns]
)

# Collect the counts for each column
non_null_columns = [
    row[0] for row in non_null_counts.collect()[0].asDict().items() if row[1] > 0
]

# Select only the non-null columns
filtered_df = df.select(*non_null_columns).persist()

## bring studyId, variantId, beta from Gwas and pValue
gwasComplete = filtered_df.join(
    credible.selectExpr(
        "studyLocusId", "studyId", "variantId", "beta as betaGwas", "pValueExponent"
    ),
    on="studyLocusId",
    how="left",
).persist()

resolvedColoc = (
    (
        newColoc.withColumnRenamed("geneId", "targetId")
        .join(
            gwasComplete.withColumnRenamed("studyLocusId", "leftStudyLocusId"),
            on=["leftStudyLocusId", "targetId"],
            how="inner",
        )
        .join(  ### propagated using parent terms
            diseases.selectExpr(
                "id as diseaseId", "name", "parents", "therapeuticAreas"
            ),
            on="diseaseId",
            how="left",
        )
        .withColumn(
            "diseaseId",
            F.explode_outer(F.concat(F.array(F.col("diseaseId")), F.col("parents"))),
        )
        .drop("parents", "oldDiseaseId")
    )
    .withColumn(
        "colocDoE",
        F.when(
            F.col("rightStudyType").isin(
                ["eqtl", "pqtl", "tuqtl", "sceqtl", "sctuqtl"]
            ),
            F.when(
                (F.col("betaGwas") > 0) & (F.col("betaRatioSignAverage") > 0),
                F.lit("GoF_risk"),
            )
            .when(
                (F.col("betaGwas") > 0) & (F.col("betaRatioSignAverage") < 0),
                F.lit("LoF_risk"),
            )
            .when(
                (F.col("betaGwas") < 0) & (F.col("betaRatioSignAverage") > 0),
                F.lit("LoF_protect"),
            )
            .when(
                (F.col("betaGwas") < 0) & (F.col("betaRatioSignAverage") < 0),
                F.lit("GoF_protect"),
            ),
        ).when(
            F.col("rightStudyType").isin(
                ["sqtl", "scsqtl"]
            ),  ### opposite directionality than sqtl
            F.when(
                (F.col("betaGwas") > 0) & (F.col("betaRatioSignAverage") > 0),
                F.lit("LoF_risk"),
            )
            .when(
                (F.col("betaGwas") > 0) & (F.col("betaRatioSignAverage") < 0),
                F.lit("GoF_risk"),
            )
            .when(
                (F.col("betaGwas") < 0) & (F.col("betaRatioSignAverage") > 0),
                F.lit("GoF_protect"),
            )
            .when(
                (F.col("betaGwas") < 0) & (F.col("betaRatioSignAverage") < 0),
                F.lit("LoF_protect"),
            ),
        ),
    )
    .persist()
)

path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"

datasource_filter = [
    "ot_genetics_portal",
    "gwas_credible_sets",
    "gene_burden",
    "eva",
    "eva_somatic",
    "gene2phenotype",
    "orphanet",
    "cancer_gene_census",
    "intogen",
    "impc",
    "chembl",
]

assessment, evidences, actionType, oncolabel = temporary_directionOfEffect(
    path, datasource_filter
)

drugApproved = (
    spark.read.parquet("gs://ot-team/irene/l2g/validation/chembl_w_flags")
    .drop("clinicalTrialId", "isComplex")
    .withColumn(
        "isApproved",
        F.when(F.col("isApproved") == "true", F.lit(1)).otherwise(F.lit(0)),
    )
    .distinct()
)


def comparisons_df_iterative(disdic, projectId):
    toAnalysis = [(key, value) for key, value in disdic.items() if value == projectId]
    schema = StructType(
        [
            StructField("comparison", StringType(), True),
            StructField("comparisonType", StringType(), True),
        ]
    )

    comparisons = spark.createDataFrame(toAnalysis, schema=schema)
    ### include all the columns as predictor

    predictions = spark.createDataFrame(
        data=[
            ("Phase4", "clinical"),
            ("Phase>=3", "clinical"),
            ("Phase>=2", "clinical"),
            ("Phase>=1", "clinical"),
            # ("nPhase4", "clinical"),
            # ("nPhase>=3", "clinical"),
            # ("nPhase>=2", "clinical"),
            # ("nPhase>=1", "clinical"),
            ("approved", "clinical"),
            # ("PhaseT", "clinical"),
        ]
    )
    return comparisons.join(predictions, how="full").collect()


full_data = spark.createDataFrame(
    data=[
        ("yes", "yes"),
        ("yes", "no"),
        ("no", "yes"),
        ("no", "no"),
    ],
    schema=StructType(
        [
            StructField("prediction", StringType(), True),
            StructField("comparison", StringType(), True),
        ]
    ),
)

chembl_indication = (
    discrepancifier(
        (
            assessment.filter(
                (F.col("datasourceId") == "chembl")
                & (F.col("homogenized") != "noEvaluable")
            ).join(
                drugApproved.filter(F.col("isApproved") == 1),
                on=["targetId", "diseaseId", "drugId"],
                how="left",
            )
        )
        .withColumn(
            "approvedDrug",
            F.max(F.col("isApproved")).over(
                Window.partitionBy("targetId", "diseaseId", "drugId")
            ),
        )
        .groupBy(
            "targetId",
            "diseaseId",
            "studyId",
            "drugId",
            "clinicalPhase",
            "approvedDrug",
        )
        .pivot("homogenized")
        .count()
    )
    .withColumnRenamed("studyId", "clinicalStudyId")
    .filter(F.col("coherencyDiagonal") == "coherent")
    .drop(
        "coherencyDiagonal", "coherencyOneCell", "noEvaluable", "GoF_risk", "LoF_risk"
    )
    .withColumnRenamed("GoF_protect", "drugGoF_protect")
    .withColumnRenamed("LoF_protect", "drugLoF_protect")
)

new_benchmark = (
    (
        resolvedColoc.filter(F.col("betaGwas") < 0)
        .join(  ### select just GWAS giving protection
            chembl_indication, on=["targetId", "diseaseId"], how="inner"
        )
        .withColumn(
            "AgreeDrug",
            F.when(
                (F.col("drugGoF_protect").isNotNull())
                & (F.col("colocDoE") == "GoF_protect"),
                F.lit("yes"),
            )
            .when(
                (F.col("drugLoF_protect").isNotNull())
                & (F.col("colocDoE") == "LoF_protect"),
                F.lit("yes"),
            )
            .otherwise(F.lit("no")),
        )
    ).filter(
        F.col("name") != "COVID-19"
    )  #### remove COVID-19 associations
).join(biosample.select("biosampleId", "biosampleName"), on="biosampleId", how="left")

import numpy as np


def convertTuple(tup):
    st = ",".join(map(str, tup))
    return st


from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from pyspark.sql.types import *

#### columsn to groupBy - introduce a dictionary for trying different list
partitionByPValue = ["targetId", "diseaseId", "rightStudyType"]

#### description of groups to make groupBy
## partitionForPhase
groupPhase1 = ["targetId", "diseaseId"]
groupPhase2 = ["targetId", "diseaseId", "clinicalPhase"]
groupPhase3 = ["targetId", "diseaseId", "clinicalStudyId"]
groupPhase4 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "clinicalStudyId",
]  ## will take each clinical phase and study
groupPhase5 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "clinicalStudyId",
    "drugId",
]  ## will take max clin phase for DrugId, clinical phase and clnical study
groupPhase6 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "drugId",
]  ## will take each max clinical phase per drugId
groupPhase7 = ["targetId", "diseaseId", "drugId"]

group_phases = {
    "targetIdDiseaseID": groupPhase1,
    "clinicalPhase": groupPhase2,
    "clinicalStudy": groupPhase3,
    "clinicalPhase&ClinicalStudy": groupPhase4,
    "clinicalPhaseClinicalStudyIdDrugId": groupPhase5,
    "clinicalPhaseDrugId": groupPhase6,
    "drugId": groupPhase7,
}

## partitionForRows
groupPhase1 = ["targetId", "diseaseId", "clinicalPhase"]
groupPhase2 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "clinicalStudyId",
]  ## will take each clinical phase and study
groupPhase3 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "clinicalStudyId",
    "drugId",
]  ## will take max clin phase for DrugId, clinical phase and clnical study
groupPhase4 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "drugId",
]  ## will take each max clinical phase per drugId
# same groups but containing approved
groupPhase5 = ["targetId", "diseaseId", "clinicalPhase", "approved_l"]
groupPhase6 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "clinicalStudyId",
    "approved_l",
]  ## will take each clinical phase and study
groupPhase7 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "clinicalStudyId",
    "drugId",
    "approved_l",
]  ## will take max clin phase for DrugId, clinical phase and clnical study
groupPhase8 = [
    "targetId",
    "diseaseId",
    "clinicalPhase",
    "drugId",
    "approved_l",
]  ## will take each max clinical phase per drugId

group_rows = {
    "clinicalPhase": groupPhase1,
    "clinicalPhaseClinicalStudyId": groupPhase2,
    "clinicalPhaseClinicalStudyIdDrugId": groupPhase3,
    "clinicalPhaseDrugId": groupPhase4,
    "clinicalPhaseApproved_l": groupPhase5,
    "clinicalPhaseClinicalStudyIdApproved_l": groupPhase6,
    "clinicalPhaseClinicalStudyIdDrugIdApproved_l": groupPhase7,
    "clinicalPhaseDrugIdApproved_l": groupPhase8,
}

value_analysis = ["pqtl", "eqtl", "sqtl"]
### define the window to order for taking pValue
window_spec = Window.partitionBy(*partitionByPValue).orderBy(
    F.col("pValueExponent").asc()
)
results = []

for value in value_analysis:
    # Iterate over group mapping
    for group_phase, group_phase_columns in group_phases.items():
        print("making group phases:", group_phase_columns)
        for group_row, group_rows_columns in group_rows.items():
            print("making group rows:", group_rows_columns)
            x = value
            print(value, x)

            if "approved_l" in group_rows_columns:

                pre = (
                    new_benchmark.withColumn(
                        "approved_l",
                        F.when(F.col("approvedDrug") == 1, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .filter(F.col("clinicalStudyId").isNotNull())
                    .withColumn(
                        "agree_lowestPval",
                        F.first("AgreeDrug", ignorenulls=True).over(
                            window_spec
                        ),  ### ignoring null values
                    )
                    .withColumn(
                        "clinicalPhase",  ### no longer maxclinphase for T-D
                        F.max("clinicalPhase").over(
                            Window.partitionBy(*group_phase_columns)
                        ),
                    )
                    .withColumn(
                        "approved",  ### no longer maxclinphase for T-D
                        F.max("approvedDrug").over(
                            Window.partitionBy(*group_phase_columns)
                        ),
                    )
                    .groupBy(*group_rows_columns)
                    .pivot("rightStudyType")
                    .agg(F.collect_set("agree_lowestPVal"))
                    .withColumn(
                        "Phase4",
                        F.when(F.col("clinicalPhase") == 4, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        "Phase>=3",
                        F.when(F.col("clinicalPhase") >= 3, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        "Phase>=2",
                        F.when(F.col("clinicalPhase") >= 2, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        "Phase>=1",
                        F.when(F.col("clinicalPhase") >= 1, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        "approved",
                        F.when(F.col("approved_l") == 1, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        x,
                        F.when(
                            F.array_contains(F.col(x), "yes"), F.lit("yes")
                        ).otherwise(F.lit("no")),
                    )
                )
                for phase in ["Phase4", "approved"]:
                    print("groupby for", phase)
                    pre1 = (
                        pre.select(
                            F.col(x).alias("comparison"),
                            F.col(phase).alias("prediction"),
                        )
                        .join(full_data, on=["prediction", "comparison"], how="outer")
                        .groupBy("comparison")
                        .pivot("prediction")
                        .count()
                        .select(F.col("comparison").alias(x), "yes", "no")
                        .sort(F.col(x).desc())
                    )
                    array1 = np.delete(
                        (pre1).fillna(0).toPandas().to_numpy(),
                        [0],
                        1,
                    )
                    total = np.sum(array1)
                    res_npPhaseX = np.array(array1, dtype=int)
                    resX = convertTuple(
                        fisher_exact(res_npPhaseX, alternative="two-sided")
                    )
                    resx_CI = convertTuple(
                        odds_ratio(res_npPhaseX).confidence_interval(
                            confidence_level=0.95
                        )
                    )
                    print(
                        round(float(resX.split(",")[0]), 2),
                        float(resX.split(",")[1]),
                        round(float(resx_CI.split(",")[0]), 2),
                        round(float(resx_CI.split(",")[1]), 2),
                    )
                    print("\n")
                    results.append(
                        [
                            partitionByPValue,
                            group_phase,
                            group_row,
                            phase,
                            x,
                            round(float(resX.split(",")[0]), 2),  ## OR
                            float(resX.split(",")[1]),  ## pValue
                            round(float(resx_CI.split(",")[0]), 2),  ## Low CI
                            round(float(resx_CI.split(",")[1]), 2),  ## High CI
                            total,
                            array1,
                        ]
                    )
            else:
                pre = (
                    new_benchmark.withColumn(
                        "approved_l",
                        F.when(F.col("approvedDrug") == 1, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .filter(F.col("clinicalStudyId").isNotNull())
                    .withColumn(
                        "agree_lowestPval",
                        F.first("AgreeDrug", ignorenulls=True).over(window_spec),
                    )
                    .withColumn(
                        "clinicalPhase",  ### no longer maxclinphase for T-D
                        F.max("clinicalPhase").over(
                            Window.partitionBy(*group_phase_columns)
                        ),
                    )
                    .groupBy(*group_rows_columns)
                    .pivot("rightStudyType")
                    .agg(F.collect_set("agree_lowestPVal"))
                    .withColumn(
                        "Phase4",
                        F.when(F.col("clinicalPhase") == 4, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        "Phase>=3",
                        F.when(F.col("clinicalPhase") >= 3, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        "Phase>=2",
                        F.when(F.col("clinicalPhase") >= 2, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        "Phase>=1",
                        F.when(F.col("clinicalPhase") >= 1, F.lit("yes")).otherwise(
                            F.lit("no")
                        ),
                    )
                    .withColumn(
                        x,
                        F.when(
                            F.array_contains(F.col(x), "yes"), F.lit("yes")
                        ).otherwise(F.lit("no")),
                    )
                )
                for phase in ["Phase4"]:
                    print("groupby for", phase)
                    # pre.groupBy(x).pivot(phase).count().select(x,"yes","no").sort(F.col(x).desc()).show()
                    pre1 = (
                        pre.select(
                            F.col(x).alias("comparison"),
                            F.col(phase).alias("prediction"),
                        )
                        .join(full_data, on=["prediction", "comparison"], how="outer")
                        .groupBy("comparison")
                        .pivot("prediction")
                        .count()
                        .select(F.col("comparison").alias(x), "yes", "no")
                        .sort(F.col(x).desc())
                    )
                    array1 = np.delete(
                        (pre1).fillna(0).toPandas().to_numpy(),
                        [0],
                        1,
                    )
                    total = np.sum(array1)
                    res_npPhaseX = np.array(array1, dtype=int)
                    resX = convertTuple(
                        fisher_exact(res_npPhaseX, alternative="two-sided")
                    )
                    resx_CI = convertTuple(
                        odds_ratio(res_npPhaseX).confidence_interval(
                            confidence_level=0.95
                        )
                    )
                    print(
                        round(float(resX.split(",")[0]), 2),
                        float(resX.split(",")[1]),
                        round(float(resx_CI.split(",")[0]), 2),
                        round(float(resx_CI.split(",")[1]), 2),
                    )
                    print("\n")
                    results.append(
                        [
                            partitionByPValue,
                            group_phase,
                            group_row,
                            phase,
                            x,
                            round(float(resX.split(",")[0]), 2),  ## OR
                            float(resX.split(",")[1]),  ## pValue
                            round(float(resx_CI.split(",")[0]), 2),  ## Low CI
                            round(float(resx_CI.split(",")[1]), 2),  ## High CI
                            total,
                            array1,
                        ]
                    )
df = pd.DataFrame(
    results,
    columns=[
        "partitionByPValue",
        "partitionForPhase",
        "groupByForRows",
        "phase",
        "x",
        "OR",
        "pValue",
        "LowCI",
        "HighCI",
        "total",
        "array",
    ],
)
print("created dataframe")
df.to_csv("gs://ot-team/jroldan/analysis/allCombinations_ignoringnulls_good.csv")
print("dataframe saved")

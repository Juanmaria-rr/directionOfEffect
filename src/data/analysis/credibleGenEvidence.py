#### 11.12.2024
#######
##########     ATENTION
#### change code to work with generated dataframe instead of reading the parquet

"""
This scripts run Odds ratio analysis for DoE and 
genetic information on drug clinical success

"""
from functions import (
    discrepancifier,
    build_gwasResolvedColoc,
    temporary_directionOfEffect,
    spreadSheetFormatter,
)

from pyspark.sql.functions import regexp_extract  # type: ignore
from pyspark.sql import SparkSession, Window  # type: ignore
import pyspark.sql.functions as F
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from pyspark.sql.types import StructType, StructField, StringType, IntegerType
from pyspark.sql.types import (
    StructType,
    StructField,
    ArrayType,
    DoubleType,
    DecimalType,
    StringType,
    FloatType,
)
from datetime import datetime


spark = SparkSession.builder.getOrCreate()
c = datetime.now()
today_date = str(date.today())
print("spark session created at", c)

print("Analysis started on " + today_date + " at ", c)
"""
#coloc = spark.read.parquet(
#    "gs://genetics-portal-dev-data/22.09.1/outputs/v2d_coloc"
#).filter(F.col("right_type") != "gwas")
"""

#### make the dataset from stopped clin trials
### read supplementary table 9
""" ### just showing how i did the dataset
st9 = spark.read.csv("/Users/juanr/Downloads/ST9.csv", sep=",", header=True)
st9.filter(
    (F.col("clinicalStatus").isin(["Terminated", "Withdrawn", "Suspended"]))
    & (F.col("prediction") == "Negative")
).groupBy(
    "targetId", "diseaseId", "clinicalStatus", "prediction"
).count().toPandas().to_csv(
    "targetDiseaseStoppedNegative.csv"
)
"""
### target-diseases terminated&withdrawal in clin trials
terminated = spark.read.csv(
    "gs://ot-team/jroldan/analysis/targetDiseaseStoppedNegative.csv",
    sep=",",
    header=True,
).drop("_c0", "Withdrawn")

path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"

evidences = (
    spark.read.parquet(f"{path}evidence").filter(
        F.col("datasourceId").isin(
            [
                "ot_genetics_portal",
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
        )
    )
    # .persist()
)
ot_genetics = evidences.filter(F.col("datasourceId") == "ot_genetics_portal")

#### Now load sources of data to generate credible_set_OT_genetics evidences and associations.

target = spark.read.parquet(f"{path}targets/")

diseases = spark.read.parquet(f"{path}diseases/")

gwasResolvedColoc = build_gwasResolvedColoc(path)

#### take the direction from the lowest p value
window_spec = Window.partitionBy("targetId", "diseaseId").orderBy(
    F.col("pValueExponent").asc()
)
gwasCredibleAssoc = (
    gwasResolvedColoc.withColumn(
        "homogenized", F.first("colocDoE", ignorenulls=True).over(window_spec)
    )  ## added 30.01.2025
    .select("targetId", "diseaseId", "homogenized")
    .withColumn(
        "homogenized",
        F.when(F.col("homogenized").isNull(), F.lit("noEvaluable")).otherwise(
            F.col("homogenized")
        ),
    )
)

### datasource_filter for temporaryDoE
datasource_filter = [
    "ot_genetics_portal",
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

prueba_assessment, evidences, actionType, oncolabel = temporary_directionOfEffect(
    path, datasource_filter
)

print("Moving to step 2")

columns_chembl = ["LoF_protect", "GoF_protect"]
columns_dataset = ["LoF_protect", "GoF_protect", "LoF_risk", "GoF_risk", "evidenceDif"]
columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect"]
terms = ["noEvaluable", "bivalent_risk", "null", "dispar"]

taDf = spark.createDataFrame(
    data=[
        ("MONDO_0045024", "cell proliferation disorder", "Oncology"),
        ("EFO_0005741", "infectious disease", "Other"),
        ("OTAR_0000014", "pregnancy or perinatal disease", "Other"),
        ("EFO_0005932", "animal disease", "Other"),
        ("MONDO_0024458", "disease of visual system", "Other"),
        ("EFO_0000319", "cardiovascular disease", "Other"),
        ("EFO_0009605", "pancreas disease", "Other"),
        ("EFO_0010282", "gastrointestinal disease", "Other"),
        ("OTAR_0000017", "reproductive system or breast disease", "Other"),
        ("EFO_0010285", "integumentary system disease", "Other"),
        ("EFO_0001379", "endocrine system disease", "Other"),
        ("OTAR_0000010", "respiratory or thoracic disease", "Other"),
        ("EFO_0009690", "urinary system disease", "Other"),
        ("OTAR_0000006", "musculoskeletal or connective tissue disease", "Other"),
        ("MONDO_0021205", "disease of ear", "Other"),
        ("EFO_0000540", "immune system disease", "Other"),
        ("EFO_0005803", "hematologic disease", "Other"),
        ("EFO_0000618", "nervous system disease", "Other"),
        ("MONDO_0002025", "psychiatric disorder", "Other"),
        ("MONDO_0024297", "nutritional or metabolic disease", "Other"),
        ("OTAR_0000018", "genetic, familial or congenital disease", "Other"),
        ("OTAR_0000009", "injury, poisoning or other complication", "Other"),
        ("EFO_0000651", "phenotype", "Other"),
        ("EFO_0001444", "measurement", "Other"),
        ("GO_0008150", "biological process", "Other"),
    ],
    schema=StructType(
        [
            StructField("taId", StringType(), True),
            StructField("taLabel", StringType(), True),
            StructField("taLabelSimple", StringType(), True),
        ]
    ),
).withColumn("taRank", F.monotonically_increasing_id())

### give us a classification of Oncology VS non oncology
wByDisease = Window.partitionBy("diseaseId")  #### checked 31.05.2023
diseaseTA = (
    diseases.withColumn("taId", F.explode("therapeuticAreas"))
    .select(F.col("id").alias("diseaseId"), "taId", "parents")
    .join(taDf, on="taId", how="left")
    .withColumn("minRank", F.min("taRank").over(wByDisease))
    .filter(F.col("taRank") == F.col("minRank"))
    .drop("taRank", "minRank")
)

#### give us propagation of diseases and list of therapeutic areas associated
diseases2 = diseases.select("id", "parents").withColumn(
    "diseaseIdPropagated",
    F.explode_outer(F.concat(F.array(F.col("id")), F.col("parents"))),
)

chembl_trials = (
    prueba_assessment.filter((F.col("datasourceId").isin(["chembl"])))
    .groupBy("targetId", "diseaseId")
    .agg(F.max(F.col("clinicalPhase")).alias("maxClinPhase"))
)

terminated_array = (
    terminated.groupBy("targetId", "diseaseId")
    .agg(F.collect_set("clinicalStatus").alias("clinicalStatus"))
    .withColumn("prediction", F.when(F.col("clinicalStatus").isNotNull(), F.lit("yes")))
)

assessment = prueba_assessment.unionByName(
    gwasCredibleAssoc.withColumn("datasourceId", F.lit("gwas_credible_set")),
    allowMissingColumns=True,
)

print("defining non propagated,propagated and analysis_drugs functions")


def analysis_nonPropagated(assessment, analysisDatasources):
    return discrepancifier(
        assessment.filter(F.col("datasourceId").isin(analysisDatasources))
        .withColumn(
            "datasources",
            F.collect_set(F.col("datasourceId")).over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .groupBy(
            "targetId",
            "diseaseId",
        )
        .pivot("homogenized")
        .agg(F.count("targetId"))
        # .persist()
    )


def analysis_propagated(assessment, analysisDatasources):
    return discrepancifier(
        assessment.filter(F.col("datasourceId").isin(analysisDatasources))
        .withColumn(
            "datasources",
            F.collect_set(F.col("datasourceId")).over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .join(
            diseases2.selectExpr("id as diseaseId", "diseaseIdPropagated"),
            on="diseaseId",
            how="left",
        )
        .withColumnRenamed("diseaseId", "oldDiseaseId")
        .withColumnRenamed("diseaseIdPropagated", "diseaseId")
        .groupBy(
            "targetId",
            "diseaseId",
        )
        .pivot("homogenized")
        .agg(F.count("targetId"))
        # .persist()
    )


chembl_ds = ["chembl"]
drugApproved = (  ### added on 30.01.2025
    spark.read.parquet("gs://ot-team/irene/l2g/validation/chembl_w_flags")
    .drop("clinicalTrialId", "isComplex")
    .withColumn(
        "isApproved",
        F.when(F.col("isApproved") == "true", F.lit(1)).otherwise(F.lit(0)),
    )
    .distinct()
)


def analysis_drugs(assessment, chembl_ds, drugApproved):
    return discrepancifier(
        assessment.filter((F.col("datasourceId").isin(chembl_ds)))
        .join(
            drugApproved.filter(F.col("isApproved") == 1),
            on=["targetId", "diseaseId", "drugId"],
            how="left",
        )
        .withColumn(
            "maxClinPhase",
            F.max(F.col("clinicalPhase")).over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .withColumn(
            "approvedDrug",
            F.max(F.col("isApproved")).over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .withColumn(
            "approved",
            F.when(F.col("approvedDrug") == 1, F.lit("yes")).otherwise(F.lit("no")),
        )
        .groupBy("targetId", "diseaseId", "maxClinPhase", "approved")
        .pivot("homogenized")
        .agg(F.count("targetId"))
        .persist()
    )


analysis_chembl = analysis_drugs(assessment, chembl_ds, drugApproved)

#######
## include here the analysis
#######

analysisDatasources = []

print("defining full_analysis_propagation")


def full_analysis_propagation(
    assessment, analysisDatasources, analysis_chembl, terminated_array, diseaseTA
):
    return (
        analysis_propagated(assessment, analysisDatasources)
        .join(
            analysis_chembl.selectExpr(
                "targetId",
                "diseaseId",
                "maxClinPhase",
                "approved",
                "coherencyDiagonal as coherencyDiagonal_ch",
                "coherencyOneCell as coherencyOneCell_ch",
                "LoF_protect as LoF_protect_ch",
                "GoF_protect as GoF_protect_ch",
            ),
            on=["targetId", "diseaseId"],
            how="right",
        )
        #### Should remove the coherencyDiagonal.isNotNull()
        .withColumn(
            "geneticEvidence",
            F.when(
                F.col("coherencyDiagonal").isNotNull(), F.lit("hasGeneticEvidence")
            ).otherwise(F.lit("noGeneticEvidence")),
        )
        # .filter(F.col("coherencyDiagonal_ch").isNotNull())
        .withColumn(
            "diagonalAgreeWithDrugs",
            F.when(
                (F.col("coherencyDiagonal_ch") == "coherent")
                & (F.col("coherencyDiagonal") == "coherent"),
                F.when(
                    (F.col("LoF_protect_ch").isNotNull())
                    & (
                        F.col("GoF_risk").isNotNull() | F.col("LoF_protect").isNotNull()
                    ),
                    F.lit("coherent"),
                )
                .when(
                    F.col("GoF_protect_ch").isNotNull()
                    & (
                        F.col("LoF_risk").isNotNull() | F.col("GoF_protect").isNotNull()
                    ),
                    F.lit("coherent"),
                )
                .otherwise(F.lit("dispar")),
            ),
        )
        .withColumn(
            "oneCellAgreeWithDrugs",
            F.when(
                (F.col("coherencyOneCell_ch") == "coherent")
                & (F.col("coherencyOneCell") == "coherent"),
                F.when(
                    (F.col("LoF_protect_ch").isNotNull())
                    & (
                        (F.col("LoF_protect").isNotNull())
                        & (F.col("LoF_risk").isNull())
                        & (F.col("GoF_protect").isNull())
                        & (F.col("GoF_risk").isNull())
                    ),
                    F.lit("coherent"),
                )
                .when(
                    (F.col("GoF_protect_ch").isNotNull())
                    & (
                        (F.col("GoF_protect").isNotNull())
                        & (F.col("LoF_risk").isNull())
                        & (F.col("LoF_protect").isNull())
                        & (F.col("GoF_risk").isNull())
                    ),
                    F.lit("coherent"),
                )
                .otherwise(F.lit("dispar")),
            ),
        )
        .withColumn(
            "Phase4",
            F.when(F.col("maxClinPhase") == 4, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "Phase>=3",
            F.when(F.col("maxClinPhase") >= 3, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "Phase>=2",
            F.when(F.col("maxClinPhase") >= 2, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "Phase>=1",
            F.when(F.col("maxClinPhase") >= 1, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "Phase0",
            F.when(F.col("maxClinPhase") == 0, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "approved",
            F.when(F.col("approved") == "yes", F.lit("yes")).otherwise(F.lit("no")),
        )
        .join(terminated_array, on=["targetId", "diseaseId"], how="left")
        .withColumn(
            "PhaseT",
            F.when(F.col("prediction") == "yes", F.lit("yes")).otherwise(F.lit("no")),
        )
        .join(
            diseaseTA.select("diseaseId", "taLabelSimple"), on="diseaseId", how="left"
        )
        .withColumn(
            "hasGeneticEvidence",
            F.when(
                F.col("geneticEvidence") == "hasGeneticEvidence", F.lit("yes")
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "diagonalYes",
            F.when(
                F.col("hasGeneticEvidence") == "yes",
                F.when(F.col("diagonalAgreeWithDrugs") == "coherent", F.lit("yes"))
                .when(F.col("diagonalAgreeWithDrugs") == "dispar", F.lit("no"))
                .otherwise(F.lit("no")),
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "oneCellYes",
            F.when(
                F.col("hasGeneticEvidence") == "yes",
                F.when(F.col("oneCellAgreeWithDrugs") == "coherent", F.lit("yes"))
                .when(F.col("oneCellAgreeWithDrugs") == "dispar", F.lit("no"))
                .otherwise(F.lit("no")),
            ).otherwise(F.lit("no")),
        )
        # .persist()
    )


#####
## no propag
#####
print("defining full analysis no propagation")


def full_analysis_noPropagation(
    assessment, analysisDatasources, analysis_chembl, terminated_array, diseaseTA
):
    return (
        analysis_nonPropagated(assessment, analysisDatasources)
        .join(
            analysis_chembl.selectExpr(
                "targetId",
                "diseaseId",
                "maxClinPhase",
                "approved",
                "coherencyDiagonal as coherencyDiagonal_ch",
                "coherencyOneCell as coherencyOneCell_ch",
                "LoF_protect as LoF_protect_ch",
                "GoF_protect as GoF_protect_ch",
            ),
            on=["targetId", "diseaseId"],
            how="right",
        )
        .withColumn(
            "geneticEvidence",
            F.when(
                F.col("coherencyDiagonal").isNotNull(), F.lit("hasGeneticEvidence")
            ).otherwise(F.lit("noGeneticEvidence")),
        )
        # .filter(F.col("coherencyDiagonal_ch").isNotNull())
        .withColumn(
            "diagonalAgreeWithDrugs",
            F.when(
                (F.col("coherencyDiagonal_ch") == "coherent")
                & (F.col("coherencyDiagonal") == "coherent"),
                F.when(
                    (F.col("LoF_protect_ch").isNotNull())
                    & (
                        F.col("GoF_risk").isNotNull() | F.col("LoF_protect").isNotNull()
                    ),
                    F.lit("coherent"),
                )
                .when(
                    F.col("GoF_protect_ch").isNotNull()
                    & (
                        F.col("LoF_risk").isNotNull() | F.col("GoF_protect").isNotNull()
                    ),
                    F.lit("coherent"),
                )
                .otherwise(F.lit("dispar")),
            ),
        )
        .withColumn(
            "oneCellAgreeWithDrugs",
            F.when(
                (F.col("coherencyOneCell_ch") == "coherent")
                & (F.col("coherencyOneCell") == "coherent"),
                F.when(
                    (F.col("LoF_protect_ch").isNotNull())
                    & (
                        (F.col("LoF_protect").isNotNull())
                        & (F.col("LoF_risk").isNull())
                        & (F.col("GoF_protect").isNull())
                        & (F.col("GoF_risk").isNull())
                    ),
                    F.lit("coherent"),
                )
                .when(
                    (F.col("GoF_protect_ch").isNotNull())
                    & (
                        (F.col("GoF_protect").isNotNull())
                        & (F.col("LoF_risk").isNull())
                        & (F.col("LoF_protect").isNull())
                        & (F.col("GoF_risk").isNull())
                    ),
                    F.lit("coherent"),
                )
                .otherwise(F.lit("dispar")),
            ),
        )
        .withColumn(
            "Phase4",
            F.when(F.col("maxClinPhase") == 4, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "Phase>=3",
            F.when(F.col("maxClinPhase") >= 3, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "Phase>=2",
            F.when(F.col("maxClinPhase") >= 2, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "Phase>=1",
            F.when(F.col("maxClinPhase") >= 1, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "Phase0",
            F.when(F.col("maxClinPhase") == 0, F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "approved",
            F.when(F.col("approved") == "yes", F.lit("yes")).otherwise(F.lit("no")),
        )
        .join(terminated_array, on=["targetId", "diseaseId"], how="left")
        .withColumn(
            "PhaseT",
            F.when(F.col("prediction") == "yes", F.lit("yes")).otherwise(F.lit("no")),
        )
        .join(
            diseaseTA.select("diseaseId", "taLabelSimple"), on="diseaseId", how="left"
        )
        .withColumn(
            "hasGeneticEvidence",
            F.when(
                F.col("geneticEvidence") == "hasGeneticEvidence", F.lit("yes")
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "diagonalYes",
            F.when(
                F.col("hasGeneticEvidence") == "yes",
                F.when(F.col("diagonalAgreeWithDrugs") == "coherent", F.lit("yes"))
                .when(F.col("diagonalAgreeWithDrugs") == "dispar", F.lit("no"))
                .otherwise(F.lit("no")),
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "oneCellYes",
            F.when(
                F.col("hasGeneticEvidence") == "yes",
                F.when(F.col("oneCellAgreeWithDrugs") == "coherent", F.lit("yes"))
                .when(F.col("oneCellAgreeWithDrugs") == "dispar", F.lit("no"))
                .otherwise(F.lit("no")),
            ).otherwise(F.lit("no")),
        )
        # .persist()
    )


print("moving to Step 3")

from functions import relative_success, spreadSheetFormatter, convertTuple
import re
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio, relative_risk

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
c = datetime.now()
print("starting dictionaries at", c)

#### continue here on 10.07.2024

## 1nd dictionary
dfs_dict = {}  ### checked and changed on 01.06.2023
dfs_dict_propag = {}


wocgc_list = [
    "gene_burden",
    "intogen",
    "eva",
    "eva_somatic",
    # "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
    "gwas_credible_set",
]
wCgc_list = [
    "gene_burden",
    "intogen",
    "eva",
    "eva_somatic",
    # "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
    "gwas_credible_set",
    "cancer_gene_census",
]

datasource_list = [
    "gene_burden",
    "intogen",
    "cancer_gene_census",
    "eva",
    "eva_somatic",
    "ot_genetics_portal",
    "gwas_credible_set",
    "impc",
    "orphanet",
    "gene2phenotype",
    "WOcgc",
    "wCgc",
    "somatic",
    "germline",
]

germline_list = [
    "gene_burden",
    "eva",
    # "ot_genetics_portal",
    "gwas_credible_set",
    "impc",
    "orphanet",
    "gene2phenotype",
]

somatic_list = ["intogen", "cancer_gene_census", "eva_somatic"]

# assessment = prueba_assessment.filter(F.col("datasourceId").isin(datasources_analysis))


def dataset_builder(assessment, value, analysis_chembl, terminated_array, diseaseTA):
    nonPropagated = full_analysis_noPropagation(
        assessment, value, analysis_chembl, terminated_array, diseaseTA
    )
    propagated = full_analysis_propagation(
        assessment, value, analysis_chembl, terminated_array, diseaseTA
    )
    return (
        # Non propagation
        ## All
        nonPropagated,
        ## Other
        nonPropagated.filter(F.col("taLabelSimple") == "Other"),
        ## Other&Null
        nonPropagated.filter(
            (F.col("taLabelSimple").isNull()) | (F.col("taLabelSimple") == "Other")
        ),
        ## Oncology
        nonPropagated.filter(F.col("taLabelSimple") == "Oncology"),
        # Propagation
        ## All
        propagated,
        ## Other
        propagated.filter(F.col("taLabelSimple") == "Other"),
        ## Other&Null
        propagated.filter(
            (F.col("taLabelSimple").isNull()) | (F.col("taLabelSimple") == "Other")
        ),
        ## Oncology
        propagated.filter(F.col("taLabelSimple") == "Oncology"),
    )


for value in datasource_list:
    print(value)
    if value == "WOcgc":
        (
            dfs_dict[f"df_{value}_All_original"],
            dfs_dict[f"df_{value}_Other_original"],
            dfs_dict[f"df_{value}_OtherNull_original"],
            dfs_dict[f"df_{value}_Oncology_original"],
            dfs_dict_propag[f"df_{value}_All_propag"],
            dfs_dict_propag[f"df_{value}_Other_propag"],
            dfs_dict_propag[f"df_{value}_OtherNull_propag"],
            dfs_dict_propag[f"df_{value}_Oncology_propag"],
        ) = dataset_builder(
            assessment, wocgc_list, analysis_chembl, terminated_array, diseaseTA
        )
    elif value == "wCgc":
        (
            dfs_dict[f"df_{value}_All_original"],
            dfs_dict[f"df_{value}_Other_original"],
            dfs_dict[f"df_{value}_OtherNull_original"],
            dfs_dict[f"df_{value}_Oncology_original"],
            dfs_dict_propag[f"df_{value}_All_propag"],
            dfs_dict_propag[f"df_{value}_Other_propag"],
            dfs_dict_propag[f"df_{value}_OtherNull_propag"],
            dfs_dict_propag[f"df_{value}_Oncology_propag"],
        ) = dataset_builder(
            assessment, wCgc_list, analysis_chembl, terminated_array, diseaseTA
        )
    elif value == "germline":
        (
            dfs_dict[f"df_{value}_All_original"],
            dfs_dict[f"df_{value}_Other_original"],
            dfs_dict[f"df_{value}_OtherNull_original"],
            dfs_dict[f"df_{value}_Oncology_original"],
            dfs_dict_propag[f"df_{value}_All_propag"],
            dfs_dict_propag[f"df_{value}_Other_propag"],
            dfs_dict_propag[f"df_{value}_OtherNull_propag"],
            dfs_dict_propag[f"df_{value}_Oncology_propag"],
        ) = dataset_builder(
            assessment,
            germline_list,
            analysis_chembl,
            terminated_array,
            diseaseTA,
        )

    elif value == "somatic":
        (
            dfs_dict[f"df_{value}_All_original"],
            dfs_dict[f"df_{value}_Other_original"],
            dfs_dict[f"df_{value}_OtherNull_original"],
            dfs_dict[f"df_{value}_Oncology_original"],
            dfs_dict_propag[f"df_{value}_All_propag"],
            dfs_dict_propag[f"df_{value}_Other_propag"],
            dfs_dict_propag[f"df_{value}_OtherNull_propag"],
            dfs_dict_propag[f"df_{value}_Oncology_propag"],
        ) = dataset_builder(
            assessment,
            somatic_list,
            analysis_chembl,
            terminated_array,
            diseaseTA,
        )

    else:
        (
            dfs_dict[f"df_{value}_All_original"],
            dfs_dict[f"df_{value}_Other_original"],
            dfs_dict[f"df_{value}_OtherNull_original"],
            dfs_dict[f"df_{value}_Oncology_original"],
            dfs_dict_propag[f"df_{value}_All_propag"],
            dfs_dict_propag[f"df_{value}_Other_propag"],
            dfs_dict_propag[f"df_{value}_OtherNull_propag"],
            dfs_dict_propag[f"df_{value}_Oncology_propag"],
        ) = dataset_builder(
            assessment, value, analysis_chembl, terminated_array, diseaseTA
        )


def comparisons_df() -> list:
    """Return list of all comparisons to be used in the analysis"""
    comparisons = spark.createDataFrame(
        data=[
            ("hasGeneticEvidence", "byDatatype"),
            ("diagonalYes", "byDatatype"),
            ("oneCellYes", "byDatatype"),
        ],
        schema=StructType(
            [
                StructField("comparison", StringType(), True),
                StructField("comparisonType", StringType(), True),
            ]
        ),
    )

    predictions = spark.createDataFrame(
        data=[
            ("Phase4", "clinical"),
            ("Phase>=3", "clinical"),
            ("Phase>=2", "clinical"),
            ("Phase>=1", "clinical"),
            ("PhaseT", "clinical"),
            ("approved", "clinical"),
        ]
    )
    return comparisons.join(predictions, how="full").collect()


result = []
result_st = []
result_ci = []
array2 = []
results = []


def aggregations_original(
    df,
    data,
    listado,
    comparisonColumn,
    comparisonType,
    predictionColumn,
    predictionType,
    today_date,
):

    wComparison = Window.partitionBy(comparisonColumn)
    wPrediction = Window.partitionBy(predictionColumn)
    wPredictionComparison = Window.partitionBy(comparisonColumn, predictionColumn)

    uniqIds = df.select("targetId", "diseaseId").distinct().count()

    out = (
        df.withColumn("comparisonType", F.lit(comparisonType))
        .withColumn("predictionType", F.lit(predictionType))
        .withColumn("total", F.lit(uniqIds))
        .withColumn("a", F.count("targetId").over(wPredictionComparison))
        .withColumn(
            "predictionTotal",
            F.count("targetId").over(wPrediction),
        )
        .withColumn(
            "comparisonTotal",
            F.count("targetId").over(wComparison),
        )
        .select(
            F.col(predictionColumn).alias("prediction"),
            F.col(comparisonColumn).alias("comparison"),
            "comparisonType",
            "predictionType",
            "a",
            "predictionTotal",
            "comparisonTotal",
            "total",
        )
        .filter(F.col("prediction").isNotNull())
        .filter(F.col("comparison").isNotNull())
        .distinct()
    )

    out.write.mode("overwrite").parquet(
        "gs://ot-team/jroldan/"
        + str(
            today_date
            + "_"
            + "analysis/"
            + data
            # + "_propagated"
            + "/"
            + comparisonColumn
            + "_"
            + predictionColumn
            + ".parquet"
        )
    )

    filePath = "gs://ot-team/jroldan/" + str(
        today_date
        + "_"
        + "analysis/"
        + data
        # + "_propagated"
        + "/"
        + comparisonColumn
        + "_"
        + predictionColumn
        + ".parquet"
    )
    listado.append(
        "gs://ot-team/jroldan/"
        + str(
            today_date
            + "_"
            + "analysis/"
            + data
            # + "_propagated"
            + "/"
            + comparisonColumn
            + "_"
            + predictionColumn
            + ".parquet"
        )
    )
    print(
        today_date
        + "_"
        + "analysis/"
        + data
        # + "_propagated"
        + "/"
        + comparisonColumn
        + "_"
        + predictionColumn
        + ".parquet"
    )
    c = datetime.now()
    c.strftime("%H:%M:%S")
    print(c)

    array1 = np.delete(
        out.join(full_data, on=["prediction", "comparison"], how="outer")
        .groupBy("comparison")
        .pivot("prediction")
        .agg(F.first("a"))
        .sort(F.col("comparison").desc())
        .select("comparison", "yes", "no")
        .fillna(0)
        .toPandas()
        .to_numpy(),
        [0],
        1,
    )
    total = np.sum(array1)
    res_npPhaseX = np.array(array1, dtype=int)
    resX = convertTuple(fisher_exact(res_npPhaseX, alternative="two-sided"))
    resx_CI = convertTuple(
        odds_ratio(res_npPhaseX).confidence_interval(confidence_level=0.95)
    )

    result_st.append(resX)
    result_ci.append(resx_CI)
    (rs_result, rs_ci) = relative_success(array1)

    results.append(
        [
            data,
            comparisonColumn,
            predictionColumn,
            round(float(resX.split(",")[0]), 2),
            float(resX.split(",")[1]),
            round(float(resx_CI.split(",")[0]), 2),
            round(float(resx_CI.split(",")[1]), 2),
            str(total),
            np.array(res_npPhaseX).tolist(),
            round(float(rs_result), 2),
            round(float(rs_ci[0]), 2),
            round(float(rs_ci[1]), 2),
            filePath,
        ]
    )
    return results


c = datetime.now()

print("start doing aggregations and writing")
today_date = str(date.today())
aggSetups_original = comparisons_df()
listado = []
results = []

print("starting with non-propagated aggregations at", c)

for key, df in dfs_dict.items():
    df = df.persist()
    for row in aggSetups_original:
        aggregations_original(df, key, listado, *row, today_date)
    df.unpersist()
    print(key + " df unpersisted")

print("non propagated files wroten succesfully at", c)

print("starting with non-propagated aggregations at", c)
for key, df in dfs_dict_propag.items():
    df = df.persist()
    for row in aggSetups_original:
        aggregations_original(df, key, listado, *row, today_date)
    df.unpersist()
    print(key + " df unpersisted")

print("non propagated files wroten succesfully at", c)


print("propagated files wroten succesfully at", c)
print("creating pandas dataframe with resulting rows")
df_results = pd.DataFrame(
    results,
    columns=[
        "group",
        "comparison",
        "phase",
        "OR",
        "pValue",
        "LowCI",
        "HighCI",
        "total",
        "array",
        "rs",
        "lowRs",
        "HighRs",
        "path",
    ],
)
print("created pandas dataframe")
print("converting to spark dataframe")
print("preparing dataframe")

schema = StructType(
    [
        StructField("group", StringType(), True),
        StructField("comparison", StringType(), True),
        StructField("phase", StringType(), True),
        StructField("oddsRatio", DoubleType(), True),
        StructField("pValue", DoubleType(), True),
        StructField("lowerInterval", DoubleType(), True),
        StructField("upperInterval", DoubleType(), True),
        StructField("total", StringType(), True),
        StructField("values", ArrayType(ArrayType(IntegerType())), True),
        StructField("relSuccess", DoubleType(), True),
        StructField("rsLower", DoubleType(), True),
        StructField("rsUpper", DoubleType(), True),
        StructField("path", StringType(), True),
    ]
)

print("read pattern variables")
df = spreadSheetFormatter(spark.createDataFrame(df_results, schema=schema))
print("processed spreadsheet")
print("writting the dataframe")

# Convert list of lists to DataFrame
# Regular expressions
value_pattern = r"df_([^_]+)_"  # Extracts {value}
middle_pattern = r"df_[^_]+_([^_]+)_"  # Extracts middle part (All, Other, etc.)
suffix_pattern = r"(original|propag)$"  # Extracts suffix (original or propag)

df.withColumn("datasource", regexp_extract("group", value_pattern, 1)).withColumn(
    "therArea", regexp_extract("group", middle_pattern, 1)
).withColumn("criteria", regexp_extract("group", suffix_pattern, 1)).toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_genEvidAnalysis.csv"
)

print("dataframe written \n Analysis finished")

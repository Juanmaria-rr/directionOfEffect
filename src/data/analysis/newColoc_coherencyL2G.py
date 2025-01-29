#### BUILDING THE NEW GWAS GENETIC EVIDENCE FROM COLOC
from functions import discrepancifier
from DoEAssessment import directionOfEffect
from functions import relative_success
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
import matplotlib.pyplot as plt
from datetime import date, datetime
from pyspark.sql.types import (
    StructType,
    StructField,
    ArrayType,
    DoubleType,
    DecimalType,
    StringType,
    FloatType,
    IntegerType,
)
import pandas as pd

spark = SparkSession.builder.getOrCreate()

path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"

#### Now load sources of data to generate credible_set_OT_genetics evidences and associations.

target = spark.read.parquet(f"{path}targets/")

diseases = spark.read.parquet(f"{path}diseases/")

evidences = spark.read.parquet(f"{path}evidence").filter(
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
ot_genetics = evidences.filter(F.col("datasourceId") == "ot_genetics_portal")

credibleEvidence = spark.read.parquet(f"{path}evidence").filter(
    F.col("datasourceId").isin(["gwas_credible_sets"])
)
credible = spark.read.parquet(f"{path}credibleSet")

index = spark.read.parquet(f"{path}gwasIndex")

new = spark.read.parquet(f"{path}colocalisation/coloc")

variantIndex = spark.read.parquet(f"{path}variantIndex")

biosample = spark.read.parquet(f"{path}biosample")

print("read spark files")

print("fixing scXQTL and XQTL studies")
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

print("fixed scXQTL and XQTL studies")

print("creating new coloc")

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
df = credibleEvidence.filter((F.col("datasourceId") == "gwas_credible_sets"))

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
)
print("creating new gwasResolvedColoc")

### bring directionality from QTL

gwasResolvedColoc = (
    (
        newColoc.filter(F.col("rightStudyType") != "gwas")
        .withColumnRenamed("geneId", "targetId")
        .join(
            gwasComplete.withColumnRenamed("studyLocusId", "leftStudyLocusId"),
            on=["leftStudyLocusId", "targetId"],
            how="right",
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

#### take the direction from the lowest p value
window_spec = Window.partitionBy("targetId", "diseaseId").orderBy(
    F.col("pValueExponent").asc()
)

print("creating new gwasCredibleAssoc")

### modify to include more information
gwasCredibleAssoc = (
    gwasResolvedColoc.withColumn(
        "homogenized", F.first("colocDoE", ignorenulls=True).over(window_spec)
    )
    .select(
        "targetId",
        "diseaseId",
        "homogenized",
        "leftStudyLocusId",
        "h4",
        "datasourceId",
        "resourceScore",
        "leftVariantId",
        "credibleLeftStudyType",
    )
    .withColumn(
        "homogenized",
        F.when(F.col("homogenized").isNull(), F.lit("noEvaluable")).otherwise(
            F.col("homogenized")
        ),
    )
)  ### there will be duplicates TargetId-DiseaseId because we are taking the most significant DoE

#### LOAD STUDYLOCUSID AND VARIANT DISTANCES
l2gPred = spark.read.parquet(
    "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/locusToGenePredictions"
)
l2gTable = (
    l2gPred.select("studyLocusId", "geneId", F.explode_outer("locusToGeneFeatures"))
    .filter(F.col("key").isin(["distanceFootprintMean", "distanceTssMean"]))
    .groupBy("studyLocusId", "geneId")
    .pivot("key")
    .agg(F.first("value"))
)
print("creating gwasCredibleAssocDistances")
gwasCredibleAssocDistances = gwasCredibleAssoc.join(
    l2gTable.withColumnRenamed("studyLocusId", "leftStudyLocusId").withColumnRenamed(
        "geneId", "targetId"
    ),
    on=["leftStudyLocusId", "targetId"],
    how="left",
)

print("creating analysis_chembl")
analysis_chembl = (
    discrepancifier(
        directionOfEffect(
            evidences.filter((F.col("datasourceId") == "chembl")), "24.09"
        )
        .withColumn(
            "maxClinPhase",
            F.max(F.col("clinicalPhase")).over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .groupBy("targetId", "diseaseId", "maxClinPhase")
        .pivot("homogenized")
        .agg(F.count("targetId"))
        .persist()
    )
    .filter(  ### ensure drug has annotated MoA and is coherent per Target-Disease
        ((F.col("GoF_protect").isNotNull()) | F.col("LoF_protect").isNotNull())
        & (F.col("coherencyDiagonal") == "coherent")
    )
    .selectExpr(
        "targetId",
        "diseaseId",
        "maxClinPhase",
        "coherencyDiagonal as coherencyDiagonal_ch",
        "coherencyOneCell as coherencyOneCell_ch",
        "LoF_protect as LoF_protect_ch",
        "GoF_protect as GoF_protect_ch",
    )
)

### pivot colocdoE grouping by T-D-studyLocusId-distances

value = "max_L2GScore"

list_l2g = [
    0.10,
    0.15,
    0.20,
    0.25,
    0.30,
    0.35,
    0.40,
    0.45,
    0.50,
    0.55,
    0.60,
    0.65,
    0.70,
    0.75,
    0.80,
    0.85,
    0.90,
    0.95,
]
print("creating benchmarkOT function")

dict_comb = {}


def benchmarkOT(
    dict_comb, value, gwasCredibleAssocDistances, analysis_chembl, list_l2g
):
    dict_comb = {
        "hasGeneticEvidence": f"{value}",
        "diagonalYes": f"{value}",
        "oneCellYes": f"{value}",
        "L2GAndColoc": f"{value}",
    }
    min_value = gwasCredibleAssocDistances.agg(
        F.min("distanceFootprintMean")
    ).collect()[0][
        0
    ]  ## take min value for scalating

    return (
        discrepancifier(
            gwasCredibleAssocDistances
            # .filter(F.col("h4").isNotNull()) #### not filter by this because we want to include the L2G AND Coloc question
            .withColumn(  ### take maximum L2G score per T-D
                "max_L2GScore",
                F.max("resourceScore").over(
                    Window.partitionBy("targetId", "diseaseId")
                ),
            )
            .withColumn("min_distanceFootprintMean", F.lit(min_value))
            .withColumn(
                "scaled_distanceFootprintMean",
                (F.col("distanceFootprintMean") - F.col("min_distanceFootprintMean"))
                / (1 - F.col("min_distanceFootprintMean")),
            )
            .withColumn(
                "minDistFootPrintMean",
                F.min("scaled_distanceFootprintMean").over(
                    Window.partitionBy("targetId", "diseaseId")
                ),
            )
            .withColumn(
                "minTssDistance",
                F.min("distanceTssMean").over(
                    Window.partitionBy("targetId", "diseaseId")
                ),
            )
            .groupBy(
                "targetId",
                "diseaseId",
                f"{value}",
                # "leftStudyLocusId",
            )
            .pivot("homogenized")
            .count()
        )
        .join(analysis_chembl, on=["targetId", "diseaseId"], how="right")
        .withColumn(
            "diagonalAgreeWithDrugs",
            F.when(
                (
                    (F.col("coherencyDiagonal_ch") == "coherent")
                    & (F.col("coherencyDiagonal") == "coherent")
                )
                # & (F.col("coherencyDiagonal") == "coherent")
                ,
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
                (
                    (F.col("coherencyOneCell_ch") == "coherent")
                    & (F.col("coherencyDiagonal") == "coherent")
                )
                # & (F.col("coherencyOneCell") == "coherent")
                ,
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
            # ).filter(
            #    F.col("diagonalAgreeWithDrugs").isNotNull()
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
            "diagonalYes",
            F.when(
                F.col("diagonalAgreeWithDrugs") == "coherent", F.lit("yes")
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "oneCellYes",
            F.when(
                F.col("oneCellAgreeWithDrugs") == "coherent", F.lit("yes")
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "hasGeneticEvidence",
            F.when(F.col(f"{value}").isNotNull(), F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "L2GAndColoc",
            F.when(
                (F.col(f"{value}").isNotNull())
                & (F.col("coherencyDiagonal").isin(["coherent", "dispar"])),
                F.lit("yes"),
            ).otherwise(F.lit("no")),
        )
        .select(
            ["*"]
            + (
                [  ### single columns
                    F.when(F.col(f"{value}") >= n, F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{value}>={str(n).replace('.', '_')}")
                    for n in list_l2g
                ]
            )
            + (
                [  ### column combinations for Yes/No colums Plus has DoE (any agreement)
                    F.when((F.col(a) == "yes") & (F.col(x) >= n), F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{x}>={str(n).replace('.', '_')}&{a}_combined")
                    for a, x in dict_comb.items()
                    for n in list_l2g
                ]
            )
        )
        .persist()
    )


print("creating dataframes in loop")

values = [
    "max_L2GScore",
    "minDistFootPrintMean",
    "minTssDistance",
]

datasetDict = {}
for value in values:
    if value == "max_L2GScore":
        datasetDict[f"df_l2g_original"] = benchmarkOT(
            dict_comb, value, gwasCredibleAssocDistances, analysis_chembl, list_l2g
        )
    else:
        datasetDict[f"{value}"] = benchmarkOT(
            dict_comb, value, gwasCredibleAssocDistances, analysis_chembl, list_l2g
        )


def comparisons_df(dataset) -> list:
    """Return list of all comparisons to be used in the analysis"""
    toAnalysis = dataset.columns[22:]
    dataType = ["byDatatype"] * len(toAnalysis)
    l_studies = []
    l_studies.extend([list(a) for a in zip(toAnalysis, dataType)])

    schema = StructType(
        [
            StructField("comparison", StringType(), True),
            StructField("comparisonType", StringType(), True),
        ]
    )

    comparisons = spark.createDataFrame(l_studies, schema=schema)
    ### include all the columns as predictor

    predictions = spark.createDataFrame(
        data=[
            ("Phase4", "clinical"),
            # ("Phase>=3", "clinical"),
            # ("Phase>=2", "clinical"),
            # ("Phase>=1", "clinical"),
            # ("PhaseT", "clinical"),
        ]
    )
    return comparisons.join(predictions, how="full").collect()


from functions import relative_success, spreadSheetFormatter, convertTuple

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

    wComparison = Window.partitionBy(F.col(comparisonColumn))
    wPrediction = Window.partitionBy(F.col(predictionColumn))
    wPredictionComparison = Window.partitionBy(
        F.col(comparisonColumn), F.col(predictionColumn)
    )

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
    """
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
    """
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


import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

# Initialize an empty list to store the results
result_st = []
result_ci = []
results = []


def convertTuple(tup):
    st = ",".join(map(str, tup))
    return st


print("launched function to run analysis")


listado = []
today_date = str(date.today())
for key, df_analysis in datasetDict.items():
    aggSetups_original = comparisons_df(df_analysis)
    print("corresponding dataframe key: ", key)
    df_analysis.persist()
    for row in aggSetups_original:
        print(key, value)
        aggregations_original(df_analysis, key, listado, *row, today_date)

print("finished analysis")
print("creating pandas dataframe with resulting rows")
df = pd.DataFrame(
    results,
    columns=[
        "type",
        "criteria",
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

schema = StructType(
    [
        StructField("type", StringType(), True),
        StructField("criteria", StringType(), True),
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

# Convert list of lists to DataFrame
print("created pandas dataframe")
print("converting to spark dataframe")
pattern = r"(\b[0-9]+)_(\d+\b)"  # Matches patterns like 0_1, 0_95, etc.
df = spreadSheetFormatter(spark.createDataFrame(results, schema=schema)).withColumn(
    "range", F.regexp_extract(F.col("criteria"), pattern, 0)
)
print("writting dataframe")
df.toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_newColoc_L2Gbenchmark.csv"
)

print("dataframe written \n Analysis finished")

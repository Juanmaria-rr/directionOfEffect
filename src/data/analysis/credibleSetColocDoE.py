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

spark = SparkSession.builder.getOrCreate()

path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"

target = spark.read.parquet(f"{path}targets/")

diseases = spark.read.parquet(f"{path}diseases/")

evidences = spark.read.parquet(f"{path}evidence")

credible = spark.read.parquet(f"{path}credibleSet")

index = spark.read.parquet(f"{path}gwasIndex")

new = spark.read.parquet(f"{path}colocalisation/coloc")

variantIndex = spark.read.parquet(f"{path}variantIndex")

biosample = spark.read.parquet(f"{path}biosample")

newColoc = (
    new.join(
        credible.selectExpr(  #### studyLocusId from credible set to uncover the codified variants on left side
            "studyLocusId as leftStudyLocusId",
            "StudyId as leftStudyId",
            "variantId as leftVariantId",
            "chromosome",
            # "region",
            # "beta",
            # "pValueExponent",
            "studyType",
        ).withColumnRenamed(
            "studyType", "credibleLeftStudyType"
        ),
        on="leftStudyLocusId",
        how="left",
    )
    .join(
        credible.selectExpr(  #### studyLocusId from credible set to uncover the codified variants on right side
            "studyLocusId",
            "studyId",
            "variantId as rightVariantId",
            # "chromosome",
            # "region",
            # "beta",
            # "pValueExponent",
            "studyType",
        )
        .withColumnRenamed("studyLocusId", "rightStudyLocusId")
        .withColumnRenamed("studyId", "rightStudyId")
        .withColumnRenamed("beta", "rightBeta")
        .withColumnRenamed("studyType", "credibleRightStudyType"),
        on="rightStudyLocusId",
        how="left",
    )
    .join(
        index.select(  ### bring modulated target on right side (QTL study)
            "studyId", "geneId", "projectId", "studyType", "condition", "biosampleId"
        )
        .withColumnRenamed("studyId", "rightStudyId")
        .withColumnRenamed("studyType", "indexStudyType"),
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

chemblAssoc = (
    discrepancifier(
        assessment.filter(
            (F.col("datasourceId") == "chembl")
            & (F.col("homogenized") != "noEvaluable")
        )
        .withColumn(
            "maxClinPhase",
            F.max("clinicalPhase").over(Window.partitionBy("targetId", "diseaseId")),
        )
        .groupBy("targetId", "diseaseId", "maxClinPhase")
        .pivot("homogenized")
        .count()
    )
    .filter(F.col("coherencyDiagonal") == "coherent")
    .drop(
        "coherencyDiagonal", "coherencyOneCell", "noEvaluable", "GoF_risk", "LoF_risk"
    )
    .withColumnRenamed("GoF_protect", "drugGoF_protect")
    .withColumnRenamed("LoF_protect", "drugLoF_protect")
)

benchmark = (
    (
        resolvedColoc.filter(F.col("betaGwas") < 0)
        .join(  ### select just GWAS giving protection
            chemblAssoc, on=["targetId", "diseaseId"], how="inner"
        )
        .withColumn(
            "AgreeDrug",
            F.when(
                (F.col("drugGoF_protect").isNotNull())
                & (F.col("colocDoE") == "GoF_protect"),
                F.lit("coherent"),
            )
            .when(
                (F.col("drugLoF_protect").isNotNull())
                & (F.col("colocDoE") == "LoF_protect"),
                F.lit("coherent"),
            )
            .otherwise(F.lit("dispar")),
        )
    )
    .filter(F.col("name") != "COVID-19")  #### remove COVID-19 associations
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
).join(biosample.select("biosampleId", "biosampleName"), on="biosampleId", how="left")

#### one-hot encoding column values yes/no
columns = ["colocDoE", "biosampleName", "rightStudyType", "projectId"]

# Build the distinct value-to-column mapping
distinct_values = (
    benchmark.select(
        F.explode(
            F.array(
                *[
                    F.struct(F.col(c).alias("value"), F.lit(c).alias("column"))
                    for c in columns
                ]
            )
        ).alias("exploded")
    )
    .select("exploded.value", "exploded.column")
    .distinct()
)

# Convert to dictionary: {value: column_name}
value_to_column_mapping = {
    row["value"]: row["column"] for row in distinct_values.collect()
}

# One-hot encoding
for col_name in columns:
    # Get distinct values for the column
    distinct_col_values = [
        row[0] for row in benchmark.select(col_name).distinct().collect()
    ]

    # Add a new column for each distinct value
    for val in distinct_col_values:
        new_col_name = f"{val}"
        benchmark = benchmark.withColumn(
            new_col_name,
            F.when(F.col(col_name) == val, F.lit("yes")).otherwise(F.lit("no")),
        )

### prepare list to feed the analysis


def comparisons_df(value_to_column_mapping) -> list:
    """Return list of all comparisons to be used in the analysis"""
    ## define the list to have all the columns no analyse and their comparison identification
    toAnalysis = []
    ### include key and values from columns to analysed done by dictionary
    toAnalysis = [[key, x] for key, x in value_to_column_mapping.items()]
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
print("created full_data and lists")

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from pyspark.sql.types import *


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
    results = []
    # uniqIds = df.select("targetId", "diseaseId").distinct().count()
    out = (
        df.withColumn("comparisonType", F.lit(comparisonType))
        .withColumn(
            "study",
            F.when(
                F.col(comparisonColumn) == "yes",
                F.collect_set(F.col("rightStudyId")).over(
                    Window.partitionBy(comparisonColumn)
                ),
            ),
        )
        .withColumn(
            "tissue",
            F.when(
                F.col(comparisonColumn) == "yes",
                F.collect_set(F.col("biosampleName")).over(
                    Window.partitionBy(comparisonColumn)
                ),
            ),
        )
        .withColumn("dataset", F.lit(data))
        .withColumn("predictionType", F.lit(predictionType))
        # .withColumn("total", F.lit(uniqIds))
        .withColumn("a", F.count("targetId").over(wPredictionComparison))
        .withColumn("comparisonColumn", F.lit(comparisonColumn))
        .withColumn("predictionColumnValue", F.lit(predictionColumn))
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
            "dataset",
            "comparisonColumn",
            "predictionColumnValue",
            "comparisonType",
            "predictionType",
            "a",
            "predictionTotal",
            "comparisonTotal",
            # "total",
            "study",
            "tissue",
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
            + comparisonType
            + "_"
            + predictionColumn
            + ".parquet"
        )
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
            + comparisonType
            + "_"
            + predictionColumn
            + ".parquet"
        )
    )
    path = "gs://ot-team/jroldan/" + str(
        today_date
        + "_"
        + "analysis/"
        + data
        # + "_propagated"
        + "/"
        + comparisonColumn
        + "_"
        + comparisonType
        + "_"
        + predictionColumn
        + ".parquet"
    )
    print(path)
    ### making analysis
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
    '''
    studies = (
        out.filter(F.col("comparison") == "yes")
        .select(F.flatten(F.collect_set(F.col("study"))).alias("study"))
        .collect()[0]["study"]
    )
    tissues = (
        out.filter(F.col("comparison") == "yes")
        .select(F.flatten(F.collect_set(F.col("tissue"))).alias("tissue"))
        .collect()[0]["tissue"]
    )
    '''
    results.extend(
        [
            comparisonType,
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
            # studies,
            # tissues,
            path,
        ]
    )
    return results


result = []
result_st = []
result_ci = []
array2 = []

# Initialize an empty list to store the results
results = []


def convertTuple(tup):
    st = ",".join(map(str, tup))
    return st


aggSetups_original = comparisons_df(value_to_column_mapping)
listado = []
result_all = []
today_date = str(date.today())

for row in aggSetups_original:
    results = aggregations_original(benchmark, "propagated", listado, *row, today_date)
    result_all.append(results)

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
        #StructField("studies", StringType(), True),
        #StructField("tissues", StringType(), True),
        StructField("path", StringType(), True),
    ]
)

# Convert list of lists to DataFrame
df = spreadSheetFormatter(spark.createDataFrame(result_all, schema=schema))
df.toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_credibleSetColocDoEanalysis.csv"
)

print("dataframe written \n Analysis finished")

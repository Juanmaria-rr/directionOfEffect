import time
from array import ArrayType
from functions import (
    relative_success,
    spreadSheetFormatter,
    discrepancifier,
    temporary_directionOfEffect,
)
# from stoppedTrials import terminated_td
from DoEAssessment import directionOfEffect
# from membraneTargets import target_membrane
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
spark.conf.set(
    "spark.sql.shuffle.partitions", "400"
)  # Default is 200, increase if needed


path_n='gs://open-targets-data-releases/25.03/output/'

target = spark.read.parquet(f"{path_n}target/")

diseases = spark.read.parquet(f"{path_n}disease/")

evidences = spark.read.parquet(f"{path_n}evidence")

credible = spark.read.parquet(f"{path_n}credible_set")

new = spark.read.parquet(f"{path_n}colocalisation_coloc") 

index=spark.read.parquet(f"{path_n}study/")

variantIndex = spark.read.parquet(f"{path_n}variant")

biosample = spark.read.parquet(f"{path_n}biosample")

print("loaded files")

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
            'isTransQtl'
        ),
        on="rightStudyLocusId",
        how="left",
    )
    .join(
        index.selectExpr(  ### bring modulated target on right side (QTL study)
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
    # .persist()
)

print("loaded newColoc")

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
filtered_df = df.select(*non_null_columns)  # .persist()

## bring studyId, variantId, beta from Gwas and pValue
gwasComplete = filtered_df.join(
    credible.selectExpr(
        "studyLocusId", "studyId", "variantId", "beta as betaGwas", "pValueExponent"
    ),
    on="studyLocusId",
    how="left",
)  # .persist()

print("loaded gwasComplete")

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
    ).withColumn(
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
    # .persist()
)
print("loaded resolvedColloc")

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
    path_n, datasource_filter
)

print("run temporary direction of effect")


print("built drugApproved dataset")

analysis_chembl_indication = (
    discrepancifier(
        assessment.filter((F.col("datasourceId") == "chembl"))
        .withColumn(
            "maxClinPhase",
            F.max(F.col("clinicalPhase")).over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .groupBy("targetId", "diseaseId", "maxClinPhase")
        .pivot("homogenized")
        .agg(F.count("targetId"))
    )
    .filter(F.col("coherencyDiagonal") == "coherent")
    .drop(
        "coherencyDiagonal", "coherencyOneCell", "noEvaluable", "GoF_risk", "LoF_risk"
    )
    .withColumnRenamed("GoF_protect", "drugGoF_protect")
    .withColumnRenamed("LoF_protect", "drugLoF_protect")
    # .persist()
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

print("built chemblAssoc dataset")

benchmark = (
    (
        resolvedColoc.filter(F.col("betaGwas") < 0)
        .join(  ### select just GWAS giving protection
            analysis_chembl_indication, on=["targetId", "diseaseId"], how="inner"
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


print("built benchmark dataset")

#### Analysis

#### 1 Build a dictionary with the distinct values as key and column names as value
variables_study = ["projectId", "biosampleName", "rightStudyType", "colocDoE"]

# Initialize an empty dictionary
disdic = {}

# Iterate over the list of column names
for col_name in variables_study:
    # Extract distinct values for the column
    distinct_values = benchmark.select(col_name).distinct().collect()

    # Populate the dictionary
    for row in distinct_values:
        distinct_value = row[col_name]
        if distinct_value is not None:  # Exclude None (null) values
            disdic[distinct_value] = col_name

####2 Define agregation function
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from pyspark.sql.types import *


def convertTuple(tup):
    st = ",".join(map(str, tup))
    return st


#####3 run in a function
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


#### 3 Loop over different datasets (as they will have different rows and columns)


def comparisons_df_iterative(elements):
    # toAnalysis = [(key, value) for key, value in disdic.items() if value == projectId]
    toAnalysis = [(col, "predictor") for col in elements]
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
            ('Phase>=3','clinical'),
            ('Phase>=2','clinical'),
            ('Phase>=1','clinical'),
            ("PhaseT", "clinical"),
        ]
    )
    return comparisons.join(predictions, how="full").collect()


print("load comparisons_df_iterative function")


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

rightTissue = spark.read.csv(
    'gs://ot-team/jroldan/analysis/20250526_curatedRightTissue.csv',
    header=True,
).drop("_c0")

print("loaded rightTissue dataset")

negativeTD = (
    evidences.filter(F.col("datasourceId") == "chembl")
    .select("targetId", "diseaseId", "studyStopReason", "studyStopReasonCategories")
    .filter(F.array_contains(F.col("studyStopReasonCategories"), "Negative"))
    .groupBy("targetId", "diseaseId")
    .count()
    .withColumn("stopReason", F.lit("Negative"))
    .drop("count")
)

print("built negativeTD dataset")

bench2 = benchmark.join(
    rightTissue, on=["name", "bioSampleName"], how="left"
).withColumn(
    "rightTissue",
    F.when(F.col("rightTissue1") == "yes", F.lit("yes")).otherwise(F.lit("no")),
)

print("built bench2 dataset")

###### cut from here
print("looping for variables_study")
# List of columns to analyze
variables_study = ["projectId", "biosampleName", "rightStudyType", "colocDoE"]

# Dictionary to store results
pivoted_dfs = {}

# Loop over the columns
for col in variables_study:
    window_spec = Window.partitionBy("targetId", "diseaseId", col).orderBy(
        F.col("pValueExponent").asc()
    )
    print(f"Processing: {col}")

    pivoted_df = (
        bench2.withColumn(
            "rightTissue",
            F.when(F.col("rightTissue1") == "yes", F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "agree_lowestPval",
            F.first("AgreeDrug", ignorenulls=True).over(
                window_spec
            ),  ### ignore nulls aded 29.01.2025
            #### take directionality from lowest p value
        )
        .withColumn(
            "isRightTissueSignalAgreed",
            F.collect_set(
                F.when(F.col("rightTissue") == "yes", F.col("agree_lowestPval"))
            ).over(window_spec),
        )
        .withColumn(
            "isSignalFromRightTissue",
            F.first(
                F.when(
                    F.col("AgreeDrug") == F.col("agree_lowestPval"),
                    F.col("rightTissue"),
                ),
                ignorenulls=True,
            ).over(window_spec),
        )
        .groupBy(
            "targetId",
            "diseaseId",
            "maxClinPhase",
            "rightTissue",
            "isRightTissueSignalAgreed",
            "isSignalFromRightTissue",
        )
        .pivot(col)  # Pivot the column dynamically
        .agg(F.collect_set("agree_lowestPval"))
        .join(negativeTD, on=["targetId", "diseaseId"], how="left")
        .withColumn(
            "PhaseT",
            F.when(F.col("stopReason") == "Negative", F.lit("yes")).otherwise(
                F.lit("no")
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
        .select(
            ["*"]
            + (
                [  ### single columns
                    F.when(F.array_contains(F.col(x), "yes"), F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{x}_only")
                    for x, value in [
                        (key, val) for key, val in disdic.items() if val == col
                    ]
                ]
            )
            + (
                [
                    F.when(
                        F.array_contains(F.col("isRightTissueSignalAgreed"), "yes"),
                        F.lit("yes"),
                    )
                    .otherwise(F.lit("no"))
                    .alias(f"{x}_isRightTissueSignalAgreed")
                    for x, value in [
                        (key, val) for key, val in disdic.items() if val == col
                    ]
                ]
            )
        )
    )  # Collect unique values

    # Store the DataFrame in the dictionary
    pivoted_dfs[col] = pivoted_df



result = []
result_st = []
result_ci = []
array2 = []
listado = []
result_all = []
today_date = str(date.today())
variables_study = ["projectId", "biosampleName", "rightStudyType", "colocDoE"]

##### PROJECT ID ###### 
print('working with projectId')
pivoted_dfs['projectId'].persist()
unique_values = benchmark.select('projectId').distinct().rdd.flatMap(lambda x: x).collect()
filter = len(pivoted_dfs['projectId'].drop(*unique_values).columns[12:])
print('There are ', filter, 'columns to analyse with phases')
rows = comparisons_df_iterative(pivoted_dfs['projectId'].columns[-filter:])

# If needed, now process the rest
for row in rows:
    results = aggregations_original(
        pivoted_dfs['projectId'], "propagated", listado, *row, today_date
    )
    result_all.append(results)

pivoted_dfs['projectId'].unpersist()
print('df unpersisted')

##### BIOSAMPLE NAME ###### 
print('working with biosampleName')
pivoted_dfs['biosampleName'].persist()
unique_values = benchmark.select('biosampleName').distinct().rdd.flatMap(lambda x: x).collect()
filter = len(pivoted_dfs['biosampleName'].drop(*unique_values).columns[12:])
print('There are ', filter, 'columns to analyse with phases')
rows = comparisons_df_iterative(pivoted_dfs['biosampleName'].columns[-filter:])

for row in rows:
    results = aggregations_original(
        pivoted_dfs['biosampleName'], "propagated", listado, *row, today_date
    )
    result_all.append(results)

pivoted_dfs['biosampleName'].unpersist()
print('df unpersisted')

##### RIGHTSTUDYTYPE  ###### 
print('working with rightStudyType')
pivoted_dfs['rightStudyType'].persist()
unique_values = benchmark.select('rightStudyType').distinct().rdd.flatMap(lambda x: x).collect()
filter = len(pivoted_dfs['rightStudyType'].drop(*unique_values).columns[12:])
print('There are ', filter, 'columns to analyse with phases')
rows = comparisons_df_iterative(pivoted_dfs['rightStudyType'].columns[-filter:])

for row in rows:
    results = aggregations_original(
        pivoted_dfs['rightStudyType'], "propagated", listado, *row, today_date
    )
    result_all.append(results)
pivoted_dfs['rightStudyType'].unpersist()
print('df unpersisted')

##### COLOC DOE ######
print('working with colocDoE')
pivoted_dfs['colocDoE'].persist()
unique_values = benchmark.select('colocDoE').distinct().rdd.flatMap(lambda x: x).collect()
filter = len(pivoted_dfs['colocDoE'].drop(*unique_values).columns[12:])
print('There are ', filter, 'columns to analyse with phases')
rows = comparisons_df_iterative(pivoted_dfs['colocDoE'].columns[-filter:])

for row in rows:
    results = aggregations_original(
        pivoted_dfs['colocDoE'], "propagated", listado, *row, today_date
    )
    result_all.append(results)
pivoted_dfs['colocDoE'].unpersist()
print('df unpersisted')

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
import re

# Define the list of patterns to search for
patterns = [
    "_only",
    #"_tissue",
    #"_isSignalFromRightTissue",
    "_isRightTissueSignalAgreed",
]
# Create a regex pattern to match any of the substrings
regex_pattern = "(" + "|".join(map(re.escape, patterns)) + ")"

# Convert list of lists to DataFrame
df = (
    spreadSheetFormatter(spark.createDataFrame(result_all, schema=schema))
    .withColumn(
        "prefix",
        F.regexp_replace(
            F.col("comparison"), regex_pattern + ".*", ""
        ),  # Extract part before the pattern
    )
    .withColumn(
        "suffix",
        F.regexp_extract(
            F.col("comparison"), regex_pattern, 0
        ),  # Extract the pattern itself
    )
)

### annotate projectId, tissue, qtl type and doe type:

from pyspark.sql.functions import create_map
from itertools import chain

mapping_expr=create_map([F.lit(x) for x in chain(*disdic.items())])

df_annot=df.withColumn('annotation',mapping_expr.getItem(F.col('prefix')))

df_annot.toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_credibleSetColocDoEanalysis_RightTissues.csv"
)

print("dataframe written \n Analysis finished")
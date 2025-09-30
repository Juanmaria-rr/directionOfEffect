import time
from array import ArrayType
from functions import (
    relative_success,
    spreadSheetFormatter,
    discrepancifier,
    temporary_directionOfEffect,
    buildColocData,
    gwasDataset,
    build_resolved_coloc,
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
from functools import reduce
# --- Build the SparkSession ---
# Use the .config() method to set these parameters before calling .getOrCreate()
# This ensures Spark requests the correct resources from YARN at the start.
driver_memory = "24g"                 # plenty for planning & small collects
executor_cores = 4                    # sweet spot for GC + Python workers
num_executors  = 12                   # 12 * 4 = 48 cores for executors; ~16 cores left for driver/OS
executor_memory = "32g"               # per executor heap
executor_memory_overhead = "8g"       # ~20% overhead for PySpark/Arrow/off-heap
# Totals: (32+8) * 12 = 480 GB executors + 24 GB driver ≈ 504 GB (adjust down if your hard cap is <500 GB)
# If you must stay strictly ≤ 500 GB, use executor_memory="30g", overhead="6g"  → (36 * 12) + 24 = 456 + 24 = 480 GB

shuffle_partitions   = 192            # ≈ 2–4× total cores (48) → start with 192
default_parallelism  = 192

spark = SparkSession.builder \
    .appName("MyOptimizedPySparkApp") \
    .config("spark.master", "yarn") \
    .config("spark.driver.memory", driver_memory) \
    .config("spark.executor.memory", executor_memory) \
    .config("spark.executor.cores", executor_cores) \
    .config("spark.executor.instances", num_executors) \
    .config("spark.yarn.executor.memoryOverhead", executor_memory_overhead) \
    .config("spark.sql.shuffle.partitions", shuffle_partitions) \
    .config("spark.default.parallelism", default_parallelism) \
    .getOrCreate()

print(f"SparkSession created successfully with the following configurations:")
print(f"  spark.driver.memory: {spark.conf.get('spark.driver.memory')}")
print(f"  spark.executor.memory: {spark.conf.get('spark.executor.memory')}")
print(f"  spark.executor.cores: {spark.conf.get('spark.executor.cores')}")
print(f"  spark.executor.instances: {spark.conf.get('spark.executor.instances')}")
print(f"  spark.yarn.executor.memoryOverhead: {spark.conf.get('spark.yarn.executor.memoryOverhead')}")
print(f"  spark.sql.shuffle.partitions: {spark.conf.get('spark.sql.shuffle.partitions')}")
print(f"  spark.default.parallelism: {spark.conf.get('spark.default.parallelism')}")
print(f"Spark UI available at: {spark.sparkContext.uiWebUrl}")

print(f"SparkSession created successfully with the following configurations:")
print(f"  spark.driver.memory: {spark.conf.get('spark.driver.memory')}")
print(f"  spark.executor.memory: {spark.conf.get('spark.executor.memory')}")
print(f"  spark.executor.cores: {spark.conf.get('spark.executor.cores')}")
print(f"  spark.executor.instances: {spark.conf.get('spark.executor.instances')}")
print(f"  spark.yarn.executor.memoryOverhead: {spark.conf.get('spark.yarn.executor.memoryOverhead')}")
print(f"  spark.sql.shuffle.partitions: {spark.conf.get('spark.sql.shuffle.partitions')}")
print(f"  spark.default.parallelism: {spark.conf.get('spark.default.parallelism')}")
print(f"Spark UI available at: {spark.sparkContext.uiWebUrl}")

# --- Your PySpark Code Here ---
# Now you can proceed with your data loading and processing.
# Example:
# df = spark.read.parquet("hdfs:///user/your_user/your_large_data.parquet")
# print(f"Number of rows in DataFrame: {df.count()}")
# df.groupBy("some_column").agg({"another_column": "sum"}).show()

# Remember to stop the SparkSession when you are done
# spark.stop()

path_n='gs://open-targets-data-releases/25.09/output/'

target = spark.read.parquet(f"{path_n}target/")

diseases = spark.read.parquet(f"{path_n}disease/")

evidences = spark.read.parquet(f"{path_n}evidence")

credible = spark.read.parquet(f"{path_n}credible_set")

new = spark.read.parquet(f"{path_n}colocalisation_coloc") 

index=spark.read.parquet(f"{path_n}study/")

variantIndex = spark.read.parquet(f"{path_n}variant")

biosample = spark.read.parquet(f"{path_n}biosample")

ecaviar=spark.read.parquet(f"{path_n}colocalisation_ecaviar")

all_coloc=ecaviar.unionByName(new, allowMissingColumns=True)

mecact_path = f"{path_n}drug_mechanism_of_action/" #  mechanismOfAction == old version


print("loaded files")

#### FIRST MODULE: BUILDING COLOC 
newColoc=buildColocData(all_coloc,credible,index)

print("loaded newColoc")

### SECOND MODULE: PROCESS EVIDENCES TO AVOID EXCESS OF COLUMNS 
gwasComplete = gwasDataset(evidences,credible)

print('gwasComplete loaded')
#### THIRD MODULE: INCLUDE COLOC IN THE 
resolvedColoc=build_resolved_coloc(newColoc, gwasComplete, diseases).withColumn('hasGenetics', F.lit('yes'))

datasource_filter = [
#   "ot_genetics_portal",
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


#### new part with chatgpt -- TEST

## QUESTIONS TO ANSWER:
# HAVE ECAVIAR >=0.8
# HAVE COLOC 
# HAVE COLOC >= 0.8
# HAVE COLOC + ECAVIAR >= 0.01
# HAVE COLOC >= 0.8 + ECAVIAR >= 0.01
# RIGHT JOING WITH CHEMBL 

### FIFTH MODULE: BUILDING BENCHMARK OF THE DATASET TO EXTRACT EHE ANALYSIS 

resolvedColocFiltered = resolvedColoc.filter((F.col('clpp')>=0.01) | (F.col('h4')>=0.8))

### drug mechanism of action
inhibitors = [
    "RNAI INHIBITOR",
    "NEGATIVE MODULATOR",
    "NEGATIVE ALLOSTERIC MODULATOR",
    "ANTAGONIST",
    "ANTISENSE INHIBITOR",
    "BLOCKER",
    "INHIBITOR",
    "DEGRADER",
    "INVERSE AGONIST",
    "ALLOSTERIC ANTAGONIST",
    "DISRUPTING AGENT",
]

activators = [
    "PARTIAL AGONIST",
    "ACTIVATOR",
    "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR",
    "AGONIST",
    "SEQUESTERING AGENT",  ## lost at 31.01.2025
    "STABILISER",
    # "EXOGENOUS GENE", ## added 24.06.2025
    # "EXOGENOUS PROTEIN" ## added 24.06.2025
]
mecact = spark.read.parquet(mecact_path)
actionType = (
        mecact.select(
            F.explode_outer("chemblIds").alias("drugId"),
            "actionType",
            "mechanismOfAction",
            "targets",
        )
        .select(
            F.explode_outer("targets").alias("targetId"),
            "drugId",
            "actionType",
            "mechanismOfAction",
        )
        .groupBy("targetId", "drugId")
        .agg(F.collect_set("actionType").alias("actionType2"))
    ).withColumn('nMoA', F.size(F.col('actionType2')))

analysis_chembl_indication = (
    discrepancifier(
        assessment.filter((F.col("datasourceId") == "chembl")).join(actionType, on=['targetId','drugId'], how='left')
        .withColumn(
            "maxClinPhase",
            F.max(F.col("clinicalPhase")).over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .groupBy("targetId", "diseaseId", "maxClinPhase",'actionType2')
        .pivot("homogenized")
        .agg(F.count("targetId"))
    )
    #.filter(F.col("coherencyDiagonal") == "coherent")
    .drop(
        "coherencyDiagonal", "coherencyOneCell", "noEvaluable", "GoF_risk", "LoF_risk"
    )
    .withColumnRenamed("GoF_protect", "drugGoF_protect")
    .withColumnRenamed("LoF_protect", "drugLoF_protect")
)

print("built drugApproved dataset")
benchmark = (
        resolvedColocFiltered.filter( ## .filter(F.col("betaGwas") < 0)
        F.col("name") != "COVID-19"
    )
        .join(  ### select just GWAS giving protection
            analysis_chembl_indication, on=["targetId", "diseaseId"], how="right"  ### RIGHT SIDE
        )
).join(biosample.select("biosampleId", "biosampleName"), on="biosampleId", how="left")

print("built benchmark")

#### HERE THE CODE FOR THE ANALYSIS

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
    #toAnalysis = [(key, value) for key, value in disdic.items() if value == projectId]
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
            ("Phase>=4", "clinical"),
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
    'gs://ot-team/jroldan/analysis/test_20250927_2509_biosampleName_curation_annotated.tsv',sep='\t',
    header=True,
).drop("Unnamed: 0")

print("loaded rightTissue dataset")

#### Do the join with Right Tissue annotation

benchmark_tissue=benchmark.join(rightTissue.drop('name'), on=['diseaseId','biosampleName'], how='left')

### create disdic dictionary
disdic={}

# --- Configuration for your iterative pivoting ---
group_by_columns = ['targetId', 'diseaseId','Relevant','phase4Clean','phase3Clean','phase2Clean','phase1Clean','PhaseT']
#columns_to_pivot_on = ['actionType2', 'biosampleName', 'projectId', 'rightStudyType','colocalisationMethod']
columns_to_pivot_on = ['biosampleName']
columns_to_aggregate = ['NoneCellYes', 'NdiagonalYes','hasGenetics'] # The values you want to collect in the pivoted cells
all_pivoted_dfs = {}

doe_columns=["LoF_protect", "GoF_risk", "LoF_risk", "GoF_protect"]
diagonal_lof=['LoF_protect','GoF_risk']
diagonal_gof=['LoF_risk','GoF_protect']

conditions = [
    F.when(F.col(c) == F.col("maxDoE"), F.lit(c)).otherwise(F.lit(None)) for c in doe_columns
    ]

# --- Nested Loops for Dynamic Pivoting ---


for pivot_col_name in columns_to_pivot_on:
    for agg_col_name in columns_to_aggregate:
        print(f"\n--- Creating DataFrame for Aggregation: '{agg_col_name}' and Pivot: '{pivot_col_name}' ---")
        current_col_pvalue_order_window = Window.partitionBy("targetId", "diseaseId", "maxClinPhase","Relevant", pivot_col_name).orderBy(F.col('colocalisationMethod').asc(), F.col("qtlPValueExponent").asc())
        test2=discrepancifier(benchmark_tissue.withColumn('actionType2', F.concat_ws(",", F.col("actionType2"))).withColumn('qtlColocDoE',F.first('colocDoE').over(current_col_pvalue_order_window)).groupBy(
        "targetId", "diseaseId", "hasGenetics","maxClinPhase", "drugLoF_protect", "drugGoF_protect","Relevant",pivot_col_name)
        .pivot("colocDoE")
        .count()
        .withColumnRenamed('drugLoF_protect', 'LoF_protect_ch')
        .withColumnRenamed('drugGoF_protect', 'GoF_protect_ch')).withColumn( ## .filter(F.col('coherencyDiagonal')!='noEvid')
    "arrayN", F.array(*[F.col(c) for c in doe_columns])
    ).withColumn(
        "maxDoE", F.array_max(F.col("arrayN"))
    ).withColumn("maxDoE_names", F.array(*conditions)
    ).withColumn("maxDoE_names", F.expr("filter(maxDoE_names, x -> x is not null)")
    ).withColumn(
        "NoneCellYes",
        F.when((F.col("LoF_protect_ch").isNotNull() & (F.col('GoF_protect_ch').isNull())) & (F.array_contains(F.col("maxDoE_names"), F.lit("LoF_protect")))==True, F.lit('yes'))
        .when((F.col("GoF_protect_ch").isNotNull() & (F.col('LoF_protect_ch').isNull())) & (F.array_contains(F.col("maxDoE_names"), F.lit("GoF_protect")))==True, F.lit('yes')
            ).otherwise(F.lit('no'))  # If the value is null, return null # Otherwise, check if name is in array
    ).withColumn(
        "NdiagonalYes",
        F.when((F.col("LoF_protect_ch").isNotNull() & (F.col('GoF_protect_ch').isNull())) & 
            (F.size(F.array_intersect(F.col("maxDoE_names"), F.array([F.lit(x) for x in diagonal_lof]))) > 0),
            F.lit("yes")
        ).when((F.col("GoF_protect_ch").isNotNull() & (F.col('LoF_protect_ch').isNull())) & 
            (F.size(F.array_intersect(F.col("maxDoE_names"), F.array([F.lit(x) for x in diagonal_gof]))) > 0),
            F.lit("yes")
        ).otherwise(F.lit('no'))
    ).withColumn(
        "drugCoherency",
        F.when(
            (F.col("LoF_protect_ch").isNotNull())
            & (F.col("GoF_protect_ch").isNull()), F.lit("coherent")
        )
        .when(
            (F.col("LoF_protect_ch").isNull())
            & (F.col("GoF_protect_ch").isNotNull()), F.lit("coherent")
        )
        .when(
            (F.col("LoF_protect_ch").isNotNull())
            & (F.col("GoF_protect_ch").isNotNull()), F.lit("dispar")
        )
        .otherwise(F.lit("other")),
    ).join(negativeTD, on=["targetId", "diseaseId"], how="left").withColumn(
        "PhaseT",
        F.when(F.col("stopReason") == "Negative", F.lit("yes")).otherwise(F.lit("no")),
    ).withColumn(
        "phase4Clean",
        F.when(
            (F.col("maxClinPhase") == 4) & (F.col("PhaseT") == "no"), F.lit("yes")
        ).otherwise(F.lit("no")),
    ).withColumn(
        "phase3Clean",
        F.when(
            (F.col("maxClinPhase") >= 3) & (F.col("PhaseT") == "no"), F.lit("yes")
        ).otherwise(F.lit("no")),
    ).withColumn(
        "phase2Clean",
        F.when(
            (F.col("maxClinPhase") >= 2) & (F.col("PhaseT") == "no"), F.lit("yes")
        ).otherwise(F.lit("no")),
    ).withColumn(
        "phase1Clean",
        F.when(
            (F.col("maxClinPhase") >= 1) & (F.col("PhaseT") == "no"), F.lit("yes")
        ).otherwise(F.lit("no")),
    )
            # 1. Get distinct values for the pivot column (essential for pivot())
        # This brings a small amount of data to the driver, but is necessary for the pivot schema.
        #distinct_pivot_values = [row[0] for row in test2.select(pivot_col_name).distinct().collect()]
        # print(f"Distinct values for '{pivot_col_name}': {distinct_pivot_values}")

        # 2. Perform the groupBy, pivot, and aggregate operations
        # The .pivot() function requires the list of distinct values for better performance
        # and correct schema inference.
        pivoted_df = (
            test2.groupBy(*group_by_columns)
            .pivot(pivot_col_name) # Provide distinct values distinct_pivot_values
            .agg(F.collect_set(F.col(agg_col_name))) # Collect all values into a set
            .fillna(0) # Fill cells that have no data with an empty list instead of null
        )
        # 3. Add items to dictionary to map the columns:
        # filter out None and 'null':
        datasetColumns=pivoted_df.columns
        filtered = [x for x in datasetColumns if x is not None and x != 'null']
        # using list comprehension
        for item in filtered:
            disdic[item] = pivot_col_name

        # 3. Add the 'data' literal column dynamically
        # This column indicates which aggregation column was used.
        #pivoted_df = pivoted_df.withColumn('data', F.lit(f'Drug_{agg_col_name}'))

        array_columns_to_convert = [
            field.name for field in pivoted_df.schema.fields
            if isinstance(field.dataType, ArrayType)
        ]
        print(f"Identified ArrayType columns for conversion: {array_columns_to_convert}")

        # 4. Apply the conversion logic to each identified array column
        df_after_conversion = pivoted_df # Start with the pivoted_df
        for col_to_convert in array_columns_to_convert:
            df_after_conversion = df_after_conversion.withColumn(
                col_to_convert,
                F.when(F.col(col_to_convert).isNull(), F.lit('no'))          # Handle NULLs (from pivot for no data)
                .when(F.size(F.col(col_to_convert)) == 0, F.lit('no'))       # Empty array -> 'no'
                .when(F.array_contains(F.col(col_to_convert), F.lit('yes')), F.lit('yes')) # Contains 'yes' -> 'yes'
                .when(F.array_contains(F.col(col_to_convert), F.lit('no')), F.lit('no'))   # Contains 'no' -> 'no'
                .otherwise(F.lit('no')) # Fallback for unexpected array content (e.g., ['other'], ['yes','no'])
            )

        # 4. Generate a unique name for this DataFrame and store it
        df_key = f"df_pivot_{agg_col_name.lower()}_by_{pivot_col_name.lower()}"
        all_pivoted_dfs[df_key] = df_after_conversion.withColumnRenamed( 'phase4Clean','Phase>=4'
        ).withColumnRenamed('phase3Clean','Phase>=3'
        ).withColumnRenamed('phase2Clean','Phase>=2'
        ).withColumnRenamed('phase1Clean','Phase>=1')

        # Exclude the reference columns
        other_cols = all_pivoted_dfs[df_key].columns[9:]
        # Loop over the other columns and add new ones
        for col in other_cols:
            new_col_name = f"tissueRelevant_{col}"
            all_pivoted_dfs[df_key] = all_pivoted_dfs[df_key].withColumn(
                new_col_name,
                F.when((F.col("Relevant") == "yes") & (F.col(col) == "yes"), F.lit("yes"))
                .otherwise(F.lit("no"))
            )


# --- Accessing your generated DataFrames ---
print("\n--- All generated DataFrames are stored in 'all_pivoted_dfs' dictionary ---")
print("Keys available:", all_pivoted_dfs.keys())

##### PROJECTID
biosample_keys = (
    benchmark
    .select("biosampleName")
    .distinct()
    .rdd
    .map(lambda r: r[0])
    .filter(lambda x: x is not None)  # <- remove NULLs
    .collect()
)
biosample_keys=[f"{k}_only" for k,v in disdic.items() if v == 'biosampleName']


###################################
###################################
result = []
result_st = []
result_ci = []
array2 = []
listado = []
result_all = []
today_date = str(date.today())

for key,df in all_pivoted_dfs.items():
#for key, df in list(all_pivoted_dfs.items())[:3]:
    print("Key:", key)
    print(f'working with {key}')
    parts = key.split('_by_') ### take the part of key belonging to column name
    column_name = parts[1] ### take the last part which is column name
    all_pivoted_dfs[key].persist()
    #unique_values = all_pivoted_dfs[key].drop('null').columns[7:]
    unique_values = all_pivoted_dfs[key].columns[:9] ### just the interesting columns for us 
    filtered_unique_values = [x for x in unique_values if x is not None and x != 'null']
    print('There are ', len(filtered_unique_values), 'columns to analyse with phases')
    rows = comparisons_df_iterative(filtered_unique_values)
    print(rows)
    print('printed all rows for', key)
    # If needed, now process the rest
    for row in rows:
        print('performing', row)
        results = aggregations_original(
            all_pivoted_dfs[key], key, listado, *row, today_date
        )
        result_all.append(results)
        print('results appended')
    all_pivoted_dfs[key].unpersist()
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
    ).withColumn(
    "folder",
    F.regexp_extract(F.col("path"), r"analysis/([^/]+)/", 1)
).withColumn(
    "suffix2",
    F.regexp_extract(F.col("folder"), r"df_pivot_(.+?)_by_", 1)
).withColumn(
    "prefix_type",
    F.when(
        F.col("prefix").rlike("^tissueRelevant_"),
        F.regexp_extract(F.col("prefix"), r"^(tissueRelevant)", 1)
    ).otherwise("single")
))

### annotate projectId, tissue, qtl type and doe type:

from pyspark.sql.functions import create_map
from itertools import chain

mapping_expr=create_map([F.lit(x) for x in chain(*disdic.items())])

df_annot=df.withColumn('annotation',mapping_expr.getItem(F.col('prefix')))


df_annot.toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_credibleSetColocDoEanalysis_filteredColocAndCaviarWithOthers4phasesTrue_testingRelevantTissue2.csv"
)

print("dataframe written \n Analysis finished")
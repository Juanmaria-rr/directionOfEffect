import time
from array import ArrayType
from functions import (
    relative_success,
    spreadSheetFormatter,
    discrepancifier,
    temporary_directionOfEffect,
    buildColocData,
    gwasDataset,
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

spark = SparkSession.builder.getOrCreate()
spark.conf.set(
    "spark.sql.shuffle.partitions", "400"
)  # Default is 200, increase if needed

print('joint groups')
path_n='gs://open-targets-data-releases/25.06/output/'

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

print("loaded files")

#### FIRST MODULE: BUILDING COLOC 
newColoc=buildColocData(all_coloc,credible,index)

print("loaded newColoc")

### SECOND MODULE: PROCESS EVIDENCES TO AVOID EXCESS OF COLUMNS 
gwasComplete = gwasDataset(evidences,credible)

#### THIRD MODULE: INCLUDE COLOC IN THE 
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


print("built drugApproved dataset")


#### FOURTH MODULE BUILDING CHEMBL ASSOCIATIONS - HERE TAKE CARE WITH FILTERING STEP 
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
    #.filter(F.col("coherencyDiagonal") == "coherent")
    .drop(
        "coherencyDiagonal", "coherencyOneCell", "noEvaluable", "GoF_risk", "LoF_risk"
    )
    .withColumnRenamed("GoF_protect", "drugGoF_protect")
    .withColumnRenamed("LoF_protect", "drugLoF_protect")
    # .persist()
)

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

#rightTissue = spark.read.csv(
#    'gs://ot-team/jroldan/analysis/20250526_rightTissue.csv',
#    header=True,
#).drop("_c0")

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

print("built bench2 dataset")

###### cut from here
print("looping for variables_study")

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
benchmark = (
    (
        resolvedColocFiltered.filter(F.col("betaGwas") < 0).filter(
        F.col("name") != "COVID-19"
    )
        .join(  ### select just GWAS giving protection
            analysis_chembl_indication, on=["targetId", "diseaseId"], how="right"  ### RIGHT SIDE
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
    )  #### remove COVID-19 associations
).join(biosample.select("biosampleId", "biosampleName"), on="biosampleId", how="left")

#bench2 = benchmark.join(
#    rightTissue, on=["name", "bioSampleName"], how="left"
#).withColumn(
#    "rightTissue",
#    F.when(F.col("rightTissue1") == "yes", F.lit("yes")).otherwise(F.lit("no")),
#)

print("built benchmark dataset")

## write the benchmark 
#name='benchmark'
#output_partitioned_path = f"gs://ot-team/jroldan/analysis/parquetFiles/{name}"
#benchmark.write.mode("overwrite").parquet(output_partitioned_path)
#print(f'written {name}')
#### Analysis

#### 1 Build a dictionary with the distinct values as key and column names as value
variables_study = ["projectId", "biosampleName", "rightStudyType", "colocDoE","colocalisationMethod"]

# List to hold temporary DataFrames
temp_dfs_for_union = []

# Iterate over the column names to prepare DataFrames for union
for col_name in variables_study:
    # Select the current column, alias it to 'distinct_value' for consistent schema
    # Filter out nulls, then get distinct values
    # Add a literal column with the original 'col_name'
    df_temp = (
        benchmark.select(F.col(col_name).alias("distinct_value"))
        .filter(F.col("distinct_value").isNotNull()) # Exclude None (null) values
        .distinct()
        .withColumn("column_name", F.lit(col_name))
    )
    temp_dfs_for_union.append(df_temp)

disdic = {}

if temp_dfs_for_union:
    # Union all the temporary DataFrames.
    # unionByName is crucial to handle potential schema differences (e.g., if columns have same name but different types)
    # and ensures columns are matched by name.
    combined_distinct_values_df = temp_dfs_for_union[0]
    for i in range(1, len(temp_dfs_for_union)):
        combined_distinct_values_df = combined_distinct_values_df.unionByName(temp_dfs_for_union[i])

    # Now, collect the combined distinct values.
    # This is a single collect operation on the aggregated DataFrame.
    print("Collecting combined distinct values from the cluster...")
    collected_rows = combined_distinct_values_df.collect()

    # Populate the dictionary from the collected rows
    for row in collected_rows:
        disdic[row.distinct_value] = row.column_name
else:
    print("variables_study list is empty, disdic will be empty.")


print("\nFinal disdic:", disdic)

# Assuming 'spark' session, 'benchmark' DataFrame, 'negativeTD' DataFrame, and 'disdic' dictionary are defined

# --- Step 1: Pre-compute 'hasboth' ONCE ---
# This is a shuffle, but only happens once.
print("Pre-computing 'hasboth' column...")
window_target_disease_only = Window.partitionBy('targetId', 'diseaseId')
benchmark_processed = benchmark.withColumn(
    'hasboth',
    F.size(F.collect_set('colocalisationMethod').over(window_target_disease_only))
)

# You might consider caching this intermediate result if 'benchmark' is very large
# and you have enough memory, to avoid re-reading from source if possible.
# benchmark_processed.cache() # or .persist(StorageLevel.MEMORY_AND_DISK)
# benchmark_processed.count() # Force computation if you cache

pivoted_dfs = {}

# --- Step 2: Loop for each variable_study column ---
for col_name in variables_study:
    print(f"Processing pivot for: {col_name}")

    # Define window specs for the current iteration, including 'col_name' in partition
    # (This shuffle is still per iteration, but unavoidable if 'resolvedAgreeDrug' depends on 'col_name' values)
    current_col_window_spec_qtl = Window.partitionBy("targetId", "diseaseId", col_name).orderBy(F.col("qtlPValueExponent").asc())
    current_col_pvalue_order_window = Window.partitionBy("targetId", "diseaseId", col_name).orderBy(F.col('colocalisationMethod').asc(), F.col("qtlPValueExponent").asc())

    # Calculate 'resolvedAgreeDrug' for the current 'col_name'
    # This involves a shuffle per iteration.
    temp_df_with_resolved = benchmark_processed.withColumn('resolvedAgreeDrug',
        F.when(F.col('hasboth') > 1,
            F.first(F.col('AgreeDrug'), ignorenulls=True).over(current_col_pvalue_order_window)
        ).otherwise(F.first(F.col('AgreeDrug'), ignorenulls=True).over(current_col_window_spec_qtl))
    )

    # --- Step 3: Perform the pivot and join ---
    # This is an expensive operation (shuffle, potential wide dataframe)
    pivoted_df = (
        temp_df_with_resolved
        .groupBy(
            "targetId",
            "diseaseId",
            "maxClinPhase",
        )
        .pivot(col_name) # Pivoting on values of the 'col_name' column
        .agg(F.collect_set("resolvedAgreeDrug"))
        .join(negativeTD, on=["targetId", "diseaseId"], how="left") # Ensure negativeTD is broadcast if small
    )

    # --- Step 4: Add derived columns (these are generally cheap) ---
    for phase in [1, 2, 3, 4]:
        pivoted_df = pivoted_df.withColumn(
            f"Phase>={phase}",
            F.when(F.col("maxClinPhase") >= phase, F.lit("yes")).otherwise(F.lit("no")),
        )

    pivoted_df = pivoted_df.withColumn(
        "PhaseT",
        F.when(F.col("stopReason") == "Negative", F.lit("yes")).otherwise(F.lit("no")),
    ).withColumn(
        "Phase>=4",
        F.when(
            (F.col("maxClinPhase") == 4) & (F.col("PhaseT") == "no"), F.lit("yes")
        ).otherwise(F.lit("no")),
    ).withColumn(
        "Phase>=3",
        F.when(
            (F.col("maxClinPhase") >= 3) & (F.col("PhaseT") == "no"), F.lit("yes")
        ).otherwise(F.lit("no")),
    ).withColumn(
        "Phase>=2",
        F.when(
            (F.col("maxClinPhase") >= 2) & (F.col("PhaseT") == "no"), F.lit("yes")
        ).otherwise(F.lit("no")),
    ).withColumn(
        "Phase>=1",
        F.when(
            (F.col("maxClinPhase") >= 1) & (F.col("PhaseT") == "no"), F.lit("yes")
        ).otherwise(F.lit("no")),
    )

    # Add _only columns dynamically based on disdic values matching current column
    matching_keys = [key for key, val in disdic.items() if val == col_name]

    for key in matching_keys:
        # F.col(key) assumes 'key' refers to a column that exists in pivoted_df after the pivot.
        pivoted_df = pivoted_df.withColumn(
            f"{key}_only",
            F.when(F.array_contains(F.col(key), "yes"), F.lit("yes")).otherwise(F.lit("no")),
        )

### making columns for the 

    # --- Step 5: Store result. Consider writing to GCS to break lineage if memory is an issue ---
    # This is highly recommended if 'variables_study' is very large.
    # Write to Parquet for efficient storage and schema preservation.
    # output_path = f"gs://your-bucket/temp_pivoted_results/{col_name}"
    # print(f"Writing results for {col_name} to {output_path}")
    # pivoted_df.write.mode("overwrite").parquet(output_path)
    # pivoted_dfs[col_name] = spark.read.parquet(output_path) # Read back if needed later
    # output_partitioned_path = f"gs://ot-team/jroldan/analysis/parquetFiles/pivoted_df_{col_name}"
    # pivoted_df.write.mode("overwrite").parquet(output_partitioned_path)
    # print(f"DataFrame successfully written and partitioned to {output_partitioned_path}")
    # If not writing to GCS, just store the DF in memory (be cautious for large number of DFs)

    pivoted_dfs[col_name] = pivoted_df

##### PROJECTID
project_keys=[f"{k}_only" for k,v in disdic.items() if v == 'projectId']
main=['GTEx_only', 'UKB_PPP_EUR_only']
stimulated=['Alasoo_2018_only','Cytoimmgen_only','Fairfax_2014_only','Kim-Hellmuth_2017_only','Nathan_2022_only','Nedelec_2016_only','Quach_2016_only','Randolph_2021_only','Schmiedel_2018_only']
cellLine=['CAP_only','HipSci_only','iPSCORE_only','Jerber_2021_only','PhLiPS_only','Schwartzentruber_2018_only','TwinsUK_only']

others=[item for item in project_keys if item not in main]
nonStimulated=[item for item in project_keys if item not in stimulated]
otherCellLine=[item for item in project_keys if item not in cellLine]

# First condition: any "yes" in list1
condition1 = reduce(lambda acc, col: acc | (F.col(col) == "yes"), others[1:], F.col(others[0]) == "yes")
# estimulated
condition2 = reduce(lambda acc, col: acc | (F.col(col) == "yes"), stimulated[1:], F.col(stimulated[0]) == "yes")
## non estimulated:
condition3 = reduce(lambda acc, col: acc | (F.col(col) == "yes"), nonStimulated[1:], F.col(nonStimulated[0]) == "yes")
# cellLine
condition4 = reduce(lambda acc, col: acc | (F.col(col) == "yes"), cellLine[1:], F.col(cellLine[0]) == "yes")
# non cellline
condition5 = reduce(lambda acc, col: acc | (F.col(col) == "yes"), otherCellLine[1:], F.col(otherCellLine[0]) == "yes")
# non cellline
condition6 = reduce(lambda acc, col: acc | (F.col(col) == "yes"), main[1:], F.col(main[0]) == "yes")

# Add both columns
pivoted_dfs['projectId'] = pivoted_dfs['projectId'].withColumn("othersProjectId_only", F.when(condition1, "yes").otherwise("no")) 
pivoted_dfs['projectId'] = pivoted_dfs['projectId'].withColumn("estimulated_only", F.when(condition2, "yes").otherwise("no")) 
pivoted_dfs['projectId'] = pivoted_dfs['projectId'].withColumn("nonStimulated_only", F.when(condition3, "yes").otherwise("no")) 
pivoted_dfs['projectId'] = pivoted_dfs['projectId'].withColumn("cellLine", F.when(condition4, "yes").otherwise("no")) 
pivoted_dfs['projectId'] = pivoted_dfs['projectId'].withColumn("nonCellLine", F.when(condition5, "yes").otherwise("no")) 
pivoted_dfs['projectId'] = pivoted_dfs['projectId'].withColumn("GTExUKB", F.when(condition6, "yes").otherwise("no")) 


###append to dictionary

disdic.update({'othersProjectId': 'projectId','Stimulated': 'projectId','cellLine': 'projectId', 'othersBiosampleName_only': 'biosampleName', 'otherRightStudyType':'rightStudyType'})


result = []
result_st = []
result_ci = []
array2 = []
listado = []
result_all = []
today_date = str(date.today())

##### PROJECT ID ###### 
print('working with projectId')
pivoted_dfs['projectId'].persist()
unique_values = benchmark.select('projectId').filter(F.col('projectId').isNotNull()).distinct().rdd.flatMap(lambda x: x).collect()
filter = len(pivoted_dfs['projectId'].drop(*unique_values).columns[10:])
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
    f"gs://ot-team/jroldan/analysis/{today_date}_credibleSetColocDoEanalysis_filteredColocAndCaviarWithOthers4phases_jointProjects.csv"
)

print("dataframe written \n Analysis finished")
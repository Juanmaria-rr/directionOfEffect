### --- REPRODUCIBILITY LOGGING (added, does not touch analysis logic) ---
# Writes a .log file next to this script that contains:
#   1) a full snapshot of the source code that was actually run
#   2) everything printed during that run (stdout/stderr), timestamped
# so a failed/successful run can always be traced back to the exact code that produced it.
import sys as _sys
import os as _os
from datetime import datetime as _dt

_SCRIPT_PATH = _os.path.abspath(__file__)
_LOG_DIR = _os.path.join(_os.path.dirname(_SCRIPT_PATH), "logs")
_os.makedirs(_LOG_DIR, exist_ok=True)
_RUN_TS = _dt.now().strftime("%Y-%m-%d_%H-%M-%S")
_LOG_PATH = _os.path.join(
    _LOG_DIR, f"{_os.path.splitext(_os.path.basename(_SCRIPT_PATH))[0]}_{_RUN_TS}.log"
)


class _TeeStream:
    """Duplicates writes to the original stream and to the log file."""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for s in self._streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self._streams:
            s.flush()


_log_file = open(_LOG_PATH, "w")
with open(_SCRIPT_PATH, "r") as _src:
    _log_file.write("=" * 80 + "\n")
    _log_file.write(f"SOURCE CODE SNAPSHOT: {_SCRIPT_PATH}\n")
    _log_file.write(f"Run started: {_dt.now().isoformat()}\n")
    _log_file.write("=" * 80 + "\n\n")
    _log_file.write(_src.read())
    _log_file.write("\n\n" + "=" * 80 + "\nEXECUTION OUTPUT\n" + "=" * 80 + "\n\n")
_log_file.flush()

_sys.stdout = _TeeStream(_sys.__stdout__, _log_file)
_sys.stderr = _TeeStream(_sys.__stderr__, _log_file)

print(f"📝 Logging full run (code + output) to: {_LOG_PATH}")
### --- END REPRODUCIBILITY LOGGING ---

import time
import os
from pyspark.sql import SparkSession, Window
from pyspark import SparkConf

print('new config for spark session')
print('🔧 Initializing Spark with optimized single-node configuration')

# functions.py and DoEAssessment.py both do `spark = SparkSession.builder.getOrCreate()`
# at import time with NO conf. If they're imported before we create our own session,
# THEIR bare getOrCreate() wins and silently locks in the cluster-wide defaults
# (dynamicAllocation.enabled=true, maxExecutors=10000) — which is what let YARN spin up
# several 28g executors on this single-node box, oversubscribing its memory and causing
# the OOM container kills / heartbeat timeouts seen in past runs (and confirmed live:
# trying to stop+recreate the session afterwards instead breaks YARN AM registration).
# Creating the session with our conf FIRST, before those imports, means their getOrCreate()
# simply attaches to the session we already configured correctly.

# Make sure /tmp spill dir exists
os.makedirs("/tmp/spark-temp", exist_ok=True)

# chmod 1777 /tmp/spark-temp if there are problems like: NameError: name 'os' is not defined
conf = (
    SparkConf()
    .setAppName("directionOfEffect")
    # Memory tuning
    .set("spark.driver.memory", "8g")
    .set("spark.executor.memory", "24g")
    .set("spark.executor.memoryOverhead", "4g")
    .set("spark.memory.fraction", "0.5")
    .set("spark.memory.storageFraction", "0.3")
    # CPU
    .set("spark.executor.instances", "1")
    .set("spark.executor.cores", "8")
    # This is a single-node cluster: force a single static executor so it can't
    # be oversubscribed. The cluster-wide spark-defaults.conf sets
    # dynamicAllocation.enabled=true / maxExecutors=10000, which — whenever it
    # actually takes effect — lets YARN launch multiple 28g executors that don't
    # fit together in the node's memory, causing the OOM container kills (exit 137)
    # and executor heartbeat timeouts seen in past runs.
    .set("spark.dynamicAllocation.enabled", "false")
    .set("spark.dynamicAllocation.minExecutors", "1")
    .set("spark.dynamicAllocation.maxExecutors", "1")
    # Give transient GC pauses more room before YARN/Spark declares an executor dead
    .set("spark.network.timeout", "800s")
    .set("spark.executor.heartbeatInterval", "30s")
    # Shuffle and partitions
    .set("spark.sql.shuffle.partitions", "32")
    .set("spark.default.parallelism", "32")
    # Disable AQE if unstable
    .set("spark.sql.adaptive.enabled", "false")
    # Disk spill directories
    .set("spark.local.dir", "/tmp/spark-temp")
    # Avoid broadcast join overload
    .set("spark.sql.autoBroadcastJoinThreshold", "-1")
    # GCS reliability
    .set("spark.hadoop.fs.gs.outputstream.upload.retry.max.retry.limit", "8")
    .set("spark.hadoop.fs.gs.outputstream.upload.chunk.size", "16777216")
    # increase cap maximum size for writting a file 
    .set("spark.rpc.message.maxSize", "1024")

)

spark = SparkSession.builder.config(conf=conf).getOrCreate()
spark.sparkContext.setLogLevel("WARN")
# This ensures Spark will spill large data to disk instead of crashing when memory fills.
spark.conf.set("spark.sql.execution.arrow.pyspark.enabled", "true")
spark.conf.set("spark.sql.execution.arrow.maxRecordsPerBatch", "200000")
spark.conf.set("spark.sql.shuffle.spill", "true")
spark.conf.set("spark.storage.memoryFraction", "0.3")

print("✅ Spark session started with optimized config")

# Deferred until after our SparkSession exists (see note above): both modules do a
# bare `SparkSession.builder.getOrCreate()` at import time, which now attaches to
# the session we just configured instead of creating a conflicting default one.
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
import pyspark.sql.functions as F
from datetime import datetime
from datetime import date
from pyspark.sql.types import StructType, StructField, StringType,ArrayType, IntegerType
from pyspark.sql.types import (
    StructType,
    StructField,
    DoubleType,
    DecimalType,
    StringType,
    FloatType,
)
from pyspark.sql.functions import create_map
from itertools import chain
import pandas as pd
from functools import reduce

###2### modification two: coalesce before writting
def safe_parquet_write(df, path, mode="overwrite"):
    """Write Parquet efficiently, adjusting partitions by DataFrame size."""
    row_count = df.count()
    n_partitions = 1 if row_count < 5_000_000 else max(2, df.rdd.getNumPartitions() // 4)
    print(f"🪶 Writing {row_count:,} rows → {n_partitions} partition(s) → {path}")

    (
        df.coalesce(n_partitions)
        .write.mode(mode)
        .option("compression", "snappy")
        .parquet(path)
    )

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

all_coloc=ecaviar.unionByName(new, allowMissingColumns=True).filter((F.col('clpp')>=0.01) | (F.col('h4')>=0.8))

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

resolvedColocFiltered = resolvedColoc

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

print("built drugApproved dataset")
benchmark = (
        resolvedColocFiltered.filter( ## .filter(F.col("betaGwas") < 0)
        F.col("name") != "COVID-19"
    )
        .join(  ### select just GWAS giving protection
            analysis_chembl_indication, on=["targetId", "diseaseId"], how="right"  ### RIGHT SIDE
        )
).join(biosample.select("biosampleId", "biosampleName"), on="biosampleId", how="left").join(diseaseTA, on='diseaseId',how='left')

print("built benchmark")

#### HERE THE CODE FOR THE ANALYSIS

####2 Define agregation function
#### ===============================================================
#### FIXED aggregation function — avoids duplication & adds logging
#### ===============================================================

from pyspark.sql import Window
import pyspark.sql.functions as F
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

def convertTuple(tup):
    return ",".join(map(str, tup))


def aggregations_original(
    df,
    data,
    listado,
    comparisonColumn,
    comparisonType,
    predictionColumn,
    predictionType,
    c,
    all_outputs,
    uniqIds,
):
    """
    Compute stats for one (comparison, prediction) combination.

    Returns only the *new* results as a list (no global list mutation).
    """

    local_results = []  # <-- local container to avoid duplication

    # Define windows
    wComparison = Window.partitionBy(comparisonColumn)
    wPrediction = Window.partitionBy(predictionColumn)
    wPredictionComparison = Window.partitionBy(comparisonColumn, predictionColumn)

    # uniqIds is constant for the whole `df` (doesn't depend on comparison/prediction),
    # so it's computed once per key by the caller instead of once per row here —
    # this was previously a full distinct().count() shuffle repeated 1590 times.

    # Build the aggregated output
    out = (
        df.withColumn("comparisonType", F.lit(comparisonType))
        .withColumn("predictionType", F.lit(predictionType))
        .withColumn("total", F.lit(uniqIds))
        .withColumn("a", F.count("targetId").over(wPredictionComparison))
        .withColumn("predictionTotal", F.count("targetId").over(wPrediction))
        .withColumn("comparisonTotal", F.count("targetId").over(wComparison))
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

    # Add metadata
    data_label = f"df_{data}_{comparisonColumn}_{predictionColumn}.parquet"
    out = out.withColumn("data", F.lit(data_label))
    print(f"📦  {data_label}   —   time: {c}")

    # Optionally keep spark DF reference for later (if needed)
    all_outputs.append(out)

    # --- Statistical analysis ---
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

    (rs_result, rs_ci) = relative_success(array1)

    local_results.append(
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
            data_label,
        ]
    )

    return local_results



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

print('ingesting curated dataset of diseases-biosampleName')
file_path='gs://ot-team/jroldan/test_20250927_2509_biosampleName_curatedOn_20251007.xls'
# Read Excel file with pandas
pdf = pd.read_excel(file_path)
# Convert to Spark DataFrame,
rightTissue = spark.createDataFrame(pdf).drop('double check', 'second check','Relevant').withColumnRenamed('Revised', 'Relevant')
print("loaded rightTissue dataset")

#### Do the join with Right Tissue annotation

benchmark_tissue=benchmark.join(rightTissue.drop('name'), on=['diseaseId','biosampleName'], how='left')

### create disdic dictionary
disdic={}

# --- Configuration for your iterative pivoting ---
group_by_columns = ['targetId', 'diseaseId','Relevant','phase4Clean','phase3Clean','phase2Clean','phase1Clean','PhaseT']
#columns_to_pivot_on = ['actionType2', 'biosampleName', 'projectId', 'rightStudyType','colocalisationMethod']
columns_to_pivot_on = ['biosampleName', 'rightStudyType','taLabel','projectId','colocalisationMethod']
columns_to_aggregate = ['NoneCellYes', 'NdiagonalYes','hasGenetics'] # The values you want to collect in the pivoted cells
all_pivoted_dfs = {}

doe_columns=["LoF_protect", "GoF_risk", "LoF_risk", "GoF_protect"]
diagonal_lof=['LoF_protect','GoF_risk']
diagonal_gof=['LoF_risk','GoF_protect']

conditions = [
    F.when(F.col(c) == F.col("maxDoE"), F.lit(c)).otherwise(F.lit(None)) for c in doe_columns
    ]


#### make column combinations
from pyspark.sql import functions as F
from functools import reduce
from itertools import combinations
### other cols rightStudyType


pesceQTL=['eqtl','pqtl','sceqtl']
alleQTL=['eqtl','sceqtl']
peQTL=['eqtl','pqtl']
tusQTL=['tuqtl','sqtl']

subgroupsQTL=['pesceQTL','peQTL','alleQTL','tusQTL']

anyImmuneCell=[
"T cell",
"B cell",
"natural killer cell",
"CD8-positive, alpha-beta T cell",
"CD14-low, CD16-positive monocyte",
"CD14-positive, CD16-negative classical monocyte",
"CD16-negative, CD56-bright natural killer cell, human",
"CD4-positive, alpha-beta T cell",
"CD4-positive, alpha-beta cytotoxic T cell",
"CD4-positive, alpha-beta memory T cell",
"T follicular helper cell",
"T-helper 1 cell",
"T-helper 17 cell",
"T-helper 2 cell",
"naive regulatory T cell",
"central memory CD4-positive, alpha-beta T cell",
"central memory CD8-positive, alpha-beta T cell",
#"lymphoblastoid cell line",
"effector memory CD4-positive, alpha-beta T cell",
"effector memory CD8-positive, alpha-beta T cell",
"gamma-delta T cell",
#hematopoietic precursor cell
"macrophage",
"memory B cell",
"memory regulatory T cell",
"microglial cell", ## primary immmune cells from CNS
"mucosal invariant T cell",
"neutrophil",
#"plasmablast",
"plasmacytoid dendritic cell"]

Tcells=[
"T cell",
"CD8-positive, alpha-beta T cell",
"CD4-positive, alpha-beta T cell",
"CD4-positive, alpha-beta cytotoxic T cell",
"CD4-positive, alpha-beta memory T cell",
"T follicular helper cell",
"T-helper 1 cell",
"T-helper 17 cell",
"T-helper 2 cell",
"naive regulatory T cell",
"central memory CD4-positive, alpha-beta T cell",
"central memory CD8-positive, alpha-beta T cell",
#"lymphoblastoid cell line",
"effector memory CD4-positive, alpha-beta T cell",
"effector memory CD8-positive, alpha-beta T cell",
"gamma-delta T cell",
"memory regulatory T cell",
"mucosal invariant T cell",
]

Bcells=[
"B cell",
"memory B cell",
]

others=[
"natural killer cell",
"macrophage",
"neutrophil",
"plasmacytoid dendritic cell"
]
## cols to make the union for:
immuneGroups=['anyImmuneCell','Tcells','Bcells','others']

from functools import reduce
from pyspark.sql import functions as F

def add_combination_column(df, combo, mode):
    """
    Adds a new binary column ("yes"/"no") for the combination of cols.
    Column name will be <value>_<mode>.
    """
    # If combo is a string (e.g. "Bcells"), resolve it into the actual list
    if isinstance(combo, str):
        list_name = combo  # save original string for naming
        combo = globals()[combo]  # get the actual list from globals
    else:
        raise ValueError("Expected 'combo' to be a string (name of the list)")

    combo = [str(c) for c in combo]
    combo_name = f"{list_name}_{mode}"  # use value as base of new colname

    if mode == "union":
        expr = reduce(lambda a, b: a | b, [(F.col(c) == "yes") for c in combo])
    elif mode == "inter":
        expr = reduce(lambda a, b: a & b, [(F.col(c) == "yes") for c in combo])
    elif mode == "unionTissue":
        expr = reduce(lambda a, b: a | b, [(F.col(f"tissueRelevant_{c}") == "yes") for c in combo])
    elif mode == "interTissue":
        expr = reduce(lambda a, b: a & b, [(F.col(f"tissueRelevant_{c}") == "yes") for c in combo])
    elif mode == "unionNonTissue":
        expr = reduce(lambda a, b: a | b, [(F.col(f"nonTissueRelevant_{c}") == "yes") for c in combo])
    elif mode == "interNonTissue":
        expr = reduce(lambda a, b: a & b, [(F.col(f"nonTissueRelevant_{c}") == "yes") for c in combo])
    else:
        raise ValueError("mode must be 'union', 'inter', 'unionTissue', 'interTissue','unionNonTissue' or 'interNonTissue'")

    return (
        df.withColumn(combo_name, F.when(expr, F.lit("yes")).otherwise(F.lit("no"))),
        combo_name,
    )

for pivot_col_name in columns_to_pivot_on:
    for agg_col_name in columns_to_aggregate:
        print(f"\n--- Creating DataFrame for Aggregation: '{agg_col_name}' and Pivot: '{pivot_col_name}' ---")
        # TODO: qtlColocDoE below is computed (best QTL signal per targetId/diseaseId/maxClinPhase/
        # Relevant/pivot_col_name group, lowest p-value, COLOC preferred over eCAVIAR — same criterion
        # genEvidDoE_colocAndCaviar_new.py uses globally at its "resolvedAgreeDrug"/"homogenized" step)
        # but the .pivot() below still runs on the raw "colocDoE" column, not on qtlColocDoE. Effect:
        # when a group has >1 QTL row with disagreeing colocDoE (e.g. two eQTL datasets in the same
        # biosample), both get counted here, and discrepancifier() then marks the group "dispar" purely
        # from that redundancy — job1 can't produce that kind of "dispar" because it already collapsed
        # to one signal before this point. Decide whether to switch .pivot("colocDoE") to
        # .pivot("qtlColocDoE") (would make this collapse per-group instead of globally, matching job1's
        # methodology) and re-run job2 to see the impact on NoneCellYes/NdiagonalYes proportions.
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
        for col in other_cols:
            new_col_name = f"nonTissueRelevant_{col}"
            all_pivoted_dfs[df_key] = all_pivoted_dfs[df_key].withColumn(
                new_col_name,
                F.when((F.col("Relevant") == "no") & (F.col(col) == "yes"), F.lit("yes"))
                .otherwise(F.lit("no"))
            )

        # 🚨 Guard clause: if pivot_col_name not recognized, skip it
        if pivot_col_name not in ["rightStudyType", "biosampleName"]:
            # Option 1: skip only this iteration (if in an outer loop)
            continue
            # Option 2 (if you are inside a function and want to stop entirely):
            # return

        # ✅ For recognized pivot types, run combinations
        if pivot_col_name == "rightStudyType":
            combos = subgroupsQTL
        elif pivot_col_name == "biosampleName":
            combos = immuneGroups

        for value in combos:
            for mode in ["union", "unionTissue", "unionNonTissue", "inter", "interTissue", "interNonTissue"]:
                all_pivoted_dfs[df_key], _ = add_combination_column(all_pivoted_dfs[df_key], value, mode)

# --- Accessing your generated DataFrames ---
print("\n--- All generated DataFrames are stored in 'all_pivoted_dfs' dictionary ---")
print("Keys available:", all_pivoted_dfs.keys())


###################################
# Iterating over all pivoted dataframes
###################################
listado = []
result_all = []
today_date = str(date.today())
run_tag = datetime.now().strftime("%H-%M")  # start-time suffix so a same-day re-run doesn't overwrite a previous one
all_outputs = []

for key, df in all_pivoted_dfs.items():
    print(f"\n🔹 Working on {key}")
    parts = key.split('_by_')
    column_name = parts[1]

    df.persist()
    unique_values = df.columns[9:]
    filtered_unique_values = [x for x in unique_values if x and x != "null"]

    print(f"   Found {len(filtered_unique_values)} columns to analyse with phases")
    rows = comparisons_df_iterative(filtered_unique_values)
    print(f"   → {len(rows)} comparisons to test")

    # Computed once per key instead of once per row (see aggregations_original):
    # it doesn't depend on comparisonColumn/predictionColumn, only on `df`.
    uniqIds = df.select("targetId", "diseaseId").distinct().count()

    for i, row in enumerate(rows, 1):
        print(f"      ▶ ({i}/{len(rows)}) running {row}")
        new_results = aggregations_original(
            df, key, listado, *row, today_date, all_outputs, uniqIds
        )
        result_all.extend(new_results)
        print(f"         +{len(new_results)} new results (total {len(result_all)})")

    df.unpersist()
    print(f"✅ Finished {key} and released memory\n")



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
    "tissueRelevant_",
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


def safe_parquet_write(df, path, mode="overwrite"):
    row_count = df.count()
    n_partitions = 1 if row_count < 1_000_000 else 4
    print(f"🪶 Writing {row_count:,} rows → {n_partitions} partition(s) → {path}")
    
    (
        df.coalesce(n_partitions)
        .write.mode(mode)
        .option("compression", "snappy")
        .parquet(path)
    )

mapping_expr=create_map([F.lit(x) for x in chain(*disdic.items())])

df_annot=df.withColumn('annotation',mapping_expr.getItem(F.col('prefix'))).withColumn(
    "part1",
    F.regexp_extract(F.col("path"), r"pivot_([^_]+)_by_", 1)
).withColumn(
    "part2",
    F.regexp_extract(F.col("path"), r"_by_([^_]+)_", 1)
).withColumn(
    "part3",
    F.regexp_extract(F.col("path"), r"_by_[^_]+_(.*?)_Phase", 1)
).withColumn(
    "part4",
    F.regexp_extract(F.col("path"), r"_(Phase[^.]+)\.parquet", 1)
).withColumn(
    "part5",
    F.when(F.col("part3").rlike("(?i)tissueRelevant"), F.lit("tissuerelevant"))
     .otherwise(F.lit("single"))
)

output_path = f"gs://ot-team/jroldan/analysis/{today_date}_{run_tag}_ColocCaviar_AllCombinations.parquet"

safe_parquet_write(df_annot, output_path)
print(f"✅ Final results written as parquet to {output_path}")

print('reading dataframe from parquet file')
totsv=spark.read.parquet(f"{output_path}")

totsv.toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_{run_tag}_ColocCaviar_AllCombinations.tsv", sep="\t", index=False)

print("dataframe written on tsv format \n Analysis finished")
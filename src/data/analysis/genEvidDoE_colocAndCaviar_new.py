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
    build_resolved_coloc_noPropag
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

diseases2 = diseases.select("id", "parents").withColumn(
    "diseaseIdPropagated",
    F.explode_outer(F.concat(F.array(F.col("id")), F.col("parents"))),
)


print("loaded files")

#### FIRST MODULE: BUILDING COLOC 
newColoc=buildColocData(all_coloc,credible,index)

print("loaded newColoc")

### SECOND MODULE: PROCESS EVIDENCES TO AVOID EXCESS OF COLUMNS 
gwasComplete = gwasDataset(evidences,credible)

print('gwasComplete loaded')
#### THIRD MODULE: INCLUDE COLOC IN THE 
# In here we use Coloc noPropag
resolvedColoc=build_resolved_coloc_noPropag(newColoc, gwasComplete,diseases).filter( ## .filter(F.col("betaGwas") < 0)
        F.col("name") != "COVID-19"
    ).withColumn('hasGenetics', F.lit('yes'))

resolvedColocFiltered = resolvedColoc.filter((F.col('clpp')>=0.01) | (F.col('h4')>=0.8))

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

##### FROM NOW IS THE OLD ANALYSIS OF GEN EVIDENCE 
# Define window specs for the current iteration, including 'col_name' in partition
# (This shuffle is still per iteration, but unavoidable if 'resolvedAgreeDrug' depends on 'col_name' values)
current_col_window_spec_qtl = Window.partitionBy("targetId", "diseaseId").orderBy(
    F.col("qtlPValueExponent").asc()
)
current_col_pvalue_order_window = Window.partitionBy("targetId", "diseaseId").orderBy(
    F.col("colocalisationMethod").asc(), F.col("qtlPValueExponent").asc()
)

### if the are two coloc method, say it
window_target_disease_only = Window.partitionBy("targetId", "diseaseId")
benchmark_processed = resolvedColocFiltered.withColumn(
    "hasboth",
    F.size(F.collect_set("colocalisationMethod").over(window_target_disease_only)),
)

# take the best coloc (Coloc > ecaviar) with the lowest p value
gwasCredibleAssoc_qtlPValue = benchmark_processed.withColumn(
    "resolvedAgreeDrug",
    F.when(
        F.col("hasboth") > 1,
        F.first(F.col("colocDoE"), ignorenulls=True).over(
            current_col_pvalue_order_window
        ),
    ).otherwise(
        F.first(F.col("colocDoE"), ignorenulls=True).over(current_col_window_spec_qtl)
    ),
    ).withColumn(
        "homogenized", F.col('colocDoE')
    ).withColumn(
        "homogenized",
        F.when(F.col("homogenized").isNull(), F.lit("noEvaluable")).otherwise(
            F.col("homogenized")
        ),
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


assessment_all = assessment.unionByName(
    gwasCredibleAssoc_qtlPValue.withColumn("datasourceId", F.lit("gwas_credible_set")),
    allowMissingColumns=True,
)

print("defining non propagated,propagated and analysis_drugs functions")

def analysis_nonPropagated(assessment_all, analysisDatasources):
    return discrepancifier(
        assessment_all.filter(F.col("datasourceId").isin(analysisDatasources))
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


def analysis_propagated(assessment_all, analysisDatasources):
    return discrepancifier(
        assessment_all.filter(F.col("datasourceId").isin(analysisDatasources))
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

def analysis_drugs(assessment_all, chembl_ds):
    return discrepancifier(
        assessment_all.filter((F.col("datasourceId").isin(chembl_ds))
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


analysis_chembl = analysis_drugs(assessment_all, chembl_ds)

#######
## include here the analysis
#######

analysisDatasources = []

print("defining full_analysis_propagation")

doe_columns=["LoF_protect", "GoF_risk", "LoF_risk", "GoF_protect"]
diagonal_lof=['LoF_protect','GoF_risk']
diagonal_gof=['LoF_risk','GoF_protect']

def full_analysis_propagation(
    doe_columns,assessment_all, analysisDatasources, analysis_chembl, negativeTD, diseaseTA,diagonal_lof,diagonal_gof
):
    conditions = [
    F.when(F.col(c) == F.col("maxDoE"), F.lit(c)).otherwise(F.lit(None)) for c in doe_columns
    ]
    
    return (
        analysis_propagated(assessment_all, analysisDatasources)
        .join(
            analysis_chembl.filter(F.col("coherencyDiagonal") == "coherent").selectExpr(
                "targetId",
                "diseaseId",
                "maxClinPhase",
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
        ).withColumn(
            "arrayN", F.array(*[F.col(c) for c in doe_columns])
        ).withColumn(
            "maxDoE", F.array_max(F.col("arrayN"))
        ).withColumn("maxDoE_names", F.array(*conditions)
        ).withColumn("maxDoE_names", F.expr("filter(maxDoE_names, x -> x is not null)")
        ).join(negativeTD, on=["targetId", "diseaseId"], how="left").withColumn(
        "PhaseT",
        F.when(F.col("stopReason") == "Negative", F.lit("yes")).otherwise(F.lit("no")),
        ).withColumn(
            "phase4",
            F.when(
                (F.col("maxClinPhase") == 4) & (F.col("PhaseT") == "no"), F.lit("yes")
            ).otherwise(F.lit("no")),
        ).withColumn(
            "phase>=3",
            F.when(
                (F.col("maxClinPhase") >= 3) & (F.col("PhaseT") == "no"), F.lit("yes")
            ).otherwise(F.lit("no")),
        ).withColumn(
            "phase>=2",
            F.when(
                (F.col("maxClinPhase") >= 2) & (F.col("PhaseT") == "no"), F.lit("yes")
            ).otherwise(F.lit("no")),
        ).withColumn(
            "phase>=1",
            F.when(
                (F.col("maxClinPhase") >= 1) & (F.col("PhaseT") == "no"), F.lit("yes")
            ).otherwise(F.lit("no")),
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
        .withColumn(
            "maxDoEArrayN",
            F.expr("aggregate(arrayN, 0, (acc, x) -> acc + IF(x = maxDoE, 1, 0))")
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
        )
        # .persist()
    )


#####
## no propag
#####
print("defining full analysis no propagation")


def full_analysis_noPropagation(
    doe_columns,assessment_all, analysisDatasources, analysis_chembl, negativeTD, diseaseTA,diagonal_lof,diagonal_gof
):
    conditions = [
    F.when(F.col(c) == F.col("maxDoE"), F.lit(c)).otherwise(F.lit(None)) for c in doe_columns
    ]
    return (
        analysis_nonPropagated(assessment_all, analysisDatasources)
        .join(
            analysis_chembl.filter(F.col("coherencyDiagonal") == "coherent").selectExpr(
                "targetId",
                "diseaseId",
                "maxClinPhase",
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
        ).withColumn(
            "arrayN", F.array(*[F.col(c) for c in doe_columns])
        ).withColumn(
            "maxDoE", F.array_max(F.col("arrayN"))
        ).withColumn("maxDoE_names", F.array(*conditions)
        ).withColumn("maxDoE_names", F.expr("filter(maxDoE_names, x -> x is not null)")
        ).join(negativeTD, on=["targetId", "diseaseId"], how="left").withColumn(
        "PhaseT",
        F.when(F.col("stopReason") == "Negative", F.lit("yes")).otherwise(F.lit("no")),
        ).withColumn(
            "phase4",
            F.when(
                (F.col("maxClinPhase") == 4) & (F.col("PhaseT") == "no"), F.lit("yes")
            ).otherwise(F.lit("no")),
        ).withColumn(
            "phase>=3",
            F.when(
                (F.col("maxClinPhase") >= 3) & (F.col("PhaseT") == "no"), F.lit("yes")
            ).otherwise(F.lit("no")),
        ).withColumn(
            "phase>=2",
            F.when(
                (F.col("maxClinPhase") >= 2) & (F.col("PhaseT") == "no"), F.lit("yes")
            ).otherwise(F.lit("no")),
        ).withColumn(
            "phase>=1",
            F.when(
                (F.col("maxClinPhase") >= 1) & (F.col("PhaseT") == "no"), F.lit("yes")
            ).otherwise(F.lit("no")),
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
        ).withColumn(
            "maxDoEArrayN",
            F.expr("aggregate(arrayN, 0, (acc, x) -> acc + IF(x = maxDoE, 1, 0))")
        ).withColumn(
            "NoneCellYes",
            F.when(F.col("LoF_protect_ch").isNotNull() & (F.array_contains(F.col("maxDoE_names"), F.lit("LoF_protect")))==True, F.lit('yes'))
            .when(F.col("GoF_protect_ch").isNotNull() & (F.array_contains(F.col("maxDoE_names"), F.lit("GoF_protect")))==True, F.lit('yes')
                ).otherwise(F.lit('no'))  # If the value is null, return null # Otherwise, check if name is in array
        ).withColumn(
            "NdiagonalYes",
            F.when(F.col("LoF_protect_ch").isNotNull() & 
                (F.size(F.array_intersect(F.col("maxDoE_names"), F.array([F.lit(x) for x in diagonal_lof]))) > 0),
                F.lit("yes")
            ).when(F.col("GoF_protect_ch").isNotNull() & 
                (F.size(F.array_intersect(F.col("maxDoE_names"), F.array([F.lit(x) for x in diagonal_gof]))) > 0),
                F.lit("yes")
            ).otherwise(F.lit('no'))
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
    "impc",
    "orphanet",
    "gene2phenotype",
    "gwas_credible_set",
    "cancer_gene_census",
]

datasource_list = [
    #"gene_burden",
    #"intogen",
    #"cancer_gene_census",
    #"eva",
    #"eva_somatic",
    "gwas_credible_set",
    #"impc",
    #"orphanet",
    #"gene2phenotype",
    #"WOcgc",
    #"wCgc",
    "somatic",
    "germline",
    "orpha_2_eva",
    "orpha_2",
    "orpha_2_eva_burden"
]

germline_list = [
    "gene_burden",
    "eva",
    "gwas_credible_set",
    "impc",
    "orphanet",
    "gene2phenotype",
]

### merge 'orphanet', 'gene2phenotype' and 'eva'

orpha_2_eva=[
    #"gene_burden",
    "eva",
    #"gwas_credible_set",
    #"impc",
    "orphanet",
    "gene2phenotype",
]

orpha_2=[
    #"gene_burden",
    "#eva",
    #"gwas_credible_set",
    #"impc",
    "orphanet",
    "gene2phenotype",
]

orpha_2_eva_burden=[
    "gene_burden",
    "#eva",
    #"gwas_credible_set",
    #"impc",
    "orphanet",
    "gene2phenotype",
]
somatic_list = ["intogen", "cancer_gene_census", "eva_somatic"]


# assessment = prueba_assessment.filter(F.col("datasourceId").isin(datasources_analysis))
def dataset_builder(assessment_all, value, analysis_chembl, negativeTD, diseaseTA):
    nonPropagated = full_analysis_noPropagation(
        doe_columns,assessment_all, value, analysis_chembl, negativeTD, diseaseTA,diagonal_lof,diagonal_gof
    )
    propagated = full_analysis_propagation(
        doe_columns,assessment_all, value, analysis_chembl, negativeTD, diseaseTA,diagonal_lof,diagonal_gof
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
    if value == "germline":
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
            assessment_all,
            germline_list,
            analysis_chembl,
            negativeTD,
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
            assessment_all,
            somatic_list,
            analysis_chembl,
            negativeTD,
            diseaseTA,
        )

    elif value == "orpha_2_eva":
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
            assessment_all,
            orpha_2_eva,
            analysis_chembl,
            negativeTD,
            diseaseTA,
        )

    elif value == "orpha_2":
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
            assessment_all,
            orpha_2,
            analysis_chembl,
            negativeTD,
            diseaseTA,
        )

    elif value == "orpha_2_eva_burden":
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
            assessment_all,
            orpha_2,
            analysis_chembl,
            negativeTD,
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
            dfs_dict_propag[f"df_{value}_Oncology_propag"]
        ) = dataset_builder(
            assessment_all, value, analysis_chembl, negativeTD, diseaseTA
        )






'''
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
            assessment_all, wocgc_list, analysis_chembl, negativeTD, diseaseTA
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
            assessment_all, wCgc_list, analysis_chembl, negativeTD, diseaseTA
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
            assessment_all,
            germline_list,
            analysis_chembl,
            negativeTD,
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
            assessment_all,
            somatic_list,
            analysis_chembl,
            negativeTD,
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
            dfs_dict_propag[f"df_{value}_Oncology_propag"]
        ) = dataset_builder(
            assessment_all, value, analysis_chembl, negativeTD, diseaseTA
        )
'''


def comparisons_df() -> list:
    """Return list of all comparisons to be used in the analysis"""
    comparisons = spark.createDataFrame(
        data=[
            ("hasGeneticEvidence", "byDatatype"),
            #("diagonalYes", "byDatatype"),
            #("oneCellYes", "byDatatype"),
            ("NdiagonalYes", "byDatatype"),
            ("NoneCellYes", "byDatatype"),
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
#for key, df in islice(dfs_dict.items(), 1): ## for debugging
for key, df in dfs_dict.items():
    df = df.persist()
    for row in aggSetups_original:
        aggregations_original(df, key, listado, *row, today_date)
    df.unpersist()
    print(key + " df unpersisted")

print("non propagated files wroten succesfully at", c)


print("starting with propagated aggregations at", c)
#for key, df in islice(dfs_dict_propag.items(), 1): ## for debugging
for key, df in dfs_dict_propag.items():
    df = df.persist()
    for row in aggSetups_original:
        aggregations_original(df, key, listado, *row, today_date)
    df.unpersist()
    print(key + " df unpersisted")

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

df.withColumn(
    "datasource",
    F.regexp_extract(F.col("group"), r"df_(.*?)_(All|Other|OtherNull|Oncology)_(propag|original)", 1)
).withColumn(
    "therArea",
    F.regexp_extract(F.col("group"), r"_(All|Other|OtherNull|Oncology)_", 1)
).withColumn(
    "type",
    F.regexp_extract(F.col("group"), r"_(propag|original)$", 1)
).toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_genEvidAnalysis_new_NoFileteredColocCaviar.csv"
)

print("dataframe written \n Analysis finished")
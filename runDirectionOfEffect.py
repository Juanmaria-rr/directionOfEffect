import pyspark
from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as F
from pyspark.sql import Window

from psutil import virtual_memory
from pyspark import SparkFiles
from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.functions import col


def detect_spark_memory_limit():
    """Spark does not automatically use all available memory on a machine. When working on large datasets, this may
    cause Java heap space errors, even though there is plenty of RAM available. To fix this, we detect the total amount
    of physical memory and allow Spark to use (almost) all of it."""
    mem_gib = virtual_memory().total >> 30
    return int(mem_gib * 0.9)


spark_mem_limit = detect_spark_memory_limit()
spark_conf = (
    SparkConf()
    .set("spark.driver.memory", f"{spark_mem_limit}g")
    .set("spark.executor.memory", f"{spark_mem_limit}g")
    .set("spark.driver.maxResultSize", "0")
    .set("spark.debug.maxToStringFields", "2000000000")
    .set("spark.sql.execution.arrow.maxRecordsPerBatch", "500000")
    .set("spark.sql.execution.arrow.pyspark.enabled", "true")
    .set("spark.ui.showConsoleProgress", "false")
)

spark = (
    SparkSession.builder.config(conf=spark_conf)
    .master("local[*]")
    .config("spark.driver.bindAddress", "127.0.0.1")
    .config("spark.driver.host", "localhost")
    .getOrCreate()
)


# 1# defining datasets and hardcoding some variables

### evideences datset and make an union between them:
otgenetics_evidence_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=ot_genetics_portal"
otgenetics = spark.read.parquet(otgenetics_evidence_path)
gene_burden_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=gene_burden"
gene_burden = spark.read.parquet(gene_burden_path)
eva_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=eva"
eva_germline = spark.read.parquet(eva_path)
eva_somatic_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=eva_somatic"
eva_somatic = spark.read.parquet(eva_somatic_path)
orphanet_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=orphanet"
orphanet = spark.read.parquet(orphanet_path)
g2p_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=gene2phenotype"
g2p = spark.read.parquet(g2p_path)
cgc_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=cancer_gene_census"
cgc = spark.read.parquet(cgc_path)
intogen_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=intogen"
intogen = spark.read.parquet(intogen_path)
impc_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=impc"
impc = spark.read.parquet(impc_path)
chembl_evidences = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=chembl/"
chembl = spark.read.parquet(chembl_evidences)

dfs = [
    otgenetics,
    gene_burden,
    eva_germline,
    eva_somatic,
    g2p,
    orphanet,
    cgc,
    intogen,
    impc,
    chembl,
]

allEvidences = dfs[0]
for df in dfs[1:]:
    allEvidences = allEvidences.unionByName(df, allowMissingColumns=True)

### external datasets:
target_path = (
    "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/targets/"
)
target = spark.read.parquet(target_path)
mecact_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/mechanismOfAction/"
mecact = spark.read.parquet(mecact_path)

### We manually annotated those studies using LoF or PTV variants - GeneBurden
burden_lof_path = "/Users/juanr/Desktop/directionOfEffect/geneBurden_20230117.csv"
burden_lof = spark.read.csv(burden_lof_path, header=True).withColumnRenamed(
    "statisticalMethodOverview", "stMethod"
)

var_filter_lof = [
    ### High impact variants https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
    "SO_0001589",  ## frameshit_variant
    "SO_0001587",  ## stop_gained
    "SO_0001574",  ## splice_acceptor_variant
    "SO_0001575",  ## splice_donor_variant
    "SO_0002012",  ## start_lost
    "SO_0001578",  ## stop_lost
    "SO_0001893",  ## transcript_ablation
]

gof = ["SO_0002053"]
lof = ["SO_0002054"]


## annotate TSG/oncogene/bivalent using 'hallmarks.attributes'
oncotsg_list = [
    "TSG",
    "oncogene",
    "Oncogene",
    "oncogene",
    "oncogene,TSG",
    "TSG,oncogene",
    "fusion,oncogene",
    "oncogene,fusion",
]

### define inhibitors and activators
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
]

activators = [
    "PARTIAL AGONIST",
    "ACTIVATOR",
    "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR",
    "AGONIST",
    "SEQUESTERING AGENT",
]

actionType = (
    mecact.select("chemblIds", "actionType", "mechanismOfAction", "targets")
    .select(
        F.explode_outer("chemblIds").alias("drugId2"),
        "actionType",
        "mechanismOfAction",
        "targets",
    )
    .select(
        F.explode_outer("targets").alias("targetId2"),
        "drugId2",
        "actionType",
        "mechanismOfAction",
    )
    .dropDuplicates()
)

oncolabel = (
    target.select(
        "id", "approvedSymbol", F.explode_outer(F.col("hallmarks.attributes"))
    )
    .select("id", "approvedSymbol", "col.description")
    .filter(F.col("description").isin(oncotsg_list))
    .groupBy("id", "approvedSymbol")
    .agg(F.collect_set("description").alias("description"))
    .withColumn("description_splited", F.concat_ws(",", F.col("description")))
    .withColumn(
        "TSorOncogene",
        F.when(
            (
                F.col("description_splited").rlike("ncogene")
                & F.col("description_splited").rlike("TSG")
            ),
            F.lit("bivalent"),
        )
        .when(F.col("description_splited").rlike("ncogene(\s|$)"), F.lit("oncogene"))
        .when(F.col("description_splited").rlike("TSG(\s|$)"), F.lit("TSG"))
        .otherwise(F.lit("noEvaluable")),
    )
    .withColumnRenamed("id", "target_id")
)

# 2# run the transformation of the evidences datasets used.


def directionOfEffect(
    allEvidences,
    oncolabel,
    burden_lof,
    actionType,
    var_filter_lof,
    gof,
    inhibitors,
    activators,
):

    doe = (
        allEvidences.withColumn(
            "beta", F.col("beta").cast("float")
        )  ## ot genetics & gene burden
        .withColumn(
            "OddsRatio", F.col("OddsRatio").cast("float")
        )  ## ot genetics & gene burden
        .withColumn(
            "clinicalSignificances", F.concat_ws(",", F.col("clinicalSignificances"))
        )  ### eva
        .join(oncolabel, oncolabel.target_id == F.col("targetId"), "left")  ###  cgc
        .join(
            burden_lof,
            burden_lof.stMethod == F.col("statisticalMethodOverview"),
            "left",
        )  ###  gene_burden
        .join(
            actionType,  ## chembl
            (actionType.drugId2 == F.col("drugId"))
            & (actionType.targetId2 == F.col("targetId")),
            "left",
        )
        ### variant Effect Column
        .withColumn(
            "variantEffect",
            F.when(
                F.col("datasourceId") == "ot_genetics_portal",
                F.when(
                    F.col("variantFunctionalConsequenceId").isNotNull(),
                    F.when(
                        F.col("variantFunctionalConsequenceFromQtlId").isNull(),
                        F.when(
                            F.col("variantFunctionalConsequenceId").isin(
                                var_filter_lof
                            ),
                            F.lit("LoF"),
                        )
                        .when(
                            F.col("variantFunctionalConsequenceId").isin(gof),
                            F.lit("GoF"),
                        )
                        .otherwise(F.lit("noEvaluable")),
                    )
                    ### variantFunctionalConsequenceFromQtlId
                    .when(
                        F.col("variantFunctionalConsequenceFromQtlId").isNotNull(),
                        F.when(
                            F.col("variantFunctionalConsequenceId").isin(
                                var_filter_lof
                            ),  ## when is a LoF variant
                            F.when(
                                F.col("variantFunctionalConsequenceFromQtlId")
                                == "SO_0002316",
                                F.lit("LoF"),
                            )
                            .when(
                                F.col("variantFunctionalConsequenceFromQtlId")
                                == "SO_0002315",
                                F.lit("conflict/noEvaluable"),
                            )
                            .otherwise(F.lit("LoF")),
                        ).when(
                            F.col("variantFunctionalConsequenceId").isin(var_filter_lof)
                            == False,  ## when is not a LoF, still can be a GoF
                            F.when(
                                F.col("variantFunctionalConsequenceId").isin(gof)
                                == False,  ##if not GoF
                                F.when(
                                    F.col("variantFunctionalConsequenceFromQtlId")
                                    == "SO_0002316",
                                    F.lit("LoF"),
                                )
                                .when(
                                    F.col("variantFunctionalConsequenceFromQtlId")
                                    == "SO_0002315",
                                    F.lit("GoF"),
                                )
                                .otherwise(F.lit("noEvaluable")),
                            ).when(
                                F.col("variantFunctionalConsequenceId").isin(
                                    gof
                                ),  ##if is GoF
                                F.when(
                                    F.col("variantFunctionalConsequenceFromQtlId")
                                    == "SO_0002316",
                                    F.lit("conflict/noEvaluable"),
                                ).when(
                                    F.col("variantFunctionalConsequenceFromQtlId")
                                    == "SO_0002315",
                                    F.lit("GoF"),
                                ),
                            ),
                        ),
                    ),
                ).when(
                    F.col("variantFunctionalConsequenceId").isNull(),
                    F.when(
                        F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002316",
                        F.lit("LoF"),
                    )
                    .when(
                        F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002315",
                        F.lit("GoF"),
                    )
                    .otherwise(F.lit("noEvaluable")),
                ),
            ).when(
                F.col("datasourceId") == "gene_burden",
                F.when(F.col("whatToDo") == "get", F.lit("LoF")).otherwise(
                    F.lit("noEvaluable")
                ),  ### son tambien no data las que tiene riesgo pero no se ensayan LoF o PT
            )
            #### Eva_germline
            .when(
                F.col("datasourceId") == "eva",
                #### .filter(F.col('variantFunctionalConsequenceId').isin(var_filter_lof))
                F.when(
                    F.col("variantFunctionalConsequenceId").isin(var_filter_lof),
                    F.lit("LoF"),
                ).otherwise(
                    F.lit("noEvaluable")
                ),  ### Son todas aquellas que tenen info pero no son LoF
            )
            #### Eva_somatic
            .when(
                F.col("datasourceId") == "eva_somatic",
                F.when(
                    F.col("variantFunctionalConsequenceId").isin(var_filter_lof),
                    F.lit("LoF"),
                ).otherwise(
                    F.lit("noEvaluable")
                ),  ### Son todas aquellas que tenen info pero no son patogenicas/protective  + LoF
            )
            #### G2P
            .when(
                F.col("datasourceId")
                == "gene2phenotype",  ### 6 types of variants [SO_0002318, SO_0002317, SO_0001622, SO_0002315, SO_0001566, SO_0002220]
                F.when(
                    F.col("variantFunctionalConsequenceId") == "SO_0002317",
                    F.lit("LoF"),
                )  ### absent gene product
                .when(
                    F.col("variantFunctionalConsequenceId") == "SO_0002315",
                    F.lit("GoF"),
                )  ### increased gene product level
                .otherwise(F.lit("noEvaluable")),
            )
            #### Orphanet
            .when(
                F.col("datasourceId") == "orphanet",
                F.when(
                    F.col("variantFunctionalConsequenceId") == "SO_0002054",
                    F.lit("LoF"),
                )  ### Loss of Function Variant
                .when(
                    F.col("variantFunctionalConsequenceId") == "SO_0002053",
                    F.lit("GoF"),
                )  ### Gain_of_Function Variant
                .otherwise(F.lit("noEvaluable")),
            )
            #### CGC
            .when(
                F.col("datasourceId") == "cancer_gene_census",
                F.when(F.col("TSorOncogene") == "oncogene", F.lit("GoF"))
                .when(F.col("TSorOncogene") == "TSG", F.lit("LoF"))
                .when(F.col("TSorOncogene") == "bivalent", F.lit("bivalent"))
                .otherwise("noEvaluable"),
            )
            #### intogen
            .when(
                F.col("datasourceId") == "intogen",
                F.when(
                    F.arrays_overlap(
                        F.col("mutatedSamples.functionalConsequenceId"),
                        F.array([F.lit(i) for i in (gof)]),
                    ),
                    F.lit("GoF"),
                )
                .when(
                    F.arrays_overlap(
                        F.col("mutatedSamples.functionalConsequenceId"),
                        F.array([F.lit(i) for i in (lof)]),
                    ),
                    F.lit("LoF"),
                )
                .otherwise(F.lit("noEvaluable")),
            )
            #### impc
            .when(
                F.col("datasourceId") == "impc",
                F.when(F.col("diseaseId").isNotNull(), F.lit("LoF")).otherwise(
                    F.lit("noEvaluable")
                ),
            )
            ### chembl
            .when(
                F.col("datasourceId") == "chembl",
                F.when(F.col("actionType").isin(inhibitors), F.lit("LoF"))
                .when(F.col("actionType").isin(activators), F.lit("GoF"))
                .otherwise(F.lit("noEvaluable")),
            ),
        )
        .withColumn(
            "directionOnTrait",
            ## ot genetics portal
            F.when(
                F.col("datasourceId")
                == "ot_genetics_portal",  ### the same for gene_burden
                F.when(
                    (F.col("beta").isNotNull()) & (F.col("OddsRatio").isNull()),
                    F.when(F.col("beta") > 0, F.lit("risk"))
                    .when(F.col("beta") < 0, F.lit("protect"))
                    .otherwise(F.lit("noEvaluable")),
                )
                .when(
                    (F.col("beta").isNull()) & (F.col("OddsRatio").isNotNull()),
                    F.when(F.col("OddsRatio") > 1, F.lit("risk"))
                    .when(F.col("OddsRatio") < 1, F.lit("protect"))
                    .otherwise(F.lit("noEvaluable")),
                )
                .when(
                    (F.col("beta").isNull()) & (F.col("OddsRatio").isNull()),
                    F.lit("noEvaluable"),
                )
                .when(
                    (F.col("beta").isNotNull()) & (F.col("OddsRatio").isNotNull()),
                    F.lit("conflict/noEvaluable"),
                ),
            ).when(
                F.col("datasourceId") == "gene_burden",
                F.when(
                    (F.col("beta").isNotNull()) & (F.col("OddsRatio").isNull()),
                    F.when(F.col("beta") > 0, F.lit("risk"))
                    .when(F.col("beta") < 0, F.lit("protect"))
                    .otherwise(F.lit("noEvaluable")),
                )
                .when(
                    (F.col("oddsRatio").isNotNull()) & (F.col("beta").isNull()),
                    F.when(F.col("oddsRatio") > 1, F.lit("risk"))
                    .when(F.col("oddsRatio") < 1, F.lit("protect"))
                    .otherwise(F.lit("noEvaluable")),
                )
                .when(
                    (F.col("beta").isNull()) & (F.col("oddsRatio").isNull()),
                    F.lit("noEvaluable"),
                )
                .when(
                    (F.col("beta").isNotNull()) & (F.col("oddsRatio").isNotNull()),
                    F.lit("conflict"),
                ),
            )
            ## Eva_germline
            .when(
                F.col("datasourceId") == "eva",  ### the same for eva_somatic
                F.when(
                    F.col("clinicalSignificances").rlike("(pathogenic)$"), F.lit("risk")
                )
                .when(
                    F.col("clinicalSignificances").contains("protect"), F.lit("protect")
                )
                .otherwise(
                    F.lit("noEvaluable")
                ),  ### Son todas aquellas que tenen info pero no son patogenicas/protective  + LoF
            )
            #### Eva_somatic
            .when(
                F.col("datasourceId") == "eva_somatic",
                F.when(
                    F.col("clinicalSignificances").rlike("(pathogenic)$"), F.lit("risk")
                )
                .when(
                    F.col("clinicalSignificances").contains("protect"), F.lit("protect")
                )
                .otherwise(
                    F.lit("noEvaluable")
                ),  ### Son todas aquellas que tenen info pero no son patogenicas/protective  + LoF
            )
            #### G2P
            .when(
                F.col("datasourceId") == "gene2phenotype",
                F.when(F.col("diseaseId").isNotNull(), F.lit("risk")).otherwise(
                    F.lit("noEvaluable")
                ),
            )
            #### Orphanet
            .when(
                F.col("datasourceId") == "orphanet",
                F.when(F.col("diseaseId").isNotNull(), F.lit("risk")).otherwise(
                    F.lit("noEvaluable")
                ),
            )
            #### CGC
            .when(
                F.col("datasourceId") == "cancer_gene_census",
                F.when(F.col("diseaseId").isNotNull(), F.lit("risk")).otherwise(
                    F.lit("noEvaluable")
                ),
            )
            #### intogen
            .when(
                F.col("datasourceId") == "intogen",
                F.when(F.col("diseaseId").isNotNull(), F.lit("risk")).otherwise(
                    F.lit("noEvaluable")
                ),
            )
            #### impc
            .when(
                F.col("datasourceId") == "impc",
                F.when(F.col("diseaseId").isNotNull(), F.lit("risk")).otherwise(
                    F.lit("noEvaluable")
                ),
            )
            ### chembl
            .when(
                F.col("datasourceId") == "chembl",
                F.when(F.col("diseaseId").isNotNull(), F.lit("protect")).otherwise(
                    F.lit("noEvaluable")
                ),
            ),
        )
        .withColumn(
            "homogenizedVersion",
            F.when(
                (F.col("variantEffect") == "LoF")
                & (F.col("directionOnTrait") == "risk"),
                F.lit("LoF_risk"),
            )
            .when(
                (F.col("variantEffect") == "LoF")
                & (F.col("directionOnTrait") == "protect"),
                F.lit("LoF_protect"),
            )
            .when(
                (F.col("variantEffect") == "GoF")
                & (F.col("directionOnTrait") == "risk"),
                F.lit("GoF_risk"),
            )
            .when(
                (F.col("variantEffect") == "GoF")
                & (F.col("directionOnTrait") == "protect"),
                F.lit("GoF_protect"),
            )
            .otherwise(F.lit("noEvaluable")),
        )
    )
    return doe

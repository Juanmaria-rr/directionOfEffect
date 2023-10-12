from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F

spark = SparkSession.builder.getOrCreate()

### read first passs of ETL results
evidences = spark.read.parquet(
    "gs://open-targets-pre-data-releases/ricardo/23.09/output/etl/parquet/evidence/"
)

# 1# defining datasets and hardcoding some variables
### external datasets:
target_path = "gs://open-targets-data-releases/23.09/output/etl/parquet/targets/"
target = spark.read.parquet(target_path)
mecact_path = (
    "gs://open-targets-data-releases/23.09/output/etl/parquet/mechanismOfAction/"
)
mecact = spark.read.parquet(mecact_path)

### We manually annotated those studies using LoF or PTV variants - GeneBurden
burden_lof_path = "gs://ot-team/jroldan/20230704_geneBurden_StudyInclusion.csv"
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
    "DISRUPTING AGENT",
]

activators = [
    "PARTIAL AGONIST",
    "ACTIVATOR",
    "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR",
    "AGONIST",
    "SEQUESTERING AGENT",
    "STABILISER",
]

actionType = (
    mecact.select(
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
    .groupBy("targetId2", "drugId2")
    .agg(
        F.collect_set("actionType").alias("actionType"),
    )
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


def new_directionOfEffect(
    evidences,
    oncolabel,
    burden_lof,
    actionType,
    var_filter_lof,
    gof,
    inhibitors,
    activators,
):
    
    all = evidences.filter(
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
    # Define a conditional column to represent the subset where datasourceId is "intogen"
    condition_col = F.when(F.col("datasourceId") == "intogen", 1).otherwise(0)
    # Define the Window specification partitioned by "targetId" and "diseaseId" and ordered by the condition column
    window_spec = Window.partitionBy("targetId", "diseaseId").orderBy(
        condition_col.desc()
    )

    doe = (
        all.withColumn(
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
        .withColumn("inhibitors_list", F.array([F.lit(i) for i in inhibitors]))
        .withColumn("activators_list", F.array([F.lit(i) for i in activators]))
        .withColumn("nullColumn", F.array(F.lit(None)))
        .withColumn(
            "intogenAnnot",
            F.size(
                F.flatten(
                    F.collect_set(
                        F.array_except(
                            F.col("mutatedSamples.functionalConsequenceId"),
                            F.col("nullColumn"),
                        )
                    ).over(window_spec)
                )
            ),
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
                    F.col("intogenAnnot")
                    == 1,  ## oncogene/tummor suppressor for a given trait
                    F.when(
                        F.arrays_overlap(
                            F.array_union(
                                F.col("mutatedSamples.functionalConsequenceId"),
                                F.array(),
                            ),
                            F.array([F.lit(i) for i in (gof)]),
                        ),
                        F.lit("GoF"),
                    ).when(
                        F.arrays_overlap(
                            F.array_union(
                                F.col("mutatedSamples.functionalConsequenceId"),
                                F.array(),
                            ),
                            F.array([F.lit(i) for i in (lof)]),
                        ),
                        F.lit("LoF"),
                    ),
                )
                .when(
                    F.col("intogenAnnot") > 1, F.lit("bivalentIntogen")
                )  ##oncogene & tumor suppressor for a given trait
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
                F.when(
                    F.size(
                        F.array_intersect(F.col("actionType"), F.col("inhibitors_list"))
                    )
                    >= 1,
                    F.lit("LoF"),
                )
                .when(
                    F.size(
                        F.array_intersect(F.col("actionType"), F.col("activators_list"))
                    )
                    >= 1,
                    F.lit("GoF"),
                )
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
            "homogenizedVersion_J",
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

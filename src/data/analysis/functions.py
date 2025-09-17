import pyspark.sql.functions as F


def discrepancifier(df):
    """
    detect discrepancies per row where there are the four
    DoE assessments using Null and isNotNull assessments
    """
    columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect", "noEvaluable"]

    for col in columns:
        if col not in df.columns:
            df = df.withColumn(col, F.lit(None)).persist()

    return df.withColumn(
        "coherencyDiagonal",
        F.when(
            (F.col("LoF_risk").isNull())
            & (F.col("LoF_protect").isNull())
            & (F.col("GoF_risk").isNull())
            & (F.col("GoF_protect").isNull())
            & (F.col("noEvaluable").isNull()),
            F.lit("noEvid"),
        )
        .when(
            (F.col("LoF_risk").isNull())
            & (F.col("LoF_protect").isNull())
            & (F.col("GoF_risk").isNull())
            & (F.col("GoF_protect").isNull())
            & (F.col("noEvaluable").isNotNull()),
            F.lit("EvidNotDoE"),
        )
        .when(
            (F.col("LoF_risk").isNotNull())
            | (F.col("LoF_protect").isNotNull())
            | (F.col("GoF_risk").isNotNull())
            | (F.col("GoF_protect").isNotNull()),
            F.when(
                ((F.col("GoF_risk").isNotNull()) & (F.col("LoF_risk").isNotNull())),
                F.lit("dispar"),
            )
            .when(
                ((F.col("LoF_protect").isNotNull()) & (F.col("LoF_risk").isNotNull())),
                F.lit("dispar"),
            )
            .when(
                ((F.col("GoF_protect").isNotNull()) & (F.col("GoF_risk").isNotNull())),
                F.lit("dispar"),
            )
            .when(
                (
                    (F.col("GoF_protect").isNotNull())
                    & (F.col("LoF_protect").isNotNull())
                ),
                F.lit("dispar"),
            )
            .otherwise(F.lit("coherent")),
        ),
    ).withColumn(
        "coherencyOneCell",
        F.when(
            (F.col("LoF_risk").isNull())
            & (F.col("LoF_protect").isNull())
            & (F.col("GoF_risk").isNull())
            & (F.col("GoF_protect").isNull())
            & (F.col("noEvaluable").isNull()),
            F.lit("noEvid"),
        )
        .when(
            (F.col("LoF_risk").isNull())
            & (F.col("LoF_protect").isNull())
            & (F.col("GoF_risk").isNull())
            & (F.col("GoF_protect").isNull())
            & (F.col("noEvaluable").isNotNull()),
            F.lit("EvidNotDoE"),
        )
        .when(
            (F.col("LoF_risk").isNotNull())
            | (F.col("LoF_protect").isNotNull())
            | (F.col("GoF_risk").isNotNull())
            | (F.col("GoF_protect").isNotNull()),
            F.when(
                F.col("LoF_risk").isNotNull()
                & (
                    (F.col("LoF_protect").isNull())
                    & (F.col("GoF_risk").isNull())
                    & (F.col("GoF_protect").isNull())
                ),
                F.lit("coherent"),
            )
            .when(
                F.col("GoF_risk").isNotNull()
                & (
                    (F.col("LoF_protect").isNull())
                    & (F.col("LoF_risk").isNull())
                    & (F.col("GoF_protect").isNull())
                ),
                F.lit("coherent"),
            )
            .when(
                F.col("LoF_protect").isNotNull()
                & (
                    (F.col("LoF_risk").isNull())
                    & (F.col("GoF_risk").isNull())
                    & (F.col("GoF_protect").isNull())
                ),
                F.lit("coherent"),
            )
            .when(
                F.col("GoF_protect").isNotNull()
                & (
                    (F.col("LoF_protect").isNull())
                    & (F.col("GoF_risk").isNull())
                    & (F.col("LoF_risk").isNull())
                ),
                F.lit("coherent"),
            )
            .otherwise(F.lit("dispar")),
        ),
    )


def relative_success(array1):
    from scipy.stats.contingency import relative_risk

    ### take numbers from array
    a, b = array1[0]
    c, d = array1[1]
    ####
    """
    Where zeros cause problems with computation of the relative risk or its standard error,
    0.5 is added to all cells (a, b, c, d) (Pagano & Gauvreau, 2000; Deeks & Higgins, 2010).
    """
    total_expo = a + b
    total_noExpo = c + d
    ### for cases when total_expo/total_noExpo = 0,
    ### we sum 1 to avoid errors an get at least 0 in the %
    if any(t == 0 for t in [total_expo, total_noExpo]):
        total_expo = total_expo + 1
        total_noExpo = total_noExpo + 1
        ### calculate relative success
        relative_success = relative_risk(a, total_expo, c, total_noExpo)
        ### calculate confidence intervals
        rs_ci = relative_risk(a, total_expo, c, total_noExpo).confidence_interval(
            confidence_level=0.95
        )
    else:

        ### calculate relative success
        relative_success = relative_risk(a, total_expo, c, total_noExpo)
        ### calculate confidence intervals
        rs_ci = relative_risk(a, total_expo, c, total_noExpo).confidence_interval(
            confidence_level=0.95
        )

    return relative_success.relative_risk, rs_ci


def spreadSheetFormatter(df):
    print("importing functions")
    from pyspark.sql.functions import format_number
    from pyspark.sql.types import (
        DoubleType,
    )
    print("imported functions")

    new_df = (
        df.withColumn(
            "significant",
            F.when(F.col("pValue") < 0.0001, F.lit("****"))
            .when((F.col("pValue") >= 0.0001) & (F.col("pValue") < 0.001), F.lit("***"))
            .when((F.col("pValue") >= 0.001) & (F.col("pValue") < 0.01), F.lit("**"))
            .when((F.col("pValue") >= 0.01) & (F.col("pValue") < 0.05), F.lit("*"))
            .when(F.col("pValue") >= 0.05, F.lit("ns")),
        )
        .withColumn(
            "writeFigure",
            F.concat(
                F.round(F.col("oddsRatio"), 2),
                F.lit(" "),
                F.lit("("),
                F.round(F.col("lowerInterval"), 2),
                F.lit("-"),
                F.round(F.col("upperInterval"), 2),
                F.lit(")"),
            ),
        )
        .withColumn("value_1", F.col("values").getItem(0).getItem(0))
        .withColumn("value_2", F.col("values").getItem(0).getItem(1))
        .withColumn("value_3", F.col("values").getItem(1).getItem(0))
        .withColumn("value_4", F.col("values").getItem(1).getItem(1))
        .withColumn("numerator", (F.col("value_1") + F.col("value_2")).cast("int"))
        .withColumn("denominator", (F.col("value_3") + F.col("value_4")).cast("int"))
        .withColumn("pValue", F.col("pValue").cast(DoubleType()))
        .withColumn(
            "pValue_formatted",
            F.when(
                F.col("pValue") < 0.0001, F.col("pValue").cast("string")
            )  # Check if value matches scientific notation range
            # .when(F.col("pValue") < 0.05, F.format_number(F.col("pValue"), 2))
            .otherwise(F.format_number(F.col("pValue"), 4)),
        )
        .withColumn(
            "numDen", F.concat_ws("/", F.col("numerator"), F.col("denominator"))
        )
    )
    return new_df


#####


def temporary_directionOfEffect(path, datasource_filter):
    from pyspark.sql import SparkSession, Window
    from pyspark.sql import functions as F

    #spark = SparkSession.builder.getOrCreate()
    """
    Function to develop DoE assessment from OT evidences files, creating two columns:
    direction on target and direction on trait
    Args:
        path (str): Base path for evidence and target data.
        datasource_filter (list): List of values to filter the datasourceId column.
    """
    evidences = (
        spark.read.parquet(f"{path}evidence/")
        .filter(F.col("datasourceId").isin(datasource_filter))
        .persist()
    ).withColumn('variantFunctionalConsequenceFromQtlId', F.lit(None)) ### create mock column to avoid fixing temporaryDirectionOfEffect

    # Create the paths using the version variable
    target_path = f"{path}target/"
    mecact_path = f"{path}drug_mechanism_of_action/" #  mechanismOfAction == old version

    target = spark.read.parquet(target_path)
    mecact = spark.read.parquet(mecact_path)

    # Define variant filters and other lists
    var_filter_lof = [
        "SO_0001589",
        "SO_0001587",
        "SO_0001574",
        "SO_0001575",
        "SO_0002012",
        "SO_0001578",
        "SO_0001893",
    ]
    gof = ["SO_0002053"]
    lof = ["SO_0002054"]

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
        .agg(F.collect_set("actionType").alias("actionType"))
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
            .when(
                F.col("description_splited").rlike("ncogene(\s|$)"), F.lit("oncogene")
            )
            .when(F.col("description_splited").rlike("TSG(\s|$)"), F.lit("TSG"))
            .otherwise(F.lit("noEvaluable")),
        )
        .withColumnRenamed("id", "target_id")
    )
    windowSpec = Window.partitionBy("targetId", "diseaseId")

    assessment = (
        evidences.withColumn(
            "beta", F.col("beta").cast("double")
        )  ## ot genetics & gene burden
        .withColumn(
            "OddsRatio", F.col("OddsRatio").cast("double")
        )  ## ot genetics & gene burden
        .withColumn(
            "clinicalSignificances", F.concat_ws(",", F.col("clinicalSignificances"))
        )  ### eva
        .join(oncolabel, oncolabel.target_id == F.col("targetId"), "left")  ###  cgc
        .join(
            actionType,  ## chembl
            (actionType.drugId2 == F.col("drugId"))
            & (actionType.targetId2 == F.col("targetId")),
            "left",
        )
        .withColumn("inhibitors_list", F.array([F.lit(i) for i in inhibitors]))
        .withColumn("activators_list", F.array([F.lit(i) for i in activators]))
        .withColumn(
            "intogen_function",
            F.when(
                F.arrays_overlap(
                    F.col("mutatedSamples.functionalConsequenceId"),
                    F.array([F.lit(i) for i in (gof)]),
                ),
                F.lit("GoF"),
            ).when(
                F.arrays_overlap(
                    F.col("mutatedSamples.functionalConsequenceId"),
                    F.array([F.lit(i) for i in (lof)]),
                ),
                F.lit("LoF"),
            ),
            # .otherwise("nodata"),
        )
        .withColumn(
            "intogenAnnot",
            F.size(F.collect_set(F.col("intogen_function")).over(windowSpec)),
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
                F.when(F.col("targetId").isNotNull(), F.lit("LoF")).otherwise(
                    F.lit("noEvaluable")
                ),
            )
            #### Eva_germline
            .when(
                F.col("datasourceId") == "eva",
                F.when(
                    F.col("variantFunctionalConsequenceId").isin(var_filter_lof),
                    F.lit("LoF"),
                ).otherwise(F.lit("noEvaluable")),
            )
            #### Eva_somatic
            .when(
                F.col("datasourceId") == "eva_somatic",
                F.when(
                    F.col("variantFunctionalConsequenceId").isin(var_filter_lof),
                    F.lit("LoF"),
                ).otherwise(F.lit("noEvaluable")),
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
                            F.col("mutatedSamples.functionalConsequenceId"),
                            F.array([F.lit(i) for i in (gof)]),
                        ),
                        F.lit("GoF"),
                    ).when(
                        F.arrays_overlap(
                            F.col("mutatedSamples.functionalConsequenceId"),
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
                F.col("datasourceId") == "ot_genetics_portal",
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
                F.col("datasourceId") == "eva",
                F.when(
                    F.col("clinicalSignificances").rlike("(pathogenic)$"), F.lit("risk")
                )
                .when(
                    F.col("clinicalSignificances").contains("protect"), F.lit("protect")
                )
                .otherwise(F.lit("noEvaluable")),
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
                .otherwise(F.lit("noEvaluable")),
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
            "homogenized",
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
    ).persist()

    return assessment, evidences, actionType, oncolabel


""" 
    usage of previous function:
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
"""


## https://stackoverflow.com/questions/45629781/drop-if-all-entries-in-a-spark-dataframes-specific-column-is-null
## drop columns with all values = Null
def drop_fully_null_columns(df, but_keep_these=[]):
    """Drops DataFrame columns that are fully null
    (i.e. the maximum value is null)

    Arguments:
        df {spark DataFrame} -- spark dataframe
        but_keep_these {list} -- list of columns to keep without checking for nulls

    Returns:
        spark DataFrame -- dataframe with fully null columns removed
    """

    # skip checking some columns
    cols_to_check = [col for col in df.columns if col not in but_keep_these]
    if len(cols_to_check) > 0:
        # drop columns for which the max is None
        rows_with_data = (
            df.select(*cols_to_check)
            .groupby()
            .agg(*[F.max(c).alias(c) for c in cols_to_check])
            .take(1)[0]
        )
        cols_to_drop = [
            c for c, const in rows_with_data.asDict().items() if const == None
        ]
        new_df = df.drop(*cols_to_drop)

        return new_df
    else:
        return df


def convertTuple(tup):
    st = ",".join(map(str, tup))
    return st


def relative_success(array1):
    from scipy.stats.contingency import relative_risk

    ### take numbers from array
    a, b = array1[0]
    c, d = array1[1]
    ####
    """
    Where zeros cause problems with computation of the relative risk or its standard error,
    0.5 is added to all cells (a, b, c, d) (Pagano & Gauvreau, 2000; Deeks & Higgins, 2010).
    """
    total_expo = a + b
    total_noExpo = c + d
    ### for cases when total_expo/total_noExpo = 0,
    ### we sum 1 to avoid errors an get at least 0 in the %
    if any(t == 0 for t in [total_expo, total_noExpo]):
        total_expo = total_expo + 1
        total_noExpo = total_noExpo + 1
        ### calculate relative success
        relative_success = relative_risk(a, total_expo, c, total_noExpo)
        ### calculate confidence intervals
        rs_ci = relative_risk(a, total_expo, c, total_noExpo).confidence_interval(
            confidence_level=0.95
        )
    else:

        ### calculate relative success
        relative_success = relative_risk(a, total_expo, c, total_noExpo)
        ### calculate confidence intervals
        rs_ci = relative_risk(a, total_expo, c, total_noExpo).confidence_interval(
            confidence_level=0.95
        )

    return relative_success.relative_risk, rs_ci


from pyspark.sql import SparkSession
import pyspark.sql.functions as F
from pyspark.sql.types import StructType, StructField, StringType, IntegerType
import pandas as pd

path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"
spark = SparkSession.builder.getOrCreate()


def build_gwasResolvedColoc(path):

    #### Now load sources of data to generate credible_set_OT_genetics evidences and associations.

    diseases = spark.read.parquet(f"{path}diseases/")

    credibleEvidence = spark.read.parquet(f"{path}evidence").filter(
        F.col("datasourceId").isin(["gwas_credible_sets"])
    )
    credible = spark.read.parquet(f"{path}credibleSet")

    index = spark.read.parquet(f"{path}gwasIndex")

    new = spark.read.parquet(f"{path}colocalisation/coloc")

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
                F.regexp_replace(
                    F.col("studyId"), r"(_ENS.*|_ILMN.*|_X.*|_[0-9]+:.*)", ""
                ),
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
                F.explode_outer(
                    F.concat(F.array(F.col("diseaseId")), F.col("parents"))
                ),
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
    return gwasResolvedColoc


def build_gwasResolvedColoc_noPropag(path):

    #### Now load sources of data to generate credible_set_OT_genetics evidences and associations.

    #diseases = spark.read.parquet(f"{path}diseases/")

    credibleEvidence = spark.read.parquet(f"{path}evidence").filter(
        F.col("datasourceId").isin(["gwas_credible_sets"])
    )
    credible = spark.read.parquet(f"{path}credibleSet")

    index = spark.read.parquet(f"{path}gwasIndex")

    new = spark.read.parquet(f"{path}colocalisation/coloc")

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
                F.regexp_replace(
                    F.col("studyId"), r"(_ENS.*|_ILMN.*|_X.*|_[0-9]+:.*)", ""
                ),
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

    gwasResolvedColoc_noPropag = (
        (
            newColoc.filter(F.col("rightStudyType") != "gwas")
            .withColumnRenamed("geneId", "targetId")
            .join(
                gwasComplete.withColumnRenamed("studyLocusId", "leftStudyLocusId"),
                on=["leftStudyLocusId", "targetId"],
                how="right",
            )
            #.join(  ### propagated using parent terms
            #    diseases.selectExpr(
            #        "id as diseaseId", "name", "parents", "therapeuticAreas"
            #    ),
            #    on="diseaseId",
            #    how="left",
            #)
            #.withColumn(
            #    "diseaseId",
            #    F.explode_outer(
            #        F.concat(F.array(F.col("diseaseId")), F.col("parents"))
            #    ),
            #)
            #.drop("parents", "oldDiseaseId")
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
    return gwasResolvedColoc_noPropag,gwasComplete

def buildColocData(all_coloc,credible,index):
    return (
    all_coloc.join(
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
            "pValueExponent as qtlPValueExponent",
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
    ))

def gwasDataset(evidences,credible):
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
    )
    print("loaded gwasComplete")
    return gwasComplete

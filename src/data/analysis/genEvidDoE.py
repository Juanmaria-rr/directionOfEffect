"""
This scripts run Odds ratio analysis for DoE and 
genetic information on drug clinical success

"""

from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
import matplotlib.pyplot as plt
import seaborn as sns
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
from datetime import datetime


spark = SparkSession.builder.getOrCreate()
c = datetime.now()
today_date = str(date.today())
print("spark session created at", c)

print("Analysis started on " + today_date + " at ", c)
"""
#coloc = spark.read.parquet(
#    "gs://genetics-portal-dev-data/22.09.1/outputs/v2d_coloc"
#).filter(F.col("right_type") != "gwas")
"""

#### make the dataset from stopped clin trials
### read supplementary table 9
""" ### just showing how i did the dataset
st9 = spark.read.csv("/Users/juanr/Downloads/ST9.csv", sep=",", header=True)
st9.filter(
    (F.col("clinicalStatus").isin(["Terminated", "Withdrawn", "Suspended"]))
    & (F.col("prediction") == "Negative")
).groupBy(
    "targetId", "diseaseId", "clinicalStatus", "prediction"
).count().toPandas().to_csv(
    "targetDiseaseStoppedNegative.csv"
)
"""
### target-diseases terminated&withdrawal in clin trials
terminated = spark.read.csv(
    "gs://ot-team/jroldan/analysis/targetDiseaseStoppedNegative.csv",
    sep=",",
    header=True,
).drop("_c0", "Withdrawn")

evidences = (
    spark.read.parquet(
        "gs://open-targets-data-releases/24.06/output/etl/parquet/evidence"
    )
    .filter(
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
    .persist()
)
ot_genetics = evidences.filter(F.col("datasourceId") == "ot_genetics_portal")

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


# 1# Make a list of variant of interest (Sequence ontology terms) to subset data of interest.

### Bear in mind that SO works with ontology structure as: SO:XXXXXX, but databases has the SO as: SO_XXXXXX

var_filter_lof = [
    ### High impact variants https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
    "SO_0001589",  ## frameshit_variant
    "SO_0001587",  ## stop_gained
    "SO_0001574",  ## splice_acceptor_variant
    "SO_0001575",  ## splice_donor_variant
    "SO_0002012",  ## start_lost
    "SO_0001578",  ## stop_lost
    "SO_0001893",  ## transcript_ablation
    # "SO:0001889", ## transcript_amplification ## the Only HIGH impact that increase protein.
]

gof = ["SO_0002053"]
lof = ["SO_0002054"]

print("loading sources")

## Building Sequence Ontology
so_path = "gs://ot-team/jroldan/sequenceOntology_20221118.csv"
so_ontology = spark.read.csv(so_path, header=True)
building = (
    so_ontology.select(F.col("Accession"), F.col("Parents"))
    .withColumn("Parentalind", F.split(F.col("Parents"), ","))
    .withColumn("Parentalind", F.explode_outer("Parentalind"))
    .groupBy("Parentalind")
    .agg(F.collect_list(F.col("Accession")).alias("childrens"))
    .join(so_ontology, F.col("Parentalind") == so_ontology.Accession, "right")
)

## others
target_path = "gs://open-targets-data-releases/24.06/output/etl/parquet/targets/"
target = spark.read.parquet(target_path)
disease_path = "gs://open-targets-data-releases/24.06/output/etl/parquet/diseases/"
diseases = spark.read.parquet(disease_path)
dis_name = diseases.select("id", "name")
indication_path = "gs://open-targets-data-releases/24.06/output/etl/parquet/indication/"
indication = spark.read.parquet(indication_path)
mecact_path = (
    "gs://open-targets-data-releases/24.06/output/etl/parquet/mechanismOfAction/"
)
mecact = spark.read.parquet(mecact_path)

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

#### rlike('('+Keywords+')(\s|$)'
### on 03.07.2023 we add the categories:
# DISRUPTING AGENT - inhibitor
# STABILISER - activator

### Hacer el join del actionType con el chembl para sacar los mecanismos de accion.
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
    "DISRUPTING AGENT",  ## added new on 03.07.2023
]

activators = [
    "PARTIAL AGONIST",
    "ACTIVATOR",
    "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR",
    "AGONIST",
    "SEQUESTERING AGENT",
    "STABILISER",  ## added new on 03.07.2023
]

columnas = ["activator", "inhibitor"]
both = activators + inhibitors

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

windowSpec = Window.partitionBy("targetId", "diseaseId")

#### version all gene burden
prueba_assessment = (
    all.withColumn("beta", F.col("beta").cast("double"))  ## ot genetics & gene burden
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
                        F.col("variantFunctionalConsequenceId").isin(var_filter_lof),
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
                F.size(F.array_intersect(F.col("actionType"), F.col("inhibitors_list")))
                >= 1,
                F.lit("LoF"),
            )
            .when(
                F.size(F.array_intersect(F.col("actionType"), F.col("activators_list")))
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
            F.col("datasourceId") == "ot_genetics_portal",  ### the same for gene_burden
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
            F.when(F.col("clinicalSignificances").rlike("(pathogenic)$"), F.lit("risk"))
            .when(F.col("clinicalSignificances").contains("protect"), F.lit("protect"))
            .otherwise(
                F.lit("noEvaluable")
            ),  ### Son todas aquellas que tenen info pero no son patogenicas/protective  + LoF
        )
        #### Eva_somatic
        .when(
            F.col("datasourceId") == "eva_somatic",
            F.when(F.col("clinicalSignificances").rlike("(pathogenic)$"), F.lit("risk"))
            .when(F.col("clinicalSignificances").contains("protect"), F.lit("protect"))
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
        "homogenized",
        F.when(
            (F.col("variantEffect") == "LoF") & (F.col("directionOnTrait") == "risk"),
            F.lit("LoF_risk"),
        )
        .when(
            (F.col("variantEffect") == "LoF")
            & (F.col("directionOnTrait") == "protect"),
            F.lit("LoF_protect"),
        )
        .when(
            (F.col("variantEffect") == "GoF") & (F.col("directionOnTrait") == "risk"),
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


###
# dataset for coherencies and DoE
###

columns_chembl = ["LoF_protect", "GoF_protect"]
columns_dataset = ["LoF_protect", "GoF_protect", "LoF_risk", "GoF_risk", "evidenceDif"]
columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect"]
terms = ["noEvaluable", "bivalent_risk", "null", "dispar"]

genEvidDataset = (
    prueba_assessment.filter(F.col("datasourceId") != "chembl")  #### checked 31.05.2023
    .groupBy("targetId", "diseaseId")
    .agg(F.count("targetId").alias("Nr_evidences"))
    .select("targetId", "diseaseId", "Nr_evidences")
    .withColumn("geneticEvidence", F.lit("hasGeneticEvidence"))
)

coherency_toAssess_others_datasource = (  #### checked 31.05.2023
    prueba_assessment.filter(
        (F.col("homogenized").isin(columns)) & (F.col("datasourceId") != "chembl")
    )
    .groupBy("targetId", "diseaseId")
    .agg(F.collect_set("datasourceId").alias("datasourceIds"))
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

#### give us propagation of diseases and list of therapeutic areas associated
diseases2 = diseases.select("id", "parents").withColumn(
    "diseaseIdPropagated",
    F.explode_outer(F.concat(F.array(F.col("id")), F.col("parents"))),
)

chembl_trials = (
    prueba_assessment.filter((F.col("datasourceId").isin(["chembl"])))
    .groupBy("targetId", "diseaseId")
    .agg(F.max(F.col("clinicalPhase")).alias("maxClinPhase"))
)


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


terminated_array = (
    terminated.groupBy("targetId", "diseaseId")
    .agg(F.collect_set("clinicalStatus").alias("clinicalStatus"))
    .withColumn("prediction", F.when(F.col("clinicalStatus").isNotNull(), F.lit("yes")))
)


def analysis_nonPropagated(prueba_assessment, analysisDatasources):
    return discrepancifier(
        prueba_assessment.filter(F.col("datasourceId").isin(analysisDatasources))
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
        .persist()
    )


def analysis_propagated(prueba_assessment, analysisDatasources):
    return discrepancifier(
        prueba_assessment.filter(F.col("datasourceId").isin(analysisDatasources))
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
        .persist()
    )


chembl_ds = ["chembl"]


def analysis_drugs(prueba_assessment, chembl_ds):
    return discrepancifier(
        prueba_assessment.filter((F.col("datasourceId").isin(chembl_ds)))
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


analysis_chembl = analysis_drugs(prueba_assessment, chembl_ds)

#######
## include here the analysis
#######

analysisDatasources = []


def full_analysis_propagation(
    prueba_assessment, analysisDatasources, analysis_chembl, terminated_array, diseaseTA
):
    return (
        analysis_propagated(prueba_assessment, analysisDatasources)
        .join(
            analysis_chembl.selectExpr(
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
        .join(terminated_array, on=["targetId", "diseaseId"], how="left")
        .withColumn(
            "PhaseT",
            F.when(F.col("prediction") == "yes", F.lit("yes")).otherwise(F.lit("no")),
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
        .persist()
    )


#####
## no propag
#####
def full_analysis_noPropagation(
    prueba_assessment, analysisDatasources, analysis_chembl, terminated_array, diseaseTA
):
    return (
        analysis_nonPropagated(prueba_assessment, analysisDatasources)
        .join(
            analysis_chembl.selectExpr(
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
        .join(terminated_array, on=["targetId", "diseaseId"], how="left")
        .withColumn(
            "PhaseT",
            F.when(F.col("prediction") == "yes", F.lit("yes")).otherwise(F.lit("no")),
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
        .persist()
    )


########
## review until here
########

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
    "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
]
datasource_list = [
    "gene_burden",
    "intogen",
    "cancer_gene_census",
    "eva",
    "eva_somatic",
    "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
    # "chembl",
    "WOcgc",
    "somatic",
    "germline",
]

germline_list = [
    "gene_burden",
    "eva",
    "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
]

somatic_list = ["intogen", "cancer_gene_census", "eva_somatic"]

# assessment = prueba_assessment.filter(F.col("datasourceId").isin(datasources_analysis))


def dataset_builder(
    prueba_assessment, value, analysis_chembl, terminated_array, diseaseTA
):
    nonPropagated = full_analysis_noPropagation(
        prueba_assessment, value, analysis_chembl, terminated_array, diseaseTA
    )
    propagated = full_analysis_propagation(
        prueba_assessment, value, analysis_chembl, terminated_array, diseaseTA
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
            prueba_assessment, wocgc_list, analysis_chembl, terminated_array, diseaseTA
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
            prueba_assessment,
            germline_list,
            analysis_chembl,
            terminated_array,
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
            prueba_assessment,
            somatic_list,
            analysis_chembl,
            terminated_array,
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
            dfs_dict_propag[f"df_{value}_Oncology_propag"],
        ) = dataset_builder(
            prueba_assessment, value, analysis_chembl, terminated_array, diseaseTA
        )


def comparisons_df() -> list:
    """Return list of all comparisons to be used in the analysis"""
    comparisons = spark.createDataFrame(
        data=[
            ("hasGeneticEvidence", "byDatatype"),
            ("diagonalYes", "byDatatype"),
            ("oneCellYes", "byDatatype"),
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


c = datetime.now()

print("start doing aggregations and writing")
today_date = str(date.today())
aggSetups_original = comparisons_df()
listado = []

print("starting with non-propagated aggregations at", c)

for key, df in dfs_dict.items():
    df = df.persist()
    for row in aggSetups_original:
        aggregations_original(df, key, listado, *row, today_date)
    df.unpersist()
    print(key + " df unpersisted")

print("non propagated files wroten succesfully at", c)

print("starting with non-propagated aggregations at", c)
for key, df in dfs_dict_propag.items():
    df = df.persist()
    for row in aggSetups_original:
        aggregations_original(df, key, listado, *row, today_date)
    df.unpersist()
    print(key + " df unpersisted")

print("propagated files wroten succesfully at", c)

##### read files and make spreadsheet

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


# Define the lists of possible substrings

comparison_opt = ["oneCellYes", "diagonalYes", "hasGeneticEvidence"]
phase_opt = [r"Phase>=\d+", r"Phase\d+", r"PhaseT"]
propag_opt = ["propag", "original"]
# Join the lists into a single pattern for each part
ds_pattern = "|".join(datasource_list)
predict_pattern = "|".join(comparison_opt)
phase_pattern = "|".join(phase_opt)
propag_pattern = "|".join(propag_opt)

# Construct the full pattern dynamically
pattern = re.compile(
    rf"df_({ds_pattern})_(.*?)_({propag_pattern})/({predict_pattern})_({phase_pattern})\.parquet"
)
"""
    "2024-07-11_analysis/df_gene_burden_Oncology_original/oneCellYes_Phase>=1.parquet",
"""
result = []
result_st = []
result_ci = []
array2 = []

# Initialize an empty list to store the results
results = []

# Iterate over the sample strings and extract the desired substrings
for path in listado:

    match = pattern.search(path)
    array1 = np.delete(
        (
            spark.read.parquet(path)
            .join(full_data, on=["prediction", "comparison"], how="outer")
            .groupBy("comparison")
            .pivot("prediction")
            .agg(F.first("a"))
            .sort(F.col("comparison").desc())
            .select("comparison", "yes", "no")
            .fillna(0)
            .toPandas()
            .to_numpy()
        ),
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

    if match:
        datasource = match.group(1)
        therapyArea = match.group(2)
        propag = match.group(3)
        comparison = match.group(4)
        phase = match.group(5)

    else:
        datasource = "unknown"
        therapyArea = "unknown"
        propag = "unknown"
        comparison = "unknown"
        phase = "unknown"

    results.append(
        [
            datasource,
            therapyArea,
            propag,
            comparison,
            phase,
            round(float(resX.split(",")[0]), 2),
            round(float(resX.split(",")[1]), 2),
            round(float(resx_CI.split(",")[0]), 2),
            round(float(resx_CI.split(",")[1]), 2),
            str(total),
            res_npPhaseX,
            round(float(rs_result), 2),
            round(float(rs_ci[0]), 2),
            round(float(rs_ci[1]), 2),
            path,
        ]
    )
    print(path)
# Convert the results to a pandas DataFrame
df = pd.DataFrame(
    results,
    columns=[
        "datasource",
        "therapyArea",
        "propag",
        "comparison",
        "phase",
        "oddsRatio",
        "pValue",
        "lowerInterval",
        "upperInterval",
        "total",
        "values",
        "relSuccess",
        "rsLower",
        "rsUpper",
        "path",
    ],
)

print("writing csv file")
file_name = f"gs://ot-team/jroldan/analysis/{today_date}_analysis.csv"
df.to_csv(file_name)
print("csv succesfully created")

""" Build datasets with DoE evidence and assoc numbers"""

from DoEAssessment import directionOfEffect
from pyspark.sql import SparkSession
import pyspark.sql.functions as F
from pyspark.ml.feature import QuantileDiscretizer
from functions import discrepancifier
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from pyspark.sql.types import StructType, StructField, StringType, IntegerType
from pyspark.sql.types import (
    StructType,
    StructField,
    ArrayType,
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


### target-diseases terminated&withdrawal in clin trials
terminated = spark.read.csv(
    "gs://ot-team/jroldan/analysis/targetDiseaseStoppedNegative.csv",
    sep=",",
    header=True,
).drop("_c0", "Withdrawn")

path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"

#### Now load sources of data to generate credible_set_OT_genetics evidences and associations.
print("loading files from newest path")
target = spark.read.parquet(f"{path}targets/")

diseases = spark.read.parquet(f"{path}diseases/")

evidences = spark.read.parquet(f"{path}evidence")

credible = spark.read.parquet(f"{path}credibleSet")

index = spark.read.parquet(f"{path}gwasIndex")

new = spark.read.parquet(f"{path}colocalisation/coloc")

variantIndex = spark.read.parquet(f"{path}variantIndex")

biosample = spark.read.parquet(f"{path}biosample")


### Take just gwas_credible_sets
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
    .persist()
)


## bring studyId, variantId, beta from Gwas and pValue
gwasComplete = filtered_df.join(
    credible.selectExpr(
        "studyLocusId", "studyId", "variantId", "beta as betaGwas", "pValueExponent"
    ),
    on="studyLocusId",
    how="left",
)

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

#### take the direction from the lowest p value
window_spec = Window.partitionBy("targetId", "diseaseId").orderBy(
    F.col("pValueExponent").asc()
)
gwasCredibleAssoc = (
    gwasResolvedColoc.withColumn("homogenized", F.first("colocDoE").over(window_spec))
    .select("targetId", "diseaseId", "homogenized")
    .withColumn(
        "homogenized",
        F.when(F.col("homogenized").isNull(), F.lit("noEvaluable")).otherwise(
            F.col("homogenized")
        ),
    )
)

######
## direction of effect
######

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

## others
target_path = f"{path}targets/"
target = spark.read.parquet(target_path)
disease_path = f"{path}diseases/"
diseases = spark.read.parquet(disease_path)
dis_name = diseases.select("id", "name")
indication_path = f"{path}indication/"
indication = spark.read.parquet(indication_path)
mecact_path = f"{path}mechanismOfAction/"
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
            #"ot_genetics_portal",
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
diseases2 = diseases.select("id", "parents").withColumn(
    "diseaseIdPropagated",
    F.explode_outer(F.concat(F.array(F.col("id")), F.col("parents"))),
)

chembl_trials = (
    prueba_assessment.filter((F.col("datasourceId").isin(["chembl"])))
    .groupBy("targetId", "diseaseId")
    .agg(F.max(F.col("clinicalPhase")).alias("maxClinPhase"))
)


terminated_array = (
    terminated.groupBy("targetId", "diseaseId")
    .agg(F.collect_set("clinicalStatus").alias("clinicalStatus"))
    .withColumn("prediction", F.when(F.col("clinicalStatus").isNotNull(), F.lit("yes")))
)



replacement_dict = {
    "gene_burden": "GeneBurden",
    "chembl": "ChEMBL",
    "intogen": "Intogen",
    "orphanet": "Orphanet",
    "cancer_gene_census": "CancerGeneCensus",
    "eva": "EvaGermline",
    "gene2phenotype": "Gene2Phenotype",
    "eva_somatic": "EvaSomatic",
    # "ot_genetics_portal": "OtGenetics",
    "impc": "IMPC",
    "gwas_credible_set": "OtGenetics",
}

assessment = prueba_assessment.filter(F.col("datasourceId")!="ot_genetics_portal").unionByName(
    gwasCredibleAssoc.withColumn("datasourceId", F.lit("gwas_credible_set")), allowMissingColumns=True
).withColumn("directionOnTrait", 
            F.when( 
                (F.col("datasourceId")=="gwas_credible_set") 
                & (F.col("homogenized").contains("risk")), 
                    F.lit("risk")
            ).when((F.col("datasourceId")=="gwas_credible_set") 
                & (F.col("homogenized").contains("protect")), 
                    F.lit("protect")
            ).otherwise(F.col("directionOnTrait"))
        ).withColumn("variantEffect", 
            F.when( 
                (F.col("datasourceId")=="gwas_credible_set") 
                & (F.col("homogenized").contains("GoF")), 
                    F.lit("GoF")
            ).when((F.col("datasourceId")=="gwas_credible_set") 
                & (F.col("homogenized").contains("LoF")), 
                    F.lit("LoF")
            ).otherwise(F.col("directionOnTrait"))
        ).withColumn("datasourceAll", F.lit("All")
        ).withColumn("niceName", F.col("datasourceId")
        ).replace(replacement_dict, subset=["niceName"]).persist()
print("assessment file ready")
doe_sources = [
    "gwas_credible_set",
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

def datasets_numbers_trait(assessment, replacement_dict, buckets_number):
    """This function creates in a long format (suitable for R) the N of evidences and association per
    DoE section. At the end, it creates a column with the corresponding deciles of the numbers to plot their intesntiy
    The deciles are trained using the assoc to avoid underrating numbers from the assoc respecting the evidences

    prueba_assessment = (
        directionOfEffect(evidences, platform_v)
        .withColumn("datasourceAll", F.lit("All"))
        .withColumn("niceName", F.col("datasourceId"))
        .replace(replacement_dict, subset=["niceName"])
    )
    """ 

    ## direction on trait
    unpivot_trait = "stack(2, 'protect', protect, 'risk', risk) as (type, count)"
    trait = (
        assessment.replace(replacement_dict, subset=["niceName"])
        .filter(F.col("directionOnTrait") != "noEvaluable")
        .groupBy("niceName")
        .pivot("directionOnTrait")
        .count()
        .union(
            assessment.filter(F.col("directionOnTrait") != "noEvaluable")
            .groupBy("datasourceAll")
            .pivot("directionOnTrait")
            .count()
            .withColumnRenamed("datasourceAll", "niceName")
        )
        .select("niceName", F.expr(unpivot_trait))
        .fillna(0)
    ).withColumn("facet", F.lit("trait"))

    #### direction on target
    unpivot_function = "stack(2, 'gof', gof, 'lof', lof) as (type, count)"

    function = (
        assessment.filter(F.col("variantEffect") != "noEvaluable")
        .groupBy("niceName")
        .pivot("variantEffect")
        .count()
        .union(
            assessment.filter(F.col("variantEffect") != "noEvaluable")
            .groupBy("datasourceAll")
            .pivot("variantEffect")
            .count()
            .withColumnRenamed("datasourceAll", "niceName")
        )
        .select("niceName", F.expr(unpivot_function))
        .fillna(0)
    ).withColumn("facet", F.lit("function"))

    #### direction complete
    unpivot_whole = "stack(4, 'LoF_protect', LoF_protect, 'LoF_risk', LoF_risk,'GoF_protect',GoF_protect,'GoF_risk',GoF_risk) as (type, count)"

    whole = (
        assessment.filter(F.col("homogenized") != "noEvaluable")
        .groupBy("targetId", "diseaseId", "niceName", "homogenized")
        .count()
        .groupBy("niceName")
        .pivot("homogenized")
        .count()
        .union(
            assessment.filter(F.col("homogenized") != "noEvaluable")
            .groupBy("targetId", "diseaseId", "datasourceAll", "homogenized")
            .count()
            .groupBy("datasourceAll")
            .pivot("homogenized")
            .count()
            .withColumnRenamed("datasourceAll", "niceName")
        )
        .select("niceName", F.expr(unpivot_whole))
        .fillna(0)
    ).withColumn("facet", F.lit("whole"))

    qds = QuantileDiscretizer(
        numBuckets=buckets_number,
        inputCol="count",
        outputCol="deciles",
    )

    all = function.union(trait).union(whole)
    result = qds.fit(whole).transform(all)  ### train qds in Whole and transform all

    return result
print("function for doe Table ready")
bucket_number=7
print("running function with bucket numbers =", bucket_number)
result=datasets_numbers_trait(assessment, replacement_dict, bucket_number)

print("printing file")
file_name = f"gs://ot-team/jroldan/analysis/doeTableNumbers.csv"
result.toPandas().to_csv(file_name)
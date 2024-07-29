#######
## ANALYSIS FOR L2G Scores, genetic evidence and Direction of Effect
## Original, propagated and Other Vs Oncology
#######
from functions import discrepancifier
from DoEAssessment import directionOfEffect
from functions import relative_success
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
import matplotlib.pyplot as plt
from datetime import date
from pyspark.sql.types import StructType, StructField, StringType, IntegerType
from datetime import datetime


spark = SparkSession.builder.getOrCreate()
c = datetime.now()
print("spark session created at", c)


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

terminated_array = (
    terminated.groupBy("targetId", "diseaseId")
    .agg(F.collect_set("clinicalStatus").alias("clinicalStatus"))
    .withColumn("prediction", F.when(F.col("clinicalStatus").isNotNull(), F.lit("yes")))
)

### Now , filter by rank, and join with the info from Ot genetics and run the DoE.
ranking = Window.partitionBy("studyId", "variantId")
### union with the other datasources
platform_v = "24.06"

target_path = (
    f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/targets/"
)
target = spark.read.parquet(target_path)

disease_path = (
    f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/diseases/"
)
diseases = spark.read.parquet(disease_path)
mecact_path = f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/mechanismOfAction/"
mecact = spark.read.parquet(mecact_path)
evidences = spark.read.parquet(
    f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/evidence"
).filter(
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
        .otherwise(F.lit("noEvaluable")),  ####
    )
    .withColumnRenamed("id", "target_id")
)

# 2# run the transformation of the evidences datasets used.

windowSpec = Window.partitionBy("targetId", "diseaseId")

columns_chembl = ["LoF_protect", "GoF_protect"]
columns_dataset = ["LoF_protect", "GoF_protect", "LoF_risk", "GoF_risk", "evidenceDif"]
columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect"]
terms = ["noEvaluable", "bivalent_risk", "null", "dispar"]

sincgc = [
    "gene_burden",
    "intogen",
    "eva",
    "eva_somatic",
    "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
]

germline = [
    "gene_burden",
    "eva",
    "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
]

somatic = ["intogen", "cancer_gene_census", "eva_somatic"]

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
    "chembl",
    "WOcgc",
    "somatic",
    "germline",
]
#### version all gene burden
prueba_assessment = (
    directionOfEffect(evidences, platform_v)
    .withColumn(
        "rank",
        F.when(
            (F.col("datasourceId") == "ot_genetics_portal"),
            F.row_number().over(ranking.orderBy(F.col("resourceScore").desc())),
        ).otherwise(F.lit(None)),
    )
    .withColumn(
        "average",
        F.when(
            (F.col("datasourceId") == "ot_genetics_portal"),
            F.avg("resourceScore").over(ranking.orderBy(F.col("resourceScore").desc())),
        ).otherwise(F.lit(None)),
    )
    .persist()
)

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

v2g = spark.read.parquet("gs://genetics-portal-dev-data/22.09.1/outputs/v2g")
varDistToGene = v2g.select(
    F.concat_ws("_", "chr_id", "position", "ref_allele", "alt_allele").alias(
        "variantId"
    ),
    F.col("gene_id").alias("targetId"),
    "source_id",
    "d",
    "distance_score",
).filter(F.col("source_id") == "canonical_tss")

ranking = Window.partitionBy("studyId", "variantId")


#######
# Build Ot genetics dataset as supporting evidence
#######
otGenetics = (
    prueba_assessment.filter(
        F.col("datasourceId").isin(
            [
                "ot_genetics_portal",
            ]
        )
    )
    # .filter((F.col("homogenized") != "noEvaluable"))
    .join(varDistToGene, on=["variantId", "targetId"], how="left")
    .withColumn(
        "datasources",
        F.collect_set("datasourceId").over(Window.partitionBy("targetId", "diseaseId")),
    )
    .withColumn(
        "L2G_ranking",
        F.when(
            (F.col("datasourceId") == "ot_genetics_portal"),
            F.row_number().over(ranking.orderBy(F.col("resourceScore").desc())),
        ).otherwise(F.lit(None)),
    )
    .withColumn(
        "averageL2G",
        F.when(
            (F.col("datasourceId") == "ot_genetics_portal"),
            F.avg("resourceScore").over(ranking.orderBy(F.col("resourceScore").desc())),
        ).otherwise(F.lit(None)),
    )
    .withColumn(
        "averageCanonicalTSSDistance",
        F.when(
            (F.col("datasourceId") == "ot_genetics_portal"),
            F.avg("d").over(ranking.orderBy(F.col("resourceScore").desc())),
        ).otherwise(F.lit(None)),
    )
    .withColumn(
        "datasources",
        F.when(
            F.col("rank").isNull(),
            F.array_remove(F.col("datasources"), "ot_genetics_portal"),
        ).otherwise(F.col("datasources")),
    )
    .withColumn(
        "distance_ranking",
        F.when(
            (F.col("datasourceId") == "ot_genetics_portal"),
            F.row_number().over(ranking.orderBy(F.col("d").asc())),
        ).otherwise(F.lit(None)),
    )
    .withColumn(
        "ChemblL2gRanking",
        F.when(
            (F.array_contains(F.col("datasources"), "chembl"))
            & (F.array_contains(F.col("datasources"), "ot_genetics_portal")),
            F.lit(F.col("L2G_ranking")),
        ).otherwise(F.lit(None)),
    )
    .withColumn(
        "chemblDistanceRanking",
        F.when(
            (F.array_contains(F.col("datasources"), "chembl"))
            & (F.array_contains(F.col("datasources"), "ot_genetics_portal")),
            F.lit(F.col("distance_ranking")),
        ).otherwise(F.lit(None)),
    )
    .withColumn(
        "frontierValue",
        ## ot genetics portal
        F.when(
            F.col("datasourceId") == "ot_genetics_portal",  ### the same for gene_burden
            F.when(
                (F.col("beta").isNotNull()) & (F.col("OddsRatio").isNull()),
                F.when(
                    (F.col("beta") <= 0.1) & (F.col("beta") >= -0.1),
                    F.lit("limitValue"),
                ).otherwise(F.lit("noLimitValue")),
            )
            .when(
                (F.col("beta").isNull()) & (F.col("OddsRatio").isNotNull()),
                F.when(
                    (F.col("OddsRatio") <= 1.1) & (F.col("OddsRatio") >= 0.9),
                    F.lit("limitValue"),
                ).otherwise(F.lit("noLimitValue")),
            )
            .when(
                (F.col("beta").isNull()) & (F.col("OddsRatio").isNull()),
                F.lit("noValue"),
            ),
        ),
    )
).persist()

#####
# function for interpreting DoE and coherencies/discrepancies
#####

diseases2 = diseases.select("id", "parents").withColumn(
    "diseaseIdPropagated",
    F.explode_outer(F.concat(F.array(F.col("id")), F.col("parents"))),
)

analysis_chembl = discrepancifier(
    prueba_assessment.filter((F.col("datasourceId") == "chembl"))
    .withColumn(
        "maxClinPhase",
        F.max(F.col("clinicalPhase")).over(Window.partitionBy("targetId", "diseaseId")),
    )
    .groupBy("targetId", "diseaseId", "maxClinPhase")
    .pivot("homogenized")
    .agg(F.count("targetId"))
    .persist()
)

#### propag OtGenetics:
otGenetics_propag = (
    otGenetics.filter((F.col("datasourceId") == "ot_genetics_portal"))
    .join(
        diseases2.selectExpr("id as diseaseId", "diseaseIdPropagated"),
        on="diseaseId",
        how="left",
    )
    .withColumnRenamed("diseaseId", "oldDiseaseId")
    .withColumnRenamed("diseaseIdPropagated", "diseaseId")
).persist()


#### include dictionary for calling dataframes:
# max_L2GScore
# min_distance_ranking


def benchmarkOT(discrepancifier, otGenetics, metric):
    dict_comb = {}
    dict_comb = {
        "diagonalAgreeWithDrugs": f"{metric}",
        "oneCellAgreeWithDrugs": f"{metric}",
    }
    list_l2g = [
        0.1,
        0.15,
        0.2,
        0.25,
        0.3,
        0.35,
        0.4,
        0.45,
        0.5,
        0.55,
        0.6,
        0.65,
        0.7,
        0.75,
        0.8,
        0.85,
        0.9,
        0.95,
    ]
    list_dist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    return (
        discrepancifier(
            otGenetics.filter((F.col("datasourceId") == "ot_genetics_portal"))
            .withColumn(
                "min_distance_ranking",
                F.min("distance_ranking").over(
                    Window.partitionBy("targetId", "diseaseId")
                ),
            )
            .withColumn(  ### take maximum L2G score per T-D
                "max_L2GScore",
                F.max("resourceScore").over(
                    Window.partitionBy("targetId", "diseaseId")
                ),
            )
            .groupBy(
                "targetId",
                "diseaseId",
                f"{value}",
            )  ##### modifications here to include the groups of ranking/distances to TSS
            .pivot("homogenized")
            .agg(F.count("targetId"))
        )
        .selectExpr(
            "targetId",
            "diseaseId",
            f"{metric}",
            "coherencyDiagonal as coherencyDiagonal",
            "coherencyOneCell as coherencyOneCell",
            "LoF_protect",
            "GoF_protect",
            "LoF_risk",
            "GoF_risk",
        )
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
                F.col(f"{metric}").isNotNull(), F.lit("hasGeneticEvidence")
            ).otherwise(F.lit("noGeneticEvidence")),
        )
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
                F.col("diagonalAgreeWithDrugs") == "coherent", F.lit("yes")
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "oneCellYes",
            F.when(
                F.col("oneCellAgreeWithDrugs") == "coherent", F.lit("yes")
            ).otherwise(F.lit("no")),
        )
        .select(
            ["*"]
            + (
                [  ### single columns
                    F.when(F.col(f"{metric}") >= n, F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{metric}>={str(n).replace('.', '_')}")
                    for n in list_l2g
                ]
                if metric == "max_L2GScore"  # Adjust this condition as needed
                else [
                    F.when(F.col(f"{metric}") <= n, F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{metric}<={n}")
                    for n in list_dist
                ]
            )
            + (
                [  ### column combinations
                    F.when((F.col(a) == "coherent") & (F.col(x) >= n), F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{x}>={str(n).replace('.', '_')}&{a}_combined")
                    for a, x in dict_comb.items()
                    for n in list_l2g
                ]
                if metric == "max_maxL2GScore"
                else [
                    F.when((F.col(a) == "coherent") & (F.col(x) <= n), F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{x}<={str(n).replace('.', '_')}&{a}_combined")
                    for a, x in dict_comb.items()
                    for n in list_dist
                ]
            )
        )
        .persist()
    )


metric_list = ["max_L2GScore", "min_distance_ranking"]
datasetDict = {}
for value in metric_list:
    if value == "max_L2GScore":
        datasetDict[f"df_l2g_original"] = benchmarkOT(
            discrepancifier, otGenetics, value
        )
        datasetDict[f"df_l2g_propagated"] = benchmarkOT(
            discrepancifier, otGenetics_propag, value
        )
    elif value == "min_distance_ranking":
        datasetDict[f"df_distance_original"] = benchmarkOT(
            discrepancifier, otGenetics, value
        )
        datasetDict[f"df_distance_propagated"] = benchmarkOT(
            discrepancifier, otGenetics_propag, value
        )


def comparisons_df(test_propag) -> list:
    """Return list of all comparisons to be used in the analysis"""
    toAnalysis = test_propag.drop(
        "Phase4",
        "Phase>=3",
        "Phase>=2",
        "Phase>=1",
        "Phase0",
        "clinicalStatus",
        "prediction",
        "count",
        "PhaseT",
        "taLabelSimple",
    ).columns[17:]
    dataType = ["byDatatype"] * len(toAnalysis)
    l_studies = []
    l_studies.extend([list(a) for a in zip(toAnalysis, dataType)])

    schema = StructType(
        [
            StructField("comparison", StringType(), True),
            StructField("comparisonType", StringType(), True),
        ]
    )

    comparisons = spark.createDataFrame(l_studies, schema=schema)
    ### include all the columns as predictor

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
    wComparison = Window.partitionBy(F.col(comparisonColumn))
    wPrediction = Window.partitionBy(F.col(predictionColumn))
    wPredictionComparison = Window.partitionBy(
        F.col(comparisonColumn), F.col(predictionColumn)
    )

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
listado = []
#### Run analysis using loop
for key, df in datasetDict.items():
    print(key)
    df = df.persist()
    aggSetups_original = comparisons_df(df)
    for row in aggSetups_original:
        aggregations_original(df, key, listado, *row, today_date)
    df.unpersist()
    print(key + " df unpersisted")

##### read files and make spreadsheet
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio


def convertTuple(tup):
    st = ",".join(map(str, tup))
    return st


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
#### update the dictionary dfs with other columns included in the analysis
key_list = [
    "hasGeneticEvidence",
    "OneCell",
    "diagonal",
    "max_L2GScore",
    "min_distance_ranking",
]
value_list = [
    "geneticEvidence",
    "oneCellDoE",
    "diagonalDoE",
    "L2GScore",
    "TSSDistance",
]

dfs = {}


def create_dict_column(dfs, key_list, value_list):
    if len(key_list) != len(value_list):
        raise ValueError("lists of different length")
    dfs.update(zip(key_list, value_list))
    return dfs


dfs = create_dict_column(dfs, key_list, value_list)

# Define the lists of possible substrings
phase_opt = [
    "Phase4",
    "Phase>=3",
    "Phase>=2",
    "Phase>=1",
    "PhaseT",
]

result = []
result_st = []
result_ci = []
array2 = []

# Initialize an empty list to store the results
results = []

# Iterate over the sample strings and extract the desired substrings
for path in listado:
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

    # Check for comparison options
    for key, value in dfs.items():
        if key in path and "combined" in path:
            if "oneCell" in path:
                comparison = dfs[key]
                group = "combinedOneCell"
            elif "diagonal" in path:
                comparison = dfs[key]
                group = "combinedDiagonal"
        elif key in path:
            comparison = dfs[key]
            group = key

    # Check for phase options
    for substr in phase_opt:
        if substr in path:
            phase = substr

    if "original" in path:
        dimension = "original"
    elif "propag" in path:
        dimension = "propagated"
    if "l2g" in path:
        dataset = "l2gScore"
    elif "distance" in path:
        dataset = "TSSdistance"

    results.append(
        [
            group,
            dataset,
            comparison,
            dimension,
            phase,
            round(float(resX.split(",")[0]), 2),
            round(float(resX.split(",")[1]), 8),
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
        "group",
        "dataset",
        "comparison",
        "dimension",
        "phase",
        "path",
    ],
)
df = pd.DataFrame(
    results,
    columns=[
        "group",
        "dataset",
        "comparison",
        "dimension",
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

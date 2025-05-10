
#### for matrix

import time
#from array import ArrayType
from functions import (
    discrepancifier,
    temporary_directionOfEffect,
)
# from stoppedTrials import terminated_td
from DoEAssessment import directionOfEffect
# from membraneTargets import target_membrane
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
#from itertools import islice
from datetime import datetime
from datetime import date
from pyspark.sql.types import (
    StructType,
    StructField,
    DoubleType,
    StringType,
    IntegerType,
    ArrayType
)
import pandas as pd

spark = SparkSession.builder.getOrCreate()

platform_v = "25.03"

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

### ingest evidences and create gwas_credible_set evidences. 

path_n='gs://open-targets-data-releases/25.03/output/'

target = spark.read.parquet(f"{path_n}target/")

diseases = spark.read.parquet(f"{path_n}disease/")

evidences = spark.read.parquet(f"{path_n}evidence")

credible = spark.read.parquet(f"{path_n}credible_set")

new = spark.read.parquet(f"{path_n}colocalisation_coloc") 

index=spark.read.parquet(f"{path_n}study/")

variantIndex = spark.read.parquet(f"{path_n}variant")

biosample = spark.read.parquet(f"{path_n}biosample")

print("loaded files")

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
)
    # .persist()
)

print("loaded newColoc")

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
)  # .persist()

print("loaded gwasComplete")

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
    #"gwas_credible_set", remove so avoid potential duplicates as it will be incorporated later (DoE is done separately)
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

window_spec = Window.partitionBy("targetId", "diseaseId",'leftStudyId').orderBy( ### include gwas study
    F.col("pValueExponent").asc()
)
gwasCredibleAssoc = (
    resolvedColoc.withColumn(
        "homogenized", F.first("colocDoE", ignorenulls=True).over(window_spec)
    )  ## added 30.01.2025
    .select("targetId", "diseaseId",'leftStudyId', "homogenized")
    .withColumn(
        "homogenized",
        F.when(F.col("homogenized").isNull(), F.lit("noEvaluable")).otherwise(
            F.col("homogenized")
        ),
    )
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
diseases2 = diseases.select("id", "parents").withColumn(
    "diseaseIdPropagated",
    F.explode_outer(F.concat(F.array(F.col("id")), F.col("parents"))),
)

chembl_trials = (
    assessment.filter((F.col("datasourceId").isin(["chembl"])))
    .groupBy("targetId", "diseaseId")
    .agg(F.max(F.col("clinicalPhase")).alias("maxClinPhase"))
)

negativeTD = (
    evidences.filter(F.col("datasourceId") == "chembl")
    .select("targetId", "diseaseId", "studyStopReason", "studyStopReasonCategories")
    .filter(F.array_contains(F.col("studyStopReasonCategories"), "Negative"))
    .groupBy("targetId", "diseaseId")
    .count()
    .withColumn("stopReason", F.lit("Negative"))
    .drop("count")
)

assessment_all = assessment.unionByName(
    gwasCredibleAssoc.withColumn("datasourceId", F.lit("gwas_credible_set")),
    allowMissingColumns=True,
)

######## 
 #evidences_all = spark.read.parquet(
 #    f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/evidence"
 #)
###
replacement_dict = {
    "gene_burden": "GeneBurden",
    "chembl": "ChEMBL",
    "intogen": "Intogen",
    "orphanet": "Orphanet",
    "cancer_gene_census": "CancerGeneCensus",
    "eva": "EvaGermline",
    "gene2phenotype": "Gene2Phenotype",
    "eva_somatic": "EvaSomatic",
    "gwas_credible_set": "OtGenetics",
    "impc": "IMPC",
}

### take only the ones with datasources for DoE
evidences_all = assessment_all.filter(F.col("datasourceId").isin(doe_sources)).replace(
    replacement_dict, subset=["datasourceId"]
)


def coincidence_matrix(evidences_all, replacement_dict):
    """Build a coincidence matrix of target disease associations between the datasources
    (these datasources are defined previously in the "evidences file passed to this function)

    Arguments:
        evidences {spark DataFrame} -- spark dataframe
        but_keep_these {list} -- list of columns to keep without checking for nulls

    Returns:
       Three spark DataFrame.show(), content: intradatasource disparities, interdatasource disparities
       and interdatasource coherency matrix coincidence matrix  -- dataframe.show()
    """

    terms = ["noEvaluable", "bivalent_risk", "null", "dispar"]
    columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect"]

    dataset1 = (  #### unique identifier and global coherency
        discrepancifier(
            # directionOfEffect(evidences, platform_v)
            evidences_all
            .withColumn("datasourceAll", F.lit("All"))
            .withColumn("niceName", F.col("datasourceId"))
            .replace(replacement_dict, subset=["niceName"])
            .filter((F.col("homogenized")).isin(terms) == False)
            .groupBy("targetId", "diseaseId")
            .pivot("homogenized")
            .agg(F.count("targetId"))
        )
        .withColumn("id", F.monotonically_increasing_id())
        .withColumnRenamed("coherencyDiagonal", "coherency_inter")
        .withColumnRenamed("coherencyOneCell", "coherency_onecell")
    ).persist()

    ### coherency intra datasource
    dataset2 = (
        discrepancifier(
            #directionOfEffect(evidences, platform_v)
            evidences_all
            .withColumn("datasourceAll", F.lit("All"))
            .withColumn("niceName", F.col("datasourceId"))
            .replace(replacement_dict, subset=["niceName"])
            .filter((F.col("homogenized")).isin(terms) == False)
            .groupBy("targetId", "diseaseId", "datasourceId")
            .pivot("homogenized")
            .agg(F.count("targetId"))
        )
        .withColumnRenamed("coherencyDiagonal", "coherency_intra")
        .withColumnRenamed("coherencyOneCell", "coherency_intra_OneCell")
    ).persist()

    ### two diferent dataset 3 depending on coherency type
    dataset3 = (
        dataset2.filter(F.col("coherency_intra") == "coherent")
        .join(dataset1, on=["targetId", "diseaseId"], how="left")
        .drop(*columns, "count")
    ).persist()

    #### cross all rows between them:
    dataset4 = (
        dataset2.join(dataset1, on=["targetId", "diseaseId"], how="left").drop(
            *columns, "count"
        )
    ).persist()

    datasetAllShared = (
        dataset2.drop("noEvaluable")
        .join(dataset1.drop("noEvaluable"), on=["targetId", "diseaseId"], how="left")
        .drop(*columns, "count")
    ).persist()

    #### matrix of Total T-D shared between datasources
    allSharedTD = (
        datasetAllShared.groupBy("datasourceId")
        .agg(F.collect_set(F.col("id")).alias("ids"))
        .selectExpr("datasourceId as datasourceId_x", "ids as ids_x")
        .join(
            datasetAllShared.groupBy("datasourceId")
            .agg(F.collect_set(F.col("id")).alias("ids"))
            .selectExpr("datasourceId as datasourceId_y", "ids as ids_y")
        )
        .withColumn(
            "sharedTD",
            F.size(F.array_intersect(F.col("ids_x"), F.col("ids_y"))),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("sharedTD")))
    ).persist()

    #### matrix of intradatasource disparities:
    matrix_intraDispar = (
        dataset4.groupBy("datasourceId")
        .pivot("coherency_intra")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_intra",
            "dispar as dispar_intra",
        )
        .join(
            dataset4.groupBy("datasourceId")
            .pivot("coherency_intra")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_intra_y",
                "dispar as dispar_intra_y",
            )
        )
        .withColumn(
            "nrs_disparIntra",
            F.size(F.array_intersect(F.col("dispar_intra"), F.col("dispar_intra_y"))),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_disparIntra")))
    )

    ### matrix of interdatasource disparities (including intra)

    matrix_interDispar = (
        dataset4.groupBy("datasourceId")
        .pivot("coherency_inter")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_inter",
            "dispar as dispar_inter",
        )
        .join(
            dataset4.groupBy("datasourceId")
            .pivot("coherency_inter")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_inter_y",
                "dispar as dispar_inter_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("dispar_inter"), F.col("dispar_inter_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherent_inter"), F.col("coherent_inter_y"))
            ),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_disparInter")))  #### change depending on interest
    )

    ### matrix of coincidences of coherencies

    dataset_matrix_interCoherent = (
        dataset4.groupBy("datasourceId")
        .pivot("coherency_inter")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_inter",
            "dispar as dispar_inter",
        )
        .join(
            dataset4.groupBy("datasourceId")
            .pivot("coherency_inter")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_inter_y",
                "dispar as dispar_inter_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("dispar_inter"), F.col("dispar_inter_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherent_inter"), F.col("coherent_inter_y"))
            ),
        )
        .withColumn(
            "nrs_Total", (F.col("nrs_disparInter") + F.col("nrs_coherentInter"))
        )
        .withColumn(
            "interCoherencyPercentage",
            (F.col("nrs_coherentInter") / F.col("nrs_Total") * 100),
        )
    )

    matrix_interCoherent = (
        dataset_matrix_interCoherent.groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_coherentInter")))  #### change depending on interest
    )
    matrix_interCoherentPercentage = (
        dataset_matrix_interCoherent.groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(
            F.first(F.col("interCoherencyPercentage"))
        )  #### change depending on interest
    )

    #### prepare dataset for one Cell coherencies
    dataset3OneCell = (
        dataset2.filter(F.col("coherency_intra_OneCell") == "coherent")
        .join(dataset1, on=["targetId", "diseaseId"], how="left")
        .drop(*columns, "count")
    ).persist()

    #### matrix of coherencies using oneCell criteria including only intracoherent data
    dataset_matrix_oneCellCoherent_onlyIntraCoherent = (
        dataset3OneCell.groupBy("datasourceId")
        .pivot("coherency_onecell")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherentOneCell",
            "dispar as disparOneCell",
        )
        .join(
            dataset3OneCell.groupBy("datasourceId")
            .pivot("coherency_onecell")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherentOneCell_y",
                "dispar as disparOneCell_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("disparOneCell"), F.col("disparOneCell_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherentOneCell"), F.col("coherentOneCell_y"))
            ),
        )
        .withColumn(
            "nrs_Total", (F.col("nrs_disparInter") + F.col("nrs_coherentInter"))
        )
        .withColumn(
            "percentageCoherencyInter",
            (F.col("nrs_coherentInter") / F.col("nrs_Total") * 100),
        )
    ).persist()

    matrix_oneCellCoherent_onlyIntraCoherent = (
        dataset_matrix_oneCellCoherent_onlyIntraCoherent.groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_coherentInter")))  #### change depending on interest
    ).persist()

    matrix_oneCellCoherent_onlyIntraCoherentPercentage = (
        dataset_matrix_oneCellCoherent_onlyIntraCoherent.groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(
            F.first(F.col("percentageCoherencyInter"))
        )  #### change depending on interest
    ).persist()

    ### the next are matrix combining coherency criteria and intradatasource disparities:
    matrix_interDispar_onlyIntraCoherentOneCell = (
        dataset3OneCell.groupBy("datasourceId")
        .pivot("coherency_inter")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_inter",
            "dispar as dispar_inter",
        )
        .join(
            dataset3OneCell.groupBy("datasourceId")
            .pivot("coherency_inter")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_inter_y",
                "dispar as dispar_inter_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("dispar_inter"), F.col("dispar_inter_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherent_inter"), F.col("coherent_inter_y"))
            ),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_disparInter")))  #### change depending on interest
    )

    matrix_interDispar_onlyIntraCoherent = (
        dataset3.groupBy("datasourceId")
        .pivot("coherency_inter")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_inter",
            "dispar as dispar_inter",
        )
        .join(
            dataset3.groupBy("datasourceId")
            .pivot("coherency_inter")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_inter_y",
                "dispar as dispar_inter_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("dispar_inter"), F.col("dispar_inter_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherent_inter"), F.col("coherent_inter_y"))
            ),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_disparInter")))  #### change depending on interest
    )

    matrix_interCoherent_onlyIntraCoherent = (
        dataset3.groupBy("datasourceId")
        .pivot("coherency_inter")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_inter",
            "dispar as dispar_inter",
        )
        .join(
            dataset3.groupBy("datasourceId")
            .pivot("coherency_inter")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_inter_y",
                "dispar as dispar_inter_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("dispar_inter"), F.col("dispar_inter_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherent_inter"), F.col("coherent_inter_y"))
            ),
        )
        .withColumn(
            "nrs_Total", (F.col("nrs_disparInter") + F.col("nrs_coherentInter"))
        )
        .withColumn(
            "percentageCoherencyInter",
            (F.col("nrs_coherentInter") / F.col("nrs_Total") * 100),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(
            F.first(F.col("percentageCoherencyInter"))
        )  #### change depending on interest
    )

    matrix_total_onlyIntraCoherent = (
        dataset3.groupBy("datasourceId")
        .pivot("coherency_inter")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_inter",
            "dispar as dispar_inter",
        )
        .join(
            dataset3.groupBy("datasourceId")
            .pivot("coherency_inter")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_inter_y",
                "dispar as dispar_inter_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("dispar_inter"), F.col("dispar_inter_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherent_inter"), F.col("coherent_inter_y"))
            ),
        )
        .withColumn(
            "nrs_Total", (F.col("nrs_disparInter") + F.col("nrs_coherentInter"))
        )
        .withColumn(
            "percentageCoherencyInter",
            (F.col("nrs_coherentInter") / F.col("nrs_Total") * 100),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_Total")))  #### change depending on interest
    )

    return (
        allSharedTD,  ### N T-D shared between all DS
        matrix_intraDispar,  ### intradatasource disparities
        matrix_interDispar,  ### using diagonal criteria, from coherent TD per DS, dispar TD shared between DS
        matrix_interCoherent,  ### using diagonal criteria, from coherent TD per DS, coherent TD shared between DS
        matrix_interCoherentPercentage,  ### same as above but % of coherent across all TD shared when being coherent intra DS
        matrix_oneCellCoherent_onlyIntraCoherent,  ### using one cell criteria, from coherent TD per DS, coherent TD shared between DS
        matrix_oneCellCoherent_onlyIntraCoherentPercentage,  ### same as above but % of coherent TD over dispar between DS, when being coherent intra DS
        matrix_interDispar_onlyIntraCoherentOneCell,  ### using one cell criteria to check intra coherency, N of TD dispar shared between DS
        matrix_interDispar_onlyIntraCoherent,  #### using diagonal criteria to check intra coherency, N of TD dispar shared between DS
        matrix_interCoherent_onlyIntraCoherent,  #### using diagonal criteria, percentage of shared T-D between datasources that are coherent
        matrix_total_onlyIntraCoherent,  ### using diagonal criteria, sum of shared TD dispar and coherent between DS
    )

print('loading matrices')
(
    allSharedTD,
    matrix_intraDispar,
    matrix_interDispar,
    matrix_interCoherent,
    matrix_interCoherentPercentage,
    matrix_oneCellCoherent_onlyIntraCoherent,
    matrix_oneCellCoherent_onlyIntraCoherentPercentage,
    matrix_interDispar_onlyIntraCoherentOneCell,
    matrix_interDispar_onlyIntraCoherent,
    matrix_interCoherent_onlyIntraCoherent,
    matrix_total_onlyIntraCoherent
) = coincidence_matrix(evidences_all, replacement_dict)

print('start making csv and saving')

print('allSharedTD')

allSharedTD.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/allsharedTD.csv"
)
print('matrix_intraDispar')
matrix_intraDispar.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_intraDispar.csv"
)
print('matrix_interDispar')
matrix_interDispar.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_interDispar.csv"
)
print('matrix_interCoherent')
matrix_interCoherent.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_interCoherent.csv"
)
print('matrix_interCoherentPercentage')
matrix_interCoherentPercentage.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_interCoherentPercentage.csv"
)
print('matrix_oneCellCoherent_onlyIntraCoherent')
matrix_oneCellCoherent_onlyIntraCoherent.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_oneCellCoherent_onlyIntraCoherent.csv"
)
print('matrix_oneCellCoherent_onlyIntraCoherentPercentage')
matrix_oneCellCoherent_onlyIntraCoherentPercentage.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_oneCellCoherent_onlyIntraCoherentPercentage.csv"
)
print('matrix_interDispar_onlyIntraCoherentOneCell')
matrix_interDispar_onlyIntraCoherentOneCell.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_interDispar_onlyIntraCoherentOneCell.csv"
)
print('matrix_interDispar_onlyIntraCoherent')
matrix_interDispar_onlyIntraCoherent.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_interDispar_onlyIntraCoherent.csv"
)
print('matrix_interCoherent_onlyIntraCoherent')
matrix_interCoherent_onlyIntraCoherent.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_interCoherent_onlyIntraCoherent.csv"
)
print('matrix_total_onlyIntraCoherent')
matrix_total_onlyIntraCoherent.toPandas().to_csv(
    "gs://ot-team/jroldan/analysis/discrepancies/matrix_total_onlyIntraCoherent.csv"
)
print('finished!')   
#### discrepancies matrices:
from DoEAssessment import directionOfEffect
from functions import discrepancifier
from pyspark.sql import SparkSession
import pyspark.sql.functions as F

spark = SparkSession.builder.getOrCreate()

platform_v = "24.06"

doe_sources = [
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

evidences = spark.read.parquet(
    f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/evidence"
)
### take only the ones with datasources for DoE
evidences = evidences.filter(F.col("datasourceId").isin(doe_sources))

replacement_dict = {
    "gene_burden": "GeneBurden",
    "chembl": "ChEMBL",
    "intogen": "Intogen",
    "orphanet": "Orphanet",
    "cancer_gene_census": "CancerGeneCensus",
    "eva": "EvaGermline",
    "gene2phenotype": "Gene2Phenotype",
    "eva_somatic": "EvaSomatic",
    "ot_genetics_portal": "OtGenetics",
    "impc": "IMPC",
}


def coincidence_matrix(
    evidences,
):
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
            directionOfEffect(evidences, platform_v)
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
            directionOfEffect(evidences, platform_v)
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

    dataset3OneCell = (
        dataset2.filter(F.col("coherency_intra_OneCell") == "coherent")
        .join(dataset1, on=["targetId", "diseaseId"], how="left")
        .drop(*columns, "count")
    ).persist()

    #### cross all rows between them:
    dataset4 = (
        dataset2.join(dataset1, on=["targetId", "diseaseId"], how="left").drop(
            *columns, "count"
        )
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

    matrix_interCoherent = (
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
        .agg(F.first(F.col("nrs_coherentInter")))  #### change depending on interest
    )

    #### matrix of coherencies using oneCell criteria including only intracoherent data
    matrix_oneCellCoherent_onlyIntraCoherent = (
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
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(
            F.first(F.col("percentageCoherencyInter"))
        )  #### change depending on interest
    )

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
        matrix_intraDispar,
        matrix_interDispar,
        matrix_interCoherent,
        matrix_oneCellCoherent_onlyIntraCoherent,
        matrix_interDispar_onlyIntraCoherentOneCell,
        matrix_interDispar_onlyIntraCoherent,
        matrix_interCoherent_onlyIntraCoherent,
        matrix_total_onlyIntraCoherent,
    )

""" Build datasets to plot % of coherency"""

from DoEAssessment import directionOfEffect
from functions import discrepancifier
from pyspark.sql import SparkSession
import pyspark.sql.functions as F
from pyspark.ml.feature import QuantileDiscretizer

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


def dataset_percentages(evidences, platform_v, replacement_dict):
    #### include nice name
    prueba_assessment = (
        directionOfEffect(evidences, platform_v)
        .withColumn("datasourceAll", F.lit("All"))
        .withColumn("niceName", F.col("datasourceId"))
        .replace(replacement_dict, subset=["niceName"])
    )
    diagonal = (
        discrepancifier(
            prueba_assessment.filter(F.col("homogenized") != "noEvaluable")
            .groupBy("targetId", "diseaseId", "datasourceId")
            .pivot("homogenized")
            .count()
        )
        .drop("noEvaluable")
        .groupBy("datasourceId")
        .pivot("coherencyDiagonal")
        .count()
        .sort(F.col("datasourceId").desc())
        .fillna(0)
        .withColumn(
            "percentageCoherent",
            F.round(F.col("coherent") / (F.col("coherent") + F.col("dispar")) * 100, 2),
        )
        .withColumn("type", F.lit("diagonal"))
    )

    oneCell = (
        discrepancifier(
            prueba_assessment.filter(F.col("homogenized") != "noEvaluable")
            .groupBy("targetId", "diseaseId", "datasourceId")
            .pivot("homogenized")
            .count()
        )
        .drop("noEvaluable")
        .groupBy("datasourceId")
        .pivot("coherencyOneCell")
        .count()
        .sort(F.col("datasourceId").desc())
        .fillna(0)
        .withColumn(
            "percentageCoherent",
            F.round(F.col("coherent") / (F.col("coherent") + F.col("dispar")) * 100, 2),
        )
        .withColumn("type", F.lit("oneCell"))
    )
    return (
        diagonal.union(oneCell)
        .withColumn("niceName", F.col("datasourceId"))
        .replace(replacement_dict, subset=["niceName"])
    )

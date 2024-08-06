""" Build datasets with DoE evidence and assoc numbers"""

from DoEAssessment import directionOfEffect
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


def datasets_numbers(evidences, platform_v, replacement_dict):
    #### include nice name
    prueba_assessment = (
        directionOfEffect(evidences, platform_v)
        .withColumn("datasourceAll", F.lit("All"))
        .withColumn("niceName", F.col("datasourceId"))
        .replace(replacement_dict, subset=["niceName"])
    )

    ## direction on trait
    unpivot_trait = "stack(2, 'protect', protect, 'risk', risk) as (type, count)"
    trait = (
        prueba_assessment.filter(F.col("directionOnTrait") != "noEvaluable")
        .groupBy("niceName")
        .pivot("directionOnTrait")
        .count()
        .union(
            prueba_assessment.filter(F.col("directionOnTrait") != "noEvaluable")
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
        prueba_assessment.filter(F.col("variantEffect") != "noEvaluable")
        .groupBy("niceName")
        .pivot("variantEffect")
        .count()
        .union(
            prueba_assessment.filter(F.col("variantEffect") != "noEvaluable")
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
        prueba_assessment.filter(F.col("homogenized") != "noEvaluable")
        .groupBy("niceName")
        .pivot("homogenized")
        .count()
        .union(
            prueba_assessment.filter(F.col("homogenized") != "noEvaluable")
            .groupBy("datasourceAll")
            .pivot("homogenized")
            .count()
            .withColumnRenamed("datasourceAll", "niceName")
        )
        .select("niceName", F.expr(unpivot_whole))
        .fillna(0)
    ).withColumn("facet", F.lit("whole"))

    qds = QuantileDiscretizer(
        numBuckets=10,
        inputCol="count",
        outputCol="deciles",
    )

    df = trait.union(function).union(whole)
    result = qds.fit(df).transform(df).presist()

    return result

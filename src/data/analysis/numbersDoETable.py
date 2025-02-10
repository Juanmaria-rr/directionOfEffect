""" Build datasets with DoE evidence and assoc numbers"""

#from DoEAssessment import directionOfEffect
from functions import temporary_directionOfEffect,build_gwasResolvedColoc
from pyspark.sql import SparkSession
import pyspark.sql.functions as F, Window 
from pyspark.ml.feature import QuantileDiscretizer

spark = SparkSession.builder.getOrCreate()

#platform_v = "24.12"
path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"

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
    "gwas_credible_set" :"credibleOtGenetics"   
}

doe_sources = [
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

evidences = spark.read.parquet(
    f"{path}/evidence"
)
### take only the ones with datasources for DoE
evidences = evidences.filter(F.col("datasourceId").isin(doe_sources))


prueba_assessment = (
        temporary_directionOfEffect(evidences, doe_sources)
        .withColumn("datasourceAll", F.lit("All"))
    )

gwasResolvedColoc = build_gwasResolvedColoc(path)

#### take the direction from the lowest p value
window_spec = Window.partitionBy("targetId", "diseaseId").orderBy(
    F.col("pValueExponent").asc()
)
gwasCredibleAssoc = (
    gwasResolvedColoc.withColumn(
        "homogenized", F.first("colocDoE", ignorenulls=True).over(window_spec)
    )  ## added 30.01.2025
    .select("targetId", "diseaseId", "homogenized")
    .withColumn(
        "homogenized",
        F.when(F.col("homogenized").isNull(), F.lit("noEvaluable")).otherwise(
            F.col("homogenized")
        ),
    )
)

### replace     
assessment = prueba_assessment.unionByName(
    gwasCredibleAssoc.withColumn("datasourceId", F.lit("gwas_credible_set")),
    allowMissingColumns=True,
).withColumn("niceName", F.col("datasourceId")).replace(replacement_dict, subset=["niceName"])

print("asessment done")

def datasets_numbers_trait(prueba_assessment, buckets_number):
    """This function creates in a long format (suitable for R) the N of evidences and association per
    DoE section. At the end, it creates a column with the corresponding deciles of the numbers to plot their intesntiy
    The deciles are trained using the assoc to avoid underrating numbers from the assoc respecting the evidences
    """
    """ 
        usage of temporary_directionOfEffect:
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
        .groupBy("targetId", "diseaseId", "niceName", "homogenized")
        .count()
        .groupBy("niceName")
        .pivot("homogenized")
        .count()
        .union(
            prueba_assessment.filter(F.col("homogenized") != "noEvaluable")
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

print("read function")

print("executing function")

dataset=datasets_numbers_trait(prueba_assessment, 7)
print("executed function")
print("printing dataset")
dataset.toPandas().to_csv("gs://ot-team/jroldan/analysis/numbersDoeTable.csv")
print("dataset saved succesfully")
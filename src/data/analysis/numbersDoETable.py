#from DoEAssessment import directionOfEffect
from functions import temporary_directionOfEffect,build_gwasResolvedColoc_noPropag
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
from pyspark.ml.feature import QuantileDiscretizer
from datetime import date
from doeFunction import doeFunction

spark = SparkSession.builder.getOrCreate()
today_date = str(date.today())

#platform_v = "25.03"
path = 'gs://open-targets-data-releases/25.03/output/'

replacement_dict = {
    "gene_burden": "GeneBurden",
    "chembl": "ChEMBL",
    "intogen": "Intogen",
    "orphanet": "Orphanet",
    "cancer_gene_census": "CancerGeneCensus",
    "eva": "EvaGermline",
    "gene2phenotype": "Gene2Phenotype",
    "eva_somatic": "EvaSomatic",
    "gwas_credible_sets": "OtGenetics",
    "impc": "IMPC",
}

doe_sources = [
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

new = spark.read.parquet(f"{path}colocalisation_coloc") 
credible = spark.read.parquet(f"{path}credible_set")
index=spark.read.parquet(f"{path}study/")
evidences = spark.read.parquet(f"{path}evidence")
diseases = spark.read.parquet(f"{path}disease/")

### take only the ones with datasources for DoE
evidences = evidences.filter(F.col("datasourceId").isin(doe_sources))

gwasCredibleAssoc,assessment=doeFunction(new,credible,index,evidences,diseases,path)

gwasCredibleAssoc.withColumn(
        "homogenized",
        F.when(F.col("homogenized").isNull(), F.lit("noEvaluable")).otherwise(
            F.col("homogenized")
        ),
    ).withColumn("variantEffect", 
                F.when(F.col("homogenized").isin(["LoF_protect","LoF_risk"]), F.lit("LoF")
                ).when(F.col("homogenized").isin(["GoF_protect","GoF_risk"]), F.lit("GoF")
                ).otherwise(F.lit("noEvaluable"))
    ).withColumn("directionOnTrait", 
                F.when(F.col("homogenized").isin(["LoF_protect","GoF_protect"]), F.lit("protect")
                ).when(F.col("homogenized").isin(["LoF_risk","GoF_risk"]), F.lit("risk")
                ).otherwise(F.lit("noEvaluable"))
    ).select("targetId", "diseaseId", "homogenized","variantEffect","directionOnTrait")


### replace     
assessment_all = assessment.unionByName(
    gwasCredibleAssoc.withColumn("datasourceId", F.lit("gwas_credible_set")),
    allowMissingColumns=True,
    ).withColumn("niceName", F.col("datasourceId")
    ).replace(replacement_dict, subset=["niceName"]
    ).withColumn("datasourceAll", F.lit("All"))

print("asessment done")

def datasets_numbers_trait(assessment_all, buckets_number):
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
        assessment_all.filter(F.col("directionOnTrait") != "noEvaluable")
        .groupBy("niceName")
        .pivot("directionOnTrait")
        .count()
        .union(
            assessment_all.filter(F.col("directionOnTrait") != "noEvaluable")
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
        assessment_all.filter(F.col("variantEffect") != "noEvaluable")
        .groupBy("niceName")
        .pivot("variantEffect")
        .count()
        .union(
            assessment_all.filter(F.col("variantEffect") != "noEvaluable")
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
        assessment_all.filter(F.col("homogenized") != "noEvaluable")
        .groupBy("targetId", "diseaseId", "niceName", "homogenized")
        .count()
        .groupBy("niceName")
        .pivot("homogenized")
        .count()
        .union(
            assessment_all.filter(F.col("homogenized") != "noEvaluable")
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

dataset=datasets_numbers_trait(assessment_all, 11)
print("executed function")
print("printing dataset")
dataset.toPandas().to_csv(f"gs://ot-team/jroldan/analysis/{today_date}_numbersDoeTable.csv")
print("dataset saved succesfully")
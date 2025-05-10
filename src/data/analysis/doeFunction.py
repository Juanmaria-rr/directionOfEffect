'''
spark = SparkSession.builder.getOrCreate()
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
path_n='gs://open-targets-data-releases/25.03/output/'

new = spark.read.parquet(f"{path_n}colocalisation_coloc") 
credible = spark.read.parquet(f"{path_n}credible_set")
index=spark.read.parquet(f"{path_n}study/")
evidences = spark.read.parquet(f"{path_n}evidence")
diseases = spark.read.parquet(f"{path_n}disease/")
'''
def doeFunction(path_n,new,credible,index,evidences,diseases):
    #### function to get genetic associations from gwas_Credible_set 
    ### and DoE across currently used datasources
    from functions import (
    temporary_directionOfEffect,
    )

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
    assessment_all = assessment.unionByName(
    gwasCredibleAssoc.withColumn("datasourceId", F.lit("gwas_credible_set")),
    allowMissingColumns=True,
    )

    return assessment_all

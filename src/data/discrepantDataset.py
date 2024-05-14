""" From evidences """

from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
from pyspark.sql.types import StructType, StructField, StringType

spark = SparkSession.builder.getOrCreate()

target_path = "gs://open-targets-data-releases/24.03/output/etl/parquet/targets/"
target = spark.read.parquet(target_path)
disease_path = "gs://open-targets-data-releases/24.03/output/etl/parquet/diseases/"
diseases = spark.read.parquet(disease_path)

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


def discrepancifier(df):
    """
    detect discrepancies per row where there are the four
    DoE assessments using Null and isNotNull assessments
    """
    columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect"]

    for col in columns:
        if col not in df.columns:
            df = df.withColumn(col, F.lit(None)).persist()

    return df.withColumn(
        "coherencyDiagonal",
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
            ((F.col("GoF_protect").isNotNull()) & (F.col("LoF_protect").isNotNull())),
            F.lit("dispar"),
        )
        .otherwise(F.lit("coherent")),
    ).withColumn(
        "coherencyOneCell",
        F.when(
            F.lit("LoF_risk").isNotNull()
            & (
                (F.col("LoF_protect").isNull())
                & (F.col("GoF_risk").isNull())
                & (F.col("GoF_protect").isNull())
            ),
            F.lit("coherent"),
        )
        .when(
            F.lit("GoF_risk").isNotNull()
            & (
                (F.col("LoF_protect").isNull())
                & (F.col("LoF_risk").isNull())
                & (F.col("GoF_protect").isNull())
            ),
            F.lit("coherent"),
        )
        .when(
            F.lit("LoF_protect").isNotNull()
            & (
                (F.col("LoF_risk").isNull())
                & (F.col("GoF_risk").isNull())
                & (F.col("GoF_protect").isNull())
            ),
            F.lit("coherent"),
        )
        .when(
            F.lit("GoF_protect").isNotNull()
            & (
                (F.col("LoF_protect").isNull())
                & (F.col("GoF_risk").isNull())
                & (F.col("LoF_risk").isNull())
            ),
            F.lit("coherent"),
        )
        .otherwise(F.lit("dispar")),
    )


def discrepant_dataset(discrepanficier, evidences):
    """return a dataframe of T-D associations annotating discrepancies
    depending on coherency Diagonal or OneCell
    """

    datatest = discrepanficier(
        evidences.filter(F.col("homogenized") != "noEvaluable")
        .withColumn(
            "datasources",
            F.collect_set("datasourceId").over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .groupBy("targetId", "diseaseId", "diseaseFromSource", "datasources")
        .pivot("homogenized")
        .agg(F.count("targetId"))
    )

    test = datatest.join(
        diseases.selectExpr("id as diseaseId", "name", "therapeuticAreas"),
        on="diseaseId",
        how="left",
    )

    return (
        test.select(
            "*",
            F.explode_outer(F.col("therapeuticAreas")).alias("therapeuticAreas_expl"),
        )
        .join(
            taDf.selectExpr(
                "taId as therapeuticAreas_expl", "taLabel", "taLabelSimple"
            ),
            on="therapeuticAreas_expl",
            how="left",
        )
        .groupBy(
            "targetId",
            "diseaseId",
            # "diseaseFromSource",
            "name",
            "GoF_protect",
            "GoF_risk",
            "LoF_protect",
            "LoF_risk",
            "coherencyDiagonal",
            "coherencyOneCell",
            "datasources",
        )
        .agg(
            F.collect_set("therapeuticAreas_expl").alias("therapeuticAreas"),
            F.collect_set("taLabel").alias("taName"),
            F.collect_set("taLabelSimple").alias("taLabelSimple"),
        )
        .filter(
            (F.col("coherencyDiagonal") == "dispar")
            | (F.col("coherencyOneCell") == "dispar")
        )
        .join(
            target.selectExpr("id as targetId", "approvedSymbol"),
            on="targetId",
            how="left",
        )
        .persist()
    )

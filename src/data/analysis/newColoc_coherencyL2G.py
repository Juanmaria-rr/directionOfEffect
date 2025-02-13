#### BUILDING THE NEW GWAS GENETIC EVIDENCE FROM COLOC
from functions import (
    discrepancifier,
    build_gwasResolvedColoc,
    temporary_directionOfEffect,
    spreadSheetFormatter,
)
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
import matplotlib.pyplot as plt
from datetime import date, datetime
from pyspark.sql.types import (
    StructType,
    StructField,
    ArrayType,
    DoubleType,
    DecimalType,
    StringType,
    FloatType,
    IntegerType,
)
import pandas as pd

spark = SparkSession.builder.getOrCreate()

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

path = "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/"

evidences = (
    spark.read.parquet(f"{path}evidence").filter(
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
    # .persist()
)

gwasResolvedColoc = build_gwasResolvedColoc(path)

#### take the direction from the lowest p value
window_spec = Window.partitionBy("targetId", "diseaseId").orderBy(
    F.col("pValueExponent").asc()
)

print("creating new gwasCredibleAssoc")

### modify to include more information
gwasCredibleAssoc = (
    gwasResolvedColoc.withColumn(
        "homogenized", F.first("colocDoE", ignorenulls=True).over(window_spec)
    )
    .select(
        "targetId",
        "diseaseId",
        "homogenized",
        "leftStudyLocusId",
        "h4",
        "datasourceId",
        "resourceScore",
        "leftVariantId",
        "credibleLeftStudyType",
    )
    .withColumn(
        "homogenized",
        F.when(F.col("homogenized").isNull(), F.lit("noEvaluable")).otherwise(
            F.col("homogenized")
        ),
    )
)  ### there will be duplicates TargetId-DiseaseId because we are taking the most significant DoE

#### LOAD STUDYLOCUSID AND VARIANT DISTANCES
l2gPred = spark.read.parquet(
    "gs://open-targets-pre-data-releases/24.12-uo_test-3/output/etl/parquet/locusToGenePredictions"
)
l2gTable = (
    l2gPred.select("studyLocusId", "geneId", F.explode_outer("locusToGeneFeatures"))
    .filter(F.col("key").isin(["distanceFootprintMean", "distanceTssMean"]))
    .groupBy("studyLocusId", "geneId")
    .pivot("key")
    .agg(F.first("value"))
)
print("creating gwasCredibleAssocDistances")
gwasCredibleAssocDistances = gwasCredibleAssoc.join(
    l2gTable.withColumnRenamed("studyLocusId", "leftStudyLocusId").withColumnRenamed(
        "geneId", "targetId"
    ),
    on=["leftStudyLocusId", "targetId"],
    how="left",
)

print("creating analysis_chembl")

datasource_filter = ["chembl"]
preanalysis_chembl, evidences, actionType, oncolabel = temporary_directionOfEffect(
    path, datasource_filter
)
analysis_chembl = (
    discrepancifier(
        preanalysis_chembl.withColumn(
            "maxClinPhase",
            F.max(F.col("clinicalPhase")).over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .groupBy("targetId", "diseaseId", "maxClinPhase")
        .pivot("homogenized")
        .agg(F.count("targetId"))
        .persist()
    )
    .filter(  ### ensure drug has annotated MoA and is coherent per Target-Disease
        ((F.col("GoF_protect").isNotNull()) | F.col("LoF_protect").isNotNull())
        & (F.col("coherencyDiagonal") == "coherent")
    )
    .selectExpr(
        "targetId",
        "diseaseId",
        "maxClinPhase",
        "coherencyDiagonal as coherencyDiagonal_ch",
        "coherencyOneCell as coherencyOneCell_ch",
        "LoF_protect as LoF_protect_ch",
        "GoF_protect as GoF_protect_ch",
    )
)
diseases = spark.read.parquet(f"{path}diseases/")
diseases2 = diseases.select("id", "parents").withColumn(
    "diseaseIdPropagated",
    F.explode_outer(F.concat(F.array(F.col("id")), F.col("parents"))),
)

### pivot colocdoE grouping by T-D-studyLocusId-distances

value = "max_L2GScore"

list_l2g = [
    0.10,
    0.15,
    0.20,
    0.25,
    0.30,
    0.35,
    0.40,
    0.45,
    0.50,
    0.55,
    0.60,
    0.65,
    0.70,
    0.72,
    0.74,
    0.75,
    0.76,
    0.78,
    0.80,
    0.82,
    0.84,
    0.85,
    0.86,
    0.88,
    0.90,
    0.92,
    0.94,
    0.95,
    0.96,
    0.98,
]
print("creating benchmarkOT function")

dict_comb = {}


def benchmarkOT(
    dict_comb,
    value,
    gwasCredibleAssocDistances,
    analysis_chembl,
    terminated_array,
    list_l2g,
):
    dict_comb = {
        "hasGeneticEvidence": f"{value}",
        "diagonalYes": f"{value}",
        "oneCellYes": f"{value}",
        "L2GAndColoc": f"{value}",
    }
    min_value = gwasCredibleAssocDistances.agg(
        F.min("distanceFootprintMean")
    ).collect()[0][
        0
    ]  ## take min value for scalating
    dataset_original = discrepancifier(
        gwasCredibleAssocDistances
        # .filter(F.col("h4").isNotNull()) #### not filter by this because we want to include the L2G AND Coloc question
        .withColumn(  ### take maximum L2G score per T-D
            "max_L2GScore",
            F.max("resourceScore").over(Window.partitionBy("targetId", "diseaseId")),
        )
        .withColumn("min_distanceFootprintMean", F.lit(min_value))
        .withColumn(
            "scaled_distanceFootprintMean",
            (F.col("distanceFootprintMean") - F.col("min_distanceFootprintMean"))
            / (1 - F.col("min_distanceFootprintMean")),
        )
        .withColumn(
            "minDistFootPrintMean",
            F.min("scaled_distanceFootprintMean").over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .withColumn(
            "minTssDistance",
            F.min("distanceTssMean").over(Window.partitionBy("targetId", "diseaseId")),
        )
        .groupBy(
            "targetId",
            "diseaseId",
            f"{value}",
            # "leftStudyLocusId",
        )
        .pivot("homogenized")
        .count()
    )
    dataset_propag = discrepancifier(
        gwasCredibleAssocDistances
        # .filter(F.col("h4").isNotNull()) #### not filter by this because we want to include the L2G AND Coloc question
        .withColumn(  ### take maximum L2G score per T-D
            "max_L2GScore",
            F.max("resourceScore").over(Window.partitionBy("targetId", "diseaseId")),
        )
        .withColumn("min_distanceFootprintMean", F.lit(min_value))
        .withColumn(
            "scaled_distanceFootprintMean",
            (F.col("distanceFootprintMean") - F.col("min_distanceFootprintMean"))
            / (1 - F.col("min_distanceFootprintMean")),
        )
        .withColumn(
            "minDistFootPrintMean",
            F.min("scaled_distanceFootprintMean").over(
                Window.partitionBy("targetId", "diseaseId")
            ),
        )
        .withColumn(
            "minTssDistance",
            F.min("distanceTssMean").over(Window.partitionBy("targetId", "diseaseId")),
        )
        .join(
            diseases2.selectExpr("id as diseaseId", "diseaseIdPropagated"),
            on="diseaseId",
            how="left",
        )
        .withColumnRenamed("diseaseId", "oldDiseaseId")
        .withColumnRenamed("diseaseIdPropagated", "diseaseId")
        .groupBy(
            "targetId",
            "diseaseId",
            f"{value}",
            # "leftStudyLocusId",
        )
        .pivot("homogenized")
        .count()
    )

    return (
        dataset_original.join(
            analysis_chembl, on=["targetId", "diseaseId"], how="right"
        )
        .withColumn(
            "diagonalAgreeWithDrugs",
            F.when(
                (
                    (F.col("coherencyDiagonal_ch") == "coherent")
                    & (F.col("coherencyDiagonal") == "coherent")
                )
                # & (F.col("coherencyDiagonal") == "coherent")
                ,
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
                (
                    (F.col("coherencyOneCell_ch") == "coherent")
                    & (F.col("coherencyDiagonal") == "coherent")
                )
                # & (F.col("coherencyOneCell") == "coherent")
                ,
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
            # ).filter(
            #    F.col("diagonalAgreeWithDrugs").isNotNull()
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
        .withColumn(
            "hasGeneticEvidence",
            F.when(F.col(f"{value}").isNotNull(), F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "L2GAndColoc",
            F.when(
                (F.col(f"{value}").isNotNull())
                & (F.col("coherencyDiagonal").isin(["coherent", "dispar"])),
                F.lit("yes"),
            ).otherwise(F.lit("no")),
        )
        .join(terminated_array, on=["targetId", "diseaseId"], how="left")
        .withColumn(
            "PhaseT",
            F.when(F.col("prediction") == "yes", F.lit("yes")).otherwise(F.lit("no")),
        )
        .select(
            ["*"]
            + (
                [  ### single columns
                    F.when(F.col(f"{value}") >= n, F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{value}>={str(n).replace('.', '_')}")
                    for n in list_l2g
                ]
            )
            + (
                [  ### column combinations for Yes/No colums Plus has DoE (any agreement)
                    F.when((F.col(a) == "yes") & (F.col(x) >= n), F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{x}>={str(n).replace('.', '_')}&{a}_combined")
                    for a, x in dict_comb.items()
                    for n in list_l2g
                ]
            )
        ),
        dataset_propag.join(analysis_chembl, on=["targetId", "diseaseId"], how="right")
        .withColumn(
            "diagonalAgreeWithDrugs",
            F.when(
                (
                    (F.col("coherencyDiagonal_ch") == "coherent")
                    & (F.col("coherencyDiagonal") == "coherent")
                )
                # & (F.col("coherencyDiagonal") == "coherent")
                ,
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
                (
                    (F.col("coherencyOneCell_ch") == "coherent")
                    & (F.col("coherencyDiagonal") == "coherent")
                )
                # & (F.col("coherencyOneCell") == "coherent")
                ,
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
            # ).filter(
            #    F.col("diagonalAgreeWithDrugs").isNotNull()
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
        .withColumn(
            "hasGeneticEvidence",
            F.when(F.col(f"{value}").isNotNull(), F.lit("yes")).otherwise(F.lit("no")),
        )
        .withColumn(
            "L2GAndColoc",
            F.when(
                (F.col(f"{value}").isNotNull())
                & (F.col("coherencyDiagonal").isin(["coherent", "dispar"])),
                F.lit("yes"),
            ).otherwise(F.lit("no")),
        )
        .join(terminated_array, on=["targetId", "diseaseId"], how="left")
        .withColumn(
            "PhaseT",
            F.when(F.col("prediction") == "yes", F.lit("yes")).otherwise(F.lit("no")),
        )
        .select(
            ["*"]
            + (
                [  ### single columns
                    F.when(F.col(f"{value}") >= n, F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{value}>={str(n).replace('.', '_')}")
                    for n in list_l2g
                ]
            )
            + (
                [  ### column combinations for Yes/No colums Plus has DoE (any agreement)
                    F.when((F.col(a) == "yes") & (F.col(x) >= n), F.lit("yes"))
                    .otherwise(F.lit("no"))
                    .alias(f"{x}>={str(n).replace('.', '_')}&{a}_combined")
                    for a, x in dict_comb.items()
                    for n in list_l2g
                ]
            )
        ),
    )


print("creating dataframes in loop")

values = [
    "max_L2GScore",
    "minDistFootPrintMean",
    "minTssDistance",
]

datasetDict = {}
datasetDict_propag = {}

for value in values:
    if value == "max_L2GScore":
        (
            datasetDict[f"max_L2GScore_original"],
            datasetDict_propag[f"max_L2GScore_propag"],
        ) = benchmarkOT(
            dict_comb,
            value,
            gwasCredibleAssocDistances,
            analysis_chembl,
            terminated_array,
            list_l2g,
        )
    else:
        datasetDict[f"{value}_original"], datasetDict_propag[f"{value}_propag"] = (
            benchmarkOT(
                dict_comb,
                value,
                gwasCredibleAssocDistances,
                analysis_chembl,
                terminated_array,
                list_l2g,
            )
        )


def comparisons_df(dataset) -> list:
    """Return list of all comparisons to be used in the analysis"""
    toAnalysis = dataset.drop("clinicalStatus", "prediction").columns[22:]
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


from functions import relative_success, spreadSheetFormatter, convertTuple

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
    filePath = "gs://ot-team/jroldan/" + str(
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
    array1 = np.delete(
        out.join(full_data, on=["prediction", "comparison"], how="outer")
        .groupBy("comparison")
        .pivot("prediction")
        .agg(F.first("a"))
        .sort(F.col("comparison").desc())
        .select("comparison", "yes", "no")
        .fillna(0)
        .toPandas()
        .to_numpy(),
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

    results.append(
        [
            data,
            comparisonColumn,
            predictionColumn,
            round(float(resX.split(",")[0]), 2),
            float(resX.split(",")[1]),
            round(float(resx_CI.split(",")[0]), 2),
            round(float(resx_CI.split(",")[1]), 2),
            str(total),
            np.array(res_npPhaseX).tolist(),
            round(float(rs_result), 2),
            round(float(rs_ci[0]), 2),
            round(float(rs_ci[1]), 2),
            filePath,
        ]
    )
    return results


import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

# Initialize an empty list to store the results
result_st = []
result_ci = []
results = []


def convertTuple(tup):
    st = ",".join(map(str, tup))
    return st


print("launched function to run analysis")


listado = []
today_date = str(date.today())
aggSetups_original = comparisons_df(datasetDict["max_L2GScore_original"])

for key, df_analysis in datasetDict.items():
    print("corresponding dataframe key: ", key)
    df_analysis.persist()
    for row in aggSetups_original:
        print(key)
        aggregations_original(df_analysis, key, listado, *row, today_date)
    df_analysis.unpersist()
    print("unpersist df")
for key, df_analysis in datasetDict_propag.items():
    print("corresponding dataframe key: ", key)
    df_analysis.persist()
    for row in aggSetups_original:
        print(key)
        aggregations_original(df_analysis, key, listado, *row, today_date)
    df_analysis.unpersist()
    print("unpersist df")

print("finished analysis")
print("creating pandas dataframe with resulting rows")
df = pd.DataFrame(
    results,
    columns=[
        "type",
        "criteria",
        "phase",
        "OR",
        "pValue",
        "LowCI",
        "HighCI",
        "total",
        "array",
        "rs",
        "lowRs",
        "HighRs",
        "path",
    ],
)

schema = StructType(
    [
        StructField("type", StringType(), True),
        StructField("criteria", StringType(), True),
        StructField("phase", StringType(), True),
        StructField("oddsRatio", DoubleType(), True),
        StructField("pValue", DoubleType(), True),
        StructField("lowerInterval", DoubleType(), True),
        StructField("upperInterval", DoubleType(), True),
        StructField("total", StringType(), True),
        StructField("values", ArrayType(ArrayType(IntegerType())), True),
        StructField("relSuccess", DoubleType(), True),
        StructField("rsLower", DoubleType(), True),
        StructField("rsUpper", DoubleType(), True),
        StructField("path", StringType(), True),
    ]
)

# Convert list of lists to DataFrame
print("created pandas dataframe")
print("converting to spark dataframe")
pattern = r"(\b[0-9]+)_(\d+\b)"  # Matches patterns like 0_1, 0_95, etc.
df = (
    spreadSheetFormatter(spark.createDataFrame(results, schema=schema))
    .withColumn("range", F.regexp_extract(F.col("criteria"), pattern, 0))
    .withColumn(
        "range", F.regexp_replace(F.col("range"), "_", ".")  ## substitute "_" by "."
    )
    .withColumn(
        "comparison",
        F.regexp_extract(F.col("criteria"), r"&(.+)", 1),  # take string after "&"
    )
    .withColumn(
        "comparison",
        F.regexp_replace(F.col("comparison"), "_combined$", ""),  ### trims "combined"
    )
    .withColumn(
        "dimension", F.regexp_extract(F.col("your_column"), r"_(original|propag)", 1)
    )
)


print("writting dataframe")
df.toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_newColoc_L2Gbenchmark.csv"
)

print("dataframe written \n Analysis finished")

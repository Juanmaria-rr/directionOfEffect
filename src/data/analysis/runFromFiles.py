from itertools import islice
import time
from functions import relative_success, spreadSheetFormatter, convertTuple
#from array import ArrayType
from functions import (
    relative_success,
    spreadSheetFormatter,
    discrepancifier,
    temporary_directionOfEffect,
)
from functions import relative_success, spreadSheetFormatter, convertTuple
import re
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio, relative_risk
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
import gcsfs


from pyspark.sql import SparkSession
import numpy as np
spark = SparkSession.builder.getOrCreate()

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
results=[]
result = []
result_st = []
result_ci = []
array2 = []
results = []

today_date = str(date.today())

fs = gcsfs.GCSFileSystem()  # Initialize Google Cloud Storage filesystem

folder_path = "gs://ot-team/jroldan/2025-06-25_analysis/*/"
parquet_files = fs.glob(f"{folder_path}*.parquet")

for path in parquet_files:
    path = 'gs://' + path
    print(f"Reading {path}")
    df = spark.read.parquet(path)
    array1 = np.delete(
        df.join(full_data, on=["prediction", "comparison"], how="outer")
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
    import re

    # 1. Extract the element between 'analysis/' and the next '/'
    match1 = re.search(r'analysis/([^/]+)/', path)
    part1 = match1.group(1) if match1 else None

    # 2. Extract the element after the second '/' until '_Phase'
    match2 = re.search(r'[^/]+/([^/_]+)_Phase', path)
    part2 = match2.group(1) if match2 else None

    # 3. Extract the element from 'Phase' until '.parquet'
    match3 = re.search(r'(Phase[^.]+)\.parquet', path)
    part3 = match3.group(1) if match3 else None

    print("Part 1:", part1)
    print("Part 2:", part2)
    print("Part 3:", part3)

    results.append(
        [
            part1, #data,
            part2, #comparisonColumn,
            part3, #predictionColumn,
            round(float(resX.split(",")[0]), 2),
            float(resX.split(",")[1]),
            round(float(resx_CI.split(",")[0]), 2),
            round(float(resx_CI.split(",")[1]), 2),
            str(total),
            np.array(res_npPhaseX).tolist(),
            round(float(rs_result), 2),
            round(float(rs_ci[0]), 2),
            round(float(rs_ci[1]), 2),
            path,
        ]
    )
    print('finished')

df_results = pd.DataFrame(
    results,
    columns=[
        "group",
        "comparison",
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
print("created pandas dataframe")
print("converting to spark dataframe")
print("preparing dataframe")

schema = StructType(
    [
        StructField("group", StringType(), True),
        StructField("comparison", StringType(), True),
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

print("read pattern variables")
df = spreadSheetFormatter(spark.createDataFrame(df_results, schema=schema))
print("processed spreadsheet")
print("writting the dataframe")

# Convert list of lists to DataFrame
# Regular expressions
'''
value_pattern = r"df_([^_]+)_"  # Extracts {value}
middle_pattern = r"df_[^_]+_([^_]+)_"  # Extracts middle part (All, Other, etc.)
suffix_pattern = r"(original|propag)$"  # Extracts suffix (original or propag)
'''

df.withColumn(
    "datasource",
    F.regexp_extract(F.col("group"), r"df_(.*?)_(All|Other|OtherNull|Oncology)_(propag|original)", 1)
).withColumn(
    "therArea",
    F.regexp_extract(F.col("group"), r"_(All|Other|OtherNull|Oncology)_", 1)
).withColumn(
    "type",
    F.regexp_extract(F.col("group"), r"_(propag|original)$", 1)
).toPandas().to_csv(
    f"gs://ot-team/jroldan/analysis/{today_date}_genEvidAnalysis_new_NoFileteredColocCaviar.csv"
)

print("dataframe written \n Analysis finished")
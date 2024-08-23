import pyspark.sql.functions as F


def discrepancifier(df):
    """
    detect discrepancies per row where there are the four
    DoE assessments using Null and isNotNull assessments
    """
    columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect", "noEvaluable"]

    for col in columns:
        if col not in df.columns:
            df = df.withColumn(col, F.lit(None)).persist()

    return df.withColumn(
        "coherencyDiagonal",
        F.when(
            (F.col("LoF_risk").isNull())
            & (F.col("LoF_protect").isNull())
            & (F.col("GoF_risk").isNull())
            & (F.col("GoF_protect").isNull())
            & (F.col("noEvaluable").isNull()),
            F.lit("noEvid"),
        )
        .when(
            (F.col("LoF_risk").isNull())
            & (F.col("LoF_protect").isNull())
            & (F.col("GoF_risk").isNull())
            & (F.col("GoF_protect").isNull())
            & (F.col("noEvaluable").isNotNull()),
            F.lit("EvidNotDoE"),
        )
        .when(
            (F.col("LoF_risk").isNotNull())
            | (F.col("LoF_protect").isNotNull())
            | (F.col("GoF_risk").isNotNull())
            | (F.col("GoF_protect").isNotNull()),
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
                (
                    (F.col("GoF_protect").isNotNull())
                    & (F.col("LoF_protect").isNotNull())
                ),
                F.lit("dispar"),
            )
            .otherwise(F.lit("coherent")),
        ),
    ).withColumn(
        "coherencyOneCell",
        F.when(
            (F.col("LoF_risk").isNull())
            & (F.col("LoF_protect").isNull())
            & (F.col("GoF_risk").isNull())
            & (F.col("GoF_protect").isNull())
            & (F.col("noEvaluable").isNull()),
            F.lit("noEvid"),
        )
        .when(
            (F.col("LoF_risk").isNull())
            & (F.col("LoF_protect").isNull())
            & (F.col("GoF_risk").isNull())
            & (F.col("GoF_protect").isNull())
            & (F.col("noEvaluable").isNotNull()),
            F.lit("EvidNotDoE"),
        )
        .when(
            (F.col("LoF_risk").isNotNull())
            | (F.col("LoF_protect").isNotNull())
            | (F.col("GoF_risk").isNotNull())
            | (F.col("GoF_protect").isNotNull()),
            F.when(
                F.col("LoF_risk").isNotNull()
                & (
                    (F.col("LoF_protect").isNull())
                    & (F.col("GoF_risk").isNull())
                    & (F.col("GoF_protect").isNull())
                ),
                F.lit("coherent"),
            )
            .when(
                F.col("GoF_risk").isNotNull()
                & (
                    (F.col("LoF_protect").isNull())
                    & (F.col("LoF_risk").isNull())
                    & (F.col("GoF_protect").isNull())
                ),
                F.lit("coherent"),
            )
            .when(
                F.col("LoF_protect").isNotNull()
                & (
                    (F.col("LoF_risk").isNull())
                    & (F.col("GoF_risk").isNull())
                    & (F.col("GoF_protect").isNull())
                ),
                F.lit("coherent"),
            )
            .when(
                F.col("GoF_protect").isNotNull()
                & (
                    (F.col("LoF_protect").isNull())
                    & (F.col("GoF_risk").isNull())
                    & (F.col("LoF_risk").isNull())
                ),
                F.lit("coherent"),
            )
            .otherwise(F.lit("dispar")),
        ),
    )


def relative_success(array1):
    from scipy.stats.contingency import relative_risk

    ### take numbers from array
    a, b = array1[0]
    c, d = array1[1]
    ####
    """
    Where zeros cause problems with computation of the relative risk or its standard error,
    0.5 is added to all cells (a, b, c, d) (Pagano & Gauvreau, 2000; Deeks & Higgins, 2010).
    """
    total_expo = a + b
    total_noExpo = c + d
    ### for cases when total_expo/total_noExpo = 0,
    ### we sum 1 to avoid errors an get at least 0 in the %
    if any(t == 0 for t in [total_expo, total_noExpo]):
        total_expo = total_expo + 1
        total_noExpo = total_noExpo + 1
        ### calculate relative success
        relative_success = relative_risk(a, total_expo, c, total_noExpo)
        ### calculate confidence intervals
        rs_ci = relative_risk(a, total_expo, c, total_noExpo).confidence_interval(
            confidence_level=0.95
        )
    else:

        ### calculate relative success
        relative_success = relative_risk(a, total_expo, c, total_noExpo)
        ### calculate confidence intervals
        rs_ci = relative_risk(a, total_expo, c, total_noExpo).confidence_interval(
            confidence_level=0.95
        )

    return relative_success.relative_risk, rs_ci


def spreadSheetFormatter(df):
    from pyspark.sql.functions import format_number
    from pyspark.sql.types import (
        DoubleType,
    )

    new_df = (
        df.withColumn(
            "significant",
            F.when(F.col("pValue") < 0.0001, F.lit("****"))
            .when((F.col("pValue") >= 0.0001) & (F.col("pValue") < 0.001), F.lit("***"))
            .when((F.col("pValue") >= 0.001) & (F.col("pValue") < 0.01), F.lit("**"))
            .when((F.col("pValue") >= 0.01) & (F.col("pValue") < 0.05), F.lit("*"))
            .when(F.col("pValue") >= 0.05, F.lit("ns")),
        )
        .withColumn(
            "writeFigure",
            F.concat(
                F.round(F.col("oddsRatio"), 2),
                F.lit(" "),
                F.lit("("),
                F.round(F.col("lowerInterval"), 2),
                F.lit("-"),
                F.round(F.col("upperInterval"), 2),
                F.lit(")"),
            ),
        )
        .withColumn("value_1", F.col("values").getItem(0).getItem(0))
        .withColumn("value_2", F.col("values").getItem(0).getItem(1))
        .withColumn("value_3", F.col("values").getItem(1).getItem(0))
        .withColumn("value_4", F.col("values").getItem(1).getItem(1))
        .withColumn("numerator", (F.col("value_1") + F.col("value_2")).cast("int"))
        .withColumn("denominator", (F.col("value_3") + F.col("value_4")).cast("int"))
        .withColumn("pValue", F.col("pValue").cast(DoubleType()))
        .withColumn(
            "pValue_formatted",
            F.when(
                F.col("pValue") < 0.0001, F.col("pValue").cast("string")
            )  # Check if value matches scientific notation range
            # .when(F.col("pValue") < 0.05, F.format_number(F.col("pValue"), 2))
            .otherwise(F.format_number(F.col("pValue"), 4)),
        )
        .withColumn(
            "numDen", F.concat_ws("/", F.col("numerator"), F.col("denominator"))
        )
    )
    return new_df

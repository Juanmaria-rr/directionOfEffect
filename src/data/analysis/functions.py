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

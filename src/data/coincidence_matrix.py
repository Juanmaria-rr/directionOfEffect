def coincidence_matrix(
    evidences_JR,
):
    """Build a coincidence matrix of target disease associations between the datasources
    (these datasources are defined previously in the "evidences file passed to this function)

    Arguments:
        evidences {spark DataFrame} -- spark dataframe
        but_keep_these {list} -- list of columns to keep without checking for nulls

    Returns:
       Three spark DataFrame.show(), content: intradatasource disparities, interdatasource disparities
       and interdatasource coherency matrix coincidence matrix  -- dataframe.show()
    """

    terms = ["noEvaluable", "bivalent_risk", "null", "dispar"]
    columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect"]

    dataset1 = (  #### unique identifier and global coherency
        evidences_JR.filter((F.col("homogenizedVersion_J")).isin(terms) == False)
        .groupBy("targetId", "diseaseId")
        .pivot("homogenizedVersion_J")
        .agg(F.count("targetId"))
        .select(
            F.col("targetId"),
            # F.col("datasourceId"),
            F.col("diseaseId"),
            *(F.col(c).cast("int").alias(c) for c in columns)
        )
        .withColumn(
            "coherency_inter",
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
        )
        .withColumn("id", monotonically_increasing_id())
    ).persist()

    ### coherency intra datasource

    dataset2 = (  ## coherency per datasource
        evidences_JR.filter((F.col("homogenizedVersion_J")).isin(terms) == False)
        .groupBy("targetId", "diseaseId", "datasourceId")
        .pivot("homogenizedVersion_J")
        .agg(F.count("targetId"))
        .select(
            F.col("targetId"),
            F.col("datasourceId"),
            F.col("diseaseId"),
            *(F.col(c).cast("int").alias(c) for c in columns)
        )
        .withColumn(
            "coherency_intra",
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
        )
    ).persist()

    #### cross all rows between them:

    dataset4 = (
        dataset2.join(dataset1, on=["targetId", "diseaseId"], how="left")
        .drop(*columns, "count")
        .persist()
    ).persist()

    #### matrix of intradatasource disparities:

    matrix_intraDispar = (
        dataset4.groupBy("datasourceId")
        .pivot("coherency_intra")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_intra",
            "dispar as dispar_intra",
        )
        .join(
            dataset4.groupBy("datasourceId")
            .pivot("coherency_intra")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_intra_y",
                "dispar as dispar_intra_y",
            )
        )
        .withColumn(
            "nrs_disparIntra",
            F.size(F.array_intersect(F.col("dispar_intra"), F.col("dispar_intra_y"))),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_disparIntra")))
    )

    ### matrix of interdatasource disparities (including intra)

    matrix_interDispar = (
        dataset4.groupBy("datasourceId")
        .pivot("coherency_inter")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_inter",
            "dispar as dispar_inter",
        )
        .join(
            dataset4.groupBy("datasourceId")
            .pivot("coherency_inter")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_inter_y",
                "dispar as dispar_inter_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("dispar_inter"), F.col("dispar_inter_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherent_inter"), F.col("coherent_inter_y"))
            ),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_disparInter")))  #### change depending on interest
    )

    ### matrix of coincidences of coherencies

    matrix_interCoherent = (
        dataset4.groupBy("datasourceId")
        .pivot("coherency_inter")
        .agg(F.collect_set(F.col("id")))
        .selectExpr(
            "datasourceId as datasourceId_x",
            "coherent as coherent_inter",
            "dispar as dispar_inter",
        )
        .join(
            dataset4.groupBy("datasourceId")
            .pivot("coherency_inter")
            .agg(F.collect_set(F.col("id")))
            .selectExpr(
                "datasourceId as datasourceId_y",
                "coherent as coherent_inter_y",
                "dispar as dispar_inter_y",
            )
        )
        .withColumn(
            "nrs_disparInter",
            F.size(F.array_intersect(F.col("dispar_inter"), F.col("dispar_inter_y"))),
        )
        .withColumn(
            "nrs_coherentInter",
            F.size(
                F.array_intersect(F.col("coherent_inter"), F.col("coherent_inter_y"))
            ),
        )
        .groupBy("datasourceId_x")
        .pivot("datasourceId_y")
        .agg(F.first(F.col("nrs_coherentInter")))  #### change depending on interest
    )

    return (
        print("matrix intraDispar"),
        matrix_intraDispar.show(),
        print("matrix interDispar"),
        matrix_interDispar.show(),
        print("matrix coherents"),
        matrix_interCoherent.show(),
    )
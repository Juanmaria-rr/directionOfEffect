def factorsAsColumns (df: DataFrame) -> DataFrame:
    """ 
    Convert each column element into a column with a Yes/No 
    content regarding the presenc of the element in the row, without
    the need of pivoting the dataframe
    'list_columns' should contain the column names with the elements 
    to be converted into columns. For instance:

        list_columns = 
        ["diseaseId",
        "approvedSymbol"]

    """
    import pyspark.sql.functions as F

    dfs = {}
    keys = []
    columns = []
    list_columns = [
        "diseaseId",
        "maxClinPhase",
        "sampleSize",
        "approvedSymbol",
    ]
    for x in list_columns:
        a = [data[0] for data in df.select(x).distinct().collect()]
        for y in a:
            keys.append(y)
            columns.append(x)

    keys2 = [str(s).replace(".0", "") for s in keys]

    for key, value in zip(keys2, columns):
        dfs[key] = value

    return df.select(
        ["*"]
        + [
            F.when((f"{x}") == F.col(c), F.lit("yes"))
            .when(F.col(c).isNull(), F.lit(None))
            .otherwise(F.lit("no"))
            .alias(f"col_{x}")
            for x, c in dfs.items()
        ]
    ).persist()
def drop_fully_null_columns(df, but_keep_these=[]):
    """Drops DataFrame columns that are fully null
    (i.e. the maximum value is null)

    Arguments:
        df {spark DataFrame} -- spark dataframe
        but_keep_these {list} -- list of columns to keep without checking for nulls

    Returns:
        spark DataFrame -- dataframe with fully null columns removed
    """

    # skip checking some columns
    cols_to_check = [col for col in df.columns if col not in but_keep_these]
    if len(cols_to_check) > 0:
        # drop columns for which the max is None
        rows_with_data = (
            df.select(*cols_to_check)
            .groupby()
            .agg(*[F.max(c).alias(c) for c in cols_to_check])
            .take(1)[0]
        )
        cols_to_drop = [
            c for c, const in rows_with_data.asDict().items() if const == None
        ]
        new_df = df.drop(*cols_to_drop)

        return new_df
    else:
        return df

import polars as pl

NUCLEOTIDES = {"A", "G", "C", "T", "a", "g", "c", "t"}
def sanitize_df(motif_df: pl.DataFrame) -> pl.DataFrame:
    motif_df = (
      motif_df
        .with_columns(
            sequence=pl.col("sequence").str.to_lowercase(),
            sequence_of_arm=pl.col("sequence_of_arm").str.to_lowercase()
        )
        .filter(
          pl.col("sequence")
            .map_elements(lambda seq: all(n in NUCLEOTIDES for n in seq),
                          return_dtype=pl.Boolean)
            )
       )
    return motif_df
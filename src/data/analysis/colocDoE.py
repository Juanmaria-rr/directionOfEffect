""" this scripts run the analysis for comparing QTL studies, tissues together with therapy areas matched"""
from functions import relative_success, spreadSheetFormatter
from stoppedTrials import terminated_td
from DoEAssessment import directionOfEffect
from membraneTargets import target_membrane
from pyspark.sql import SparkSession, Window
import pyspark.sql.functions as F
from datetime import datetime
from datetime import date
from pyspark.sql.types import StructType, StructField, StringType, IntegerType
from pyspark.sql.types import (
    StructType,
    StructField,
    DoubleType,
    DecimalType,
    StringType,
    FloatType,
)

spark = SparkSession.builder.getOrCreate()

platform_v = "24.06"

target_path = (
    f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/targets/"
)
target = spark.read.parquet(target_path)

disease_path = (
    f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/diseases/"
)
diseases = spark.read.parquet(disease_path)

evidences = spark.read.parquet(
    f"gs://open-targets-data-releases/{platform_v}/output/etl/parquet/evidence"
)

coloc = spark.read.parquet(
    "gs://genetics-portal-dev-data/22.09.1/outputs/v2d_coloc"
).filter(F.col("right_type") != "gwas")


### laod sample sizes - look at the document
sampleSize = spark.read.csv("gs://ot-team/jroldan/colocSampleSize.csv", header=True)

terminated = terminated_td(
    spark, "gs://ot-team/jroldan/analysis/targetDiseaseStoppedNegative.csv"
)
### load QTL tissues mapped to therapy areas
onto_samples = spark.read.csv(
    "gs://ot-team/jroldan/20240112_mappedSQLtissuesGTP3_5.csv", header=True
)

#### GSEA annotation for hallmark inflamation targets
immflam_annot = (
    spark.read.json(
        "gs://ot-team/jroldan/analysis/HALLMARK_INFLAMMATORY_RESPONSE.v2023.2.Hs_edited.json"
    )
    .select(F.explode_outer("geneSymbols").alias("approvedSymbol"))
    .withColumn("isInflam", F.lit("yes"))
)

#### Build the uniqueBetas_analyse
targetType = (
    target.select("id", "approvedSymbol", F.explode_outer("targetClass"))
    .select("id", "approvedSymbol", "col.label")
    .groupBy("id", "approvedSymbol")
    .agg(F.collect_set("label").alias("label"))
    .filter(F.col("label").isNotNull())
    .selectExpr(
        "id as targetIdtargetType",
        "approvedSymbol as approvedSymbol",
        "label as targetType",
    )
    .join(immflam_annot, on="approvedSymbol", how="left")
)

### take ontology of samples
samplesOnto = (
    onto_samples.withColumn(
        "right_bio_feature", F.split(F.col("original"), " - ").getItem(0)
    )
    .withColumn(
        "therapyArea", F.split(F.col("20231207_curated_simplified"), " - ").getItem(2)
    )
    .withColumn("EFO", F.split(F.col("20231207_curated_simplified"), " - ").getItem(1))
    .withColumn(
        "right_bio_feature2",
        F.split(F.col("20231207_curated_simplified"), " - ").getItem(0),
    )
    .drop(
        "curated_simplified",
        "20231207_curated_simplified",
        "original",
        "curated",
        "_c3",
        "_c4",
    )
)

coloc2 = (
    coloc.select(
        F.concat_ws("_", "left_chrom", "left_pos", "left_ref", "left_alt").alias(
            "left_locus_id"
        ),
        F.concat_ws("_", "right_chrom", "right_pos", "right_ref", "right_alt").alias(
            "right_locus_id"
        ),
        F.col("left_study").alias("left_study_id"),
        F.col("right_study").alias("right_study_id"),
        "right_gene_id",
        "coloc_h4",
        "left_var_right_study_beta",
        "right_phenotype",
        F.col("left_type"),
        F.col("right_type"),
        F.col("right_bio_feature"),
        F.col("is_flipped"),
        "left_var_right_study_pval",
    )
    .withColumn(
        "beta_assessed",  ### diferent from sQTL and oQTL
        F.when(
            (F.col("left_var_right_study_beta") > 0)
            & (F.col("right_study_id") != "GTEx-sQTL"),
            F.lit("gof"),
        ).when(
            (F.col("left_var_right_study_beta") < 0)
            & (F.col("right_study_id") != "GTEx-sQTL"),
            F.lit("lof"),
        )
        #### for sQTL is the opposite
        .when(
            (F.col("left_var_right_study_beta") > 0)
            & (F.col("right_study_id") == "GTEx-sQTL"),
            F.lit("lof"),
        )
        .when(
            (F.col("left_var_right_study_beta") < 0)
            & (F.col("right_study_id") == "GTEx-sQTL"),
            F.lit("gof"),
        )
        .otherwise(F.lit("neutral")),
    )
    .join(samplesOnto, on=["right_bio_feature"], how="left")
    .drop("right_bio_feature")
)
### check for disparities using count of different assessment for beta for target

disparities = coloc2.groupBy("left_locus_id", "left_study_id", "right_gene_id").agg(
    F.size(F.collect_set("beta_assessed")).alias("count"),
)
### add the label of which left_locus_id,left_study_id and right_gene_id are having contradictions

coloc3 = (
    coloc2.withColumnRenamed("right_bio_feature2", "right_bio_feature")
    .join(
        disparities, on=["left_locus_id", "left_study_id", "right_gene_id"], how="left"
    )
    .persist()
)
#### Run directionOfEffect
prueba_assessment = directionOfEffect(evidences, platform_v)

## add therapeuticArea name to  diseases
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
        ("OTAR_0000020", "nutritional or metabolic disease", "Other"),
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
)
#### prepare for disease propagation, but explode at the end
diseases2 = (
    diseases.withColumn("taId", F.explode_outer("therapeuticAreas"))
    .join(taDf.drop("taLabelSimple"), on="taId", how="left")
    .groupBy("id", "parents", "name")
    .agg(
        F.collect_set("taId").alias("taId"),
        F.collect_set("taLabel").alias("diseaseTherapyAreas"),
    )
)
#### load data from ot genetics with the assessments
delimiter = ","
evidence_ot_genetics = (
    (
        prueba_assessment.filter((F.col("datasourceId").isin(["ot_genetics_portal"])))
        .groupBy("targetId", "diseaseId", "variantId", "studyId")
        .pivot("directionOnTrait")
        .count()
        .drop("noEvaluable", "conflict/noEvaluable")
        .persist()
    )
    .join(diseases2, F.col("diseaseId") == diseases2.id, "left")
    .withColumnRenamed("diseaseId", "oldDiseaseId")
    .withColumn(
        "diseaseId",
        F.explode_outer(F.concat(F.array(F.col("id")), F.col("parents"))),
    )
    .drop("id")
)

coloc_otgene = (
    coloc3.withColumnRenamed("left_study_id", "studyId")
    .withColumnRenamed("left_locus_id", "locusId")
    .withColumnRenamed("right_gene_id", "targetId")
    .join(
        evidence_ot_genetics.withColumnRenamed("variantId", "locusId"),
        on=["studyId", "locusId", "targetId"],
        how="left",
    )
    .withColumnRenamed("GoF", "GoF_OT")  # remove
    .withColumnRenamed("LoF", "LoF_OT")  # remove
).persist()

chembl_trials = (
    prueba_assessment.filter((F.col("datasourceId").isin(["chembl"])))
    .groupBy("targetId", "diseaseId")
    .agg(F.max(F.col("clinicalPhase")).alias("maxClinPhase"))
)

chembl = (
    (
        prueba_assessment.filter(
            (F.col("datasourceId").isin(["chembl"]))
            & (F.col("homogenized").isin(["noEvaluable", "dispar"]) == False)
        )
    )
    .groupBy("targetId", "diseaseId")
    .pivot("variantEffect")
    .count()
    .withColumnRenamed("LoF", "LoF_Ch")
    .withColumnRenamed("GoF", "GoF_Ch")
    .join(chembl_trials, on=["targetId", "diseaseId"], how="left")
    .persist()
)

coloc_bnch = coloc_otgene.join(chembl, on=["targetId", "diseaseId"], how="inner")
coloc_bnch2 = coloc_bnch.join(
    sampleSize.select("right_study_id", "sampleSize"), on="right_study_id", how="left"
).persist()

withEvidence = (
    coloc_bnch2.withColumn(
        "ChEMBL",
        F.when(
            (F.col("GoF_Ch").isNotNull()) & (F.col("LoF_Ch").isNotNull()),
            F.lit("gof&lof"),
        )
        .when(
            (F.col("LoF_Ch").isNotNull()) & (F.col("GoF_Ch").isNull()),
            F.lit(F.lit("lof")),
        )
        .when(
            (F.col("GoF_Ch").isNotNull()) & (F.col("LoF_Ch").isNull()),
            F.lit(F.lit("gof")),
        ),
    )
    .withColumn(
        "Coherency_chembl",
        F.when(  ### there are cases of drug with gof&lof
            (F.col("protect").isNotNull()),
            F.when(
                (F.col("beta_assessed") == "gof"),
                F.when(
                    (F.col("GoF_Ch").isNotNull()) & (F.col("LoF_Ch").isNull()),
                    F.lit("coherent"),
                ).when(
                    (F.col("LoF_Ch").isNotNull()) & (F.col("GoF_Ch").isNull()),
                    F.lit("dispar"),
                ),
            ).when(
                (F.col("beta_assessed") == "lof"),
                F.when(
                    (F.col("GoF_Ch").isNotNull()) & (F.col("LoF_Ch").isNull()),
                    F.lit("dispar"),
                ).when(
                    (F.col("LoF_Ch").isNotNull()) & (F.col("GoF_Ch").isNull()),
                    F.lit("coherent"),
                ),
            ),
        ),
    )
    .join(target.select("id", "approvedSymbol"), target.id == F.col("targetId"), "left")
    .drop("id")
    .persist()
)

uniqueBetas = withEvidence.persist()

custom_schema = StructType(
    [
        StructField("_c0", StringType(), True),
        StructField("_c1", StringType(), True),
        StructField("_c2", StringType(), True),
        StructField("_c3", DecimalType(38, 37), True),
        # Add more fields as needed
    ]
)

### windows function for adjusting pvalue using benjamini-hochberg correction
window = Window.partitionBy("studyId").orderBy("pVal")
window2 = Window.partitionBy("studyId")

### read the tissue enrichent's file
tissueEnrichment = (
    spark.read.csv(
        "gs://ot-team/jroldan/Tissue_enrichment_results/",
        sep="\t",
        schema=custom_schema,
    )
    .selectExpr(
        "_c0 as studyId",
        "_c1 as tissueEnriched",
        "_c3 as pVal",
        #### instead of removing, keeping the ones with pvalues=0 to not perturb the rank
    )
    .withColumn(
        "pVal",
        F.when(F.col("pVal") == 0e-37, F.lit(1e-37)).otherwise(F.col("pVal")),
        ### ranking column
    )
    .withColumn(
        "index",
        F.row_number().over(window),
        ### total number of rows per study column
    )
    .withColumn(
        "length",
        F.last("index").over(window2),
        ### BH adjusted p values calculation
    )
    .withColumn("adjPVal", (F.col("pVal") * F.col("length")) / F.col("index"))
)

### read the file with matching tissue enriched with therapy areas
tissueEnrichTherAreas = (
    spark.read.csv(
        "gs://ot-team/jroldan/analysis/20231207_gwasTissueEnrrichedToTherapyAreas.csv/",
        sep=",",
        header=True,
    )
    .withColumnRenamed("tissue", "tissueEnriched")
    .withColumn(
        "tissueEnrichedTherapyAreas",
        F.array(F.col("TherapyArea"), F.col("Alternative")),
    )
    .persist()
)

### prepare the file of tissue enriched with their therapy areas
studyEnrichTherArea = (
    tissueEnrichment.join(
        tissueEnrichTherAreas.select("tissueEnriched", "tissueEnrichedTherapyAreas"),
        on="tissueEnriched",
        how="left",
        ## filter by p value <0.05
    )
    ### filter by adjusted p values below 0.05
    .filter(F.col("adjPVal") < 0.05)
    .groupBy("studyId")
    .agg(
        F.array_except(
            F.flatten(F.collect_set("tissueEnrichedTherapyAreas")), F.array(F.lit(None))
        ).alias("studyTherapyArea"),
        F.collect_set("tissueEnriched").alias("tissuesEnriched"),
    )
)

#### mapping right_bio_feature - disease
schema_rbf = StructType(
    [
        StructField("remove", StringType(), True),
        StructField("right_bio_feature", StringType(), True),
        StructField("name", StringType(), True),
        StructField("matchRBFToDisease", StringType(), True),
        # Add more fields as needed
    ]
)

### read the tissue enrichent's file
tissueToDisease = (
    spark.read.csv(
        "gs://ot-team/jroldan/analysis/jroldan_analysis_20240801_newcurationRBFToDisease_GOOD.csv/",
        header=True,
    )
    .select("name", "right_bio_feature", "matchRBFToDisease")
    .distinct()
)


uniqueBetas_location = (
    (
        target_membrane(spark, target, uniqueBetas).drop(
            "mb",
            "counted",
            "HPA_membrane",
            "HPA_secreted",
            "uniprot_membrane",
            "uniprot_secreted",
            "loc",
            "location_id",
            # "result",
            "loc_id",
        )
        # .join(targetType, on=["targetId"], how="left")
    )
    .join(tissueToDisease, on=["name", "right_bio_feature"], how="left")
    .withColumn("qtlTherapyArea", F.split(F.col("therapyArea"), ","))
    .join(studyEnrichTherArea, on="studyId", how="left")
    .withColumn(
        "taIntersectQtlDisease",
        F.array_intersect(F.col("qtlTherapyArea"), F.col("diseaseTherapyAreas")),
    )
    .withColumn(
        "taIntersectQtlTissueEnriched",
        F.array_intersect(F.col("qtlTherapyArea"), F.col("studyTherapyArea")),
    )
    .withColumn(
        "taNIntersectQtlDisease",
        F.when(
            F.size(
                F.array_intersect(F.col("qtlTherapyArea"), F.col("diseaseTherapyAreas"))
            )
            >= 1,
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "taNIntersectQtlTissueEnriched",
        F.when(
            F.size(
                F.array_intersect(F.col("qtlTherapyArea"), F.col("studyTherapyArea"))
            )
            >= 1,
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .join(targetType, on="approvedSymbol", how="left")
)


"""single column comparisons"""
dfs = {}
claves = []
columns = []
list_columns = [
    "right_bio_feature",
    "right_type",
    "beta_assessed",
    "right_study_id",
    "therapyArea",
]

for x in list_columns:
    a = [data[0] for data in uniqueBetas_location.select(x).distinct().collect()]
    for y in a:
        claves.append(y)
        columns.append(x)
#### remove the ".0" that is problematic
claves2 = [str(s).replace(".0", "") for s in claves]

for key, value in zip(claves2, columns):
    dfs[key] = value

"""combined column comparisons"""
list1 = []
list2 = []
cols = uniqueBetas_location.drop("Coherency_geneBurden").columns

### make list per column to be combined
study_column = [
    str(row.right_study_id)
    for row in uniqueBetas_location.filter(F.col("right_study_id").isNotNull())
    .select("right_study_id")
    .distinct()
    .collect()
]

qtl_type = [
    str(row.right_type)
    for row in uniqueBetas_location.filter(F.col("right_type").isNotNull())
    .select("right_type")
    .distinct()
    .collect()  ## there are "Nulls" in this column
]

biofeature_column = [
    str(row.right_bio_feature)
    for row in uniqueBetas_location.filter(F.col("right_bio_feature").isNotNull())
    .select("right_bio_feature")
    .distinct()
    .collect()
]

### make a dictionary with column combinations study & qtl
for x in study_column:
    for y in qtl_type:
        list1.append(x)
        list2.append(y)

### make a dictionary with column combinations biofeature & qtl
for x in biofeature_column:
    for y in qtl_type:
        list1.append(x)
        list2.append(y)

combined_lists = [list(pair) for pair in zip(list1, list2)]


## collect all columns in an array
### drop Coherency_geneBurden because it contains "coherent" and "dispar" words
def dict_comb_comp(combined_lists, study_column, qtl_type):
    new2 = {}
    for x, n in combined_lists:
        if x in study_column:
            if n in qtl_type:
                new2.update(
                    {f"{x}&{n}": {"right_type": {f"{n}": {"right_study_id": f"{x}"}}}}
                )
    return new2


new2 = dict_comb_comp(combined_lists, study_column, qtl_type)

df_string = (
    uniqueBetas_location.drop("Coherency_geneBurden")
    .select([F.col(col_name).cast("string").alias(col_name) for col_name in cols])
    .withColumn("array", F.array(cols))
)

"""make the dataframe with column combinations"""
df_string2 = df_string.select(
    ["*"]
    + [
        F.when(
            (F.col(i) == b)
            & (F.col(c) == t)
            & (F.col("Coherency_chembl") == "coherent"),
            F.lit("yes"),
        )
        .otherwise(F.lit("no"))
        .alias(a)
        for a, x in new2.items()
        for i, z in x.items()
        for b, y in z.items()
        for c, t in y.items()
        # print(a, i, b, c, t)
    ]
).persist()

delimiter = ","  ### to convert string of propagated traits to array

uniqueBetas_analyse = (
    df_string2.select(
        ["*"]
        + [  #### only doe
            F.when(
                ((f"{x}") == F.col(c)) & (F.col("Coherency_chembl") == "coherent"),
                F.lit("yes"),
            )
            .otherwise(F.lit("no"))
            .alias(f"{x}_{c}")
            for x, c in dfs.items()
        ]
    )
    .join(terminated, on=["targetId", "diseaseId"], how="left")
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
        F.when(F.col("maxClinPhase") >= 0, F.lit("yes")).otherwise(F.lit("no")),
    )
    .withColumn(
        "allPhases",
        F.when(F.col("maxClinPhase").isNotNull(), F.lit("yes")).otherwise(F.lit("no")),
    )
    .withColumn(
        "PhaseT",
        F.when(F.col("prediction").isNotNull(), F.lit("yes")).otherwise(F.lit("no")),
    )
    .withColumn(
        "coloc>_60",
        F.when(
            (F.col("coloc_h4") >= 0.60) & (F.col("Coherency_chembl") == "coherent"),
            F.lit("yes"),
        )
        .when(
            (F.col("coloc_h4") >= 0.60) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4") < 0.60) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4").isNull()) | (F.col("Coherency_chembl").isNull()),
            F.lit("no"),
        )
        .otherwise(F.lit("no")),
    )
    .withColumn(
        "coloc>_80",
        F.when(
            (F.col("coloc_h4") >= 0.80) & (F.col("Coherency_chembl") == "coherent"),
            F.lit("yes"),
        )
        .when(
            (F.col("coloc_h4") >= 0.80) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4") < 0.80) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4").isNull()) | (F.col("Coherency_chembl").isNull()),
            F.lit("no"),
        )
        .otherwise(F.lit("no")),
    )
    .withColumn(
        "coloc>_85",
        F.when(
            (F.col("coloc_h4") >= 0.85) & (F.col("Coherency_chembl") == "coherent"),
            F.lit("yes"),
        )
        .when(
            (F.col("coloc_h4") >= 0.85) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4") < 0.85) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4").isNull()) | (F.col("Coherency_chembl").isNull()),
            F.lit("no"),
        )
        .otherwise(F.lit("no")),
    )
    .withColumn(
        "coloc>_90",
        F.when(
            (F.col("coloc_h4") >= 0.90) & (F.col("Coherency_chembl") == "coherent"),
            F.lit("yes"),
        )
        .when(
            (F.col("coloc_h4") >= 0.90) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4") < 0.90) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4").isNull()) | (F.col("Coherency_chembl").isNull()),
            F.lit("no"),
        )
        .otherwise(F.lit("no")),
    )
    .withColumn(
        "coloc>_95",
        F.when(
            (F.col("coloc_h4") >= 0.95) & (F.col("Coherency_chembl") == "coherent"),
            F.lit("yes"),
        )
        .when(
            (F.col("coloc_h4") >= 0.95) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4") < 0.95) & (F.col("Coherency_chembl") == "dispar"),
            F.lit("no"),
        )
        .when(
            (F.col("coloc_h4").isNull()) | (F.col("Coherency_chembl").isNull()),
            F.lit("no"),
        )
        .otherwise(F.lit("no")),
    )
    .withColumn(
        "secreted",
        F.when(F.col("Nr_secreted") == 1, F.lit("yes"))
        .when(F.col("Nr_secreted") == 0, F.lit("no"))
        .otherwise(None),
    )
    .withColumn(
        "matchQtlTeTherArea",
        F.when(
            (F.col("taNIntersectQtlTissueEnriched") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "matchQtlDiseaseTherArea",
        F.when(
            (F.col("taNIntersectQtlDisease") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "matchQtlRBFrelevantToDisease",
        F.when(
            (F.col("matchRBFToDisease") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "doeCoherent",
        F.when(F.col("Coherency_chembl") == "coherent", F.lit("yes")).otherwise(
            F.lit("no")
        ),
    )
    .withColumn(
        "doeNotCoherent",
        F.when(F.col("Coherency_chembl") == "coherent", F.lit("yes")).otherwise(
            F.lit("no")
        ),
    )
    .withColumn(
        "doe&matchQtlTeTherArea",
        F.when(
            (F.col("Coherency_chembl") == "coherent")
            & (F.col("taNIntersectQtlTissueEnriched") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "doe&matchQtlDiseaseTherArea",
        F.when(
            (F.col("Coherency_chembl") == "coherent")
            & (F.col("taNIntersectQtlDisease") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "doe&matchQtlRBFrelevantToDisease",
        F.when(
            (F.col("Coherency_chembl") == "coherent")
            & (F.col("matchRBFToDisease") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "doeInflam",
        F.when(
            (F.col("coherency_chembl") == "coherent") & (F.col("isInflam") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "isInflam",
        F.when(F.col("isInflam") == "yes", F.lit("yes")).otherwise(F.lit("no")),
    )
    ### columns combinations for tissue enrichment and therapy areas
    .withColumn(
        "doeInflamMatchQtlTeTherArea",
        F.when(
            (F.col("Coherency_chembl") == "coherent")
            & (F.col("taNIntersectQtlTissueEnriched") == "yes")
            & (F.col("isInflam") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "doeInflamMatchQtlDiseaseTherArea",
        F.when(
            (F.col("Coherency_chembl") == "coherent")
            & (F.col("taNIntersectQtlDisease") == "yes")
            & (F.col("isInflam") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "doeInflamMatchQtlRBFrelevantToDisease",
        F.when(
            (F.col("Coherency_chembl") == "coherent")
            & (F.col("matchRBFToDisease") == "yes")
            & (F.col("isInflam") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    #### columns considering only inflam and tissue/therapy areas
    .withColumn(
        "inflamMatchQtlTeTherArea",
        F.when(
            (F.col("taNIntersectQtlTissueEnriched") == "yes")
            & (F.col("isInflam") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "inflamMatchQtlDiseaseTherArea",
        F.when(
            (F.col("taNIntersectQtlDisease") == "yes") & (F.col("isInflam") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .withColumn(
        "inflamMatchQtlRBFrelevantToDisease",
        F.when(
            (F.col("matchRBFToDisease") == "yes") & (F.col("isInflam") == "yes"),
            F.lit("yes"),
        ).otherwise(F.lit("no")),
    )
    .drop("targetType")  ### remove nonsense columns for this
    .filter(
        (F.col("coherency_chembl").isNotNull())
        & (F.col("name") != "COVID-19")
        # | (F.col("name").isNull())
    )
    .drop(
        "array",
        "remove",
        "Phase0",
        "allPhases",
        "Terminated",
        "propagatedTraits",
        "count",
        "clinicalStatus",
        "prediction",
        "allPhases",
        "propagatedTraits",
        "targetType",
    )
    .repartition(100)
    .persist()
)

import re
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from pyspark.sql.types import *

result = []
result_st = []
result_ci = []
array2 = []

# Initialize an empty list to store the results
results = []


def convertTuple(tup):
    st = ",".join(map(str, tup))
    return st


key_list = [
    "coloc>_60",
    "coloc>_80",
    "coloc>_90",
    "coloc>_95",
    "secreted",
    "matchQtlTeTherArea",
    "matchQtlDiseaseTherArea",
    "matchQtlRBFrelevantToDisease",
    "doeCoherent",
    "doeInflam",
    "doeInflamMatchQtlTeTherArea",
    "doeInflamMatchQtlRBFrelevantToDisease",
    "inflamMatchQtlTeTherArea",
    "inflamMatchQtlDiseaseTherArea",
    "inflamMatchQtlRBFrelevantToDisease",
    "doe&matchQtlTeTherArea",
    "doe&matchQtlDiseaseTherArea",
    "doe&matchQtlRBFrelevantToDisease",
]
value_list = [
    "coloc",
    "coloc",
    "coloc",
    "coloc",
    "location",
    "rght_tissue",
    "right_tissue",
    "right_tissue",
    "doeGlobal",
    "inflammation",
    "inflammation&DoE",
    "inflammation&DoE&right_tissue",
    "inflammation&DoE&right_tissue",
    "inflammation&DoE&right_tissue",
    "inflammation&DoE&right_tissue",
    "doe&right_tissue",
    "doe&right_tissue",
    "doe&right_tissue",
]


def comparisons_df(dfs, key_list, value_list) -> list:
    """Return list of all comparisons to be used in the analysis"""

    toAnalysis = []
    toAnalysis = [[key, x] for key, x in dfs.items()]
    if len(key_list) == len(value_list):
        # Pair the elements from the two lists together
        toAnalysis.extend(
            [[key, custom_name] for key, custom_name in zip(key_list, value_list)]
        )
    else:
        print("The input lists have different lengths and cannot be paired.")

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
print("created full_data and lists")


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
    wComparison = Window.partitionBy(comparisonColumn)
    wPrediction = Window.partitionBy(predictionColumn)
    wPredictionComparison = Window.partitionBy(comparisonColumn, predictionColumn)

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
            + comparisonType
            + "_"
            + predictionColumn
            + ".parquet"
        )
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
            + comparisonType
            + "_"
            + predictionColumn
            + ".parquet"
        )
    )
    path = (
        today_date
        + "_"
        + "analysis/"
        + data
        # + "_propagated"
        + "/"
        + comparisonColumn
        + "_"
        + comparisonType
        + "_"
        + predictionColumn
        + ".parquet"
    )
    print(path)
    ### making analysis
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
    # print(array1)
    total = np.sum(array1)
    res_npPhaseX = np.array(array1, dtype=int)
    resX = convertTuple(fisher_exact(res_npPhaseX, alternative="two-sided"))
    resx_CI = convertTuple(
        odds_ratio(res_npPhaseX).confidence_interval(confidence_level=0.95)
    )

    result_st.append(resX)
    result_ci.append(resx_CI)
    (rs_result, rs_ci) = relative_success(array1)
    # print(total, resX, resx_CI, rs_result, rs_ci)

    results.append(
        [
            comparisonType,
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
            path,
        ]
    )
    return results


c = datetime.now()
c.strftime("%H:%M:%S")
# print(c)

print("start doing aggregations and writing")
today_date = str(date.today())
# today_date = "2024-07-17"
aggSetups_original = comparisons_df(dfs, key_list, value_list)
listado = []

print("starting with aggregations at", c)

for row in aggSetups_original:
    results = aggregations_original(
        uniqueBetas_analyse, "propagated", listado, *row, today_date
    )
print("aggregations and analysis done")
print("building dataframe")
"""
df = pd.DataFrame(
    results,
    columns=[
        "group",
        "comparison",
        "phase",
        "oddsRatio",
        "pValue",
        "lowerInterval",
        "upperInterval",
        "total",
        "values",
        "relSuccess",
        "rsLower",
        "rsUpper",
        "path",
    ],
)
"""
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

# Convert list of lists to DataFrame
df = spreadSheetFormatter(spark.createDataFrame(results, schema=schema))
df.toPandas().to_csv(f"gs://ot-team/jroldan/analysis/{today_date}_colocDoEanalysis.csv")

# 1# Make a list of variant of interest (Sequence ontology terms) to subset data of interest.

### Bear in mind that SO works with ontology structure as: SO:XXXXXX, but databases has the SO as: SO_XXXXXX

var_filter_lof = [
    ### High impact variants https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
    "SO_0001589",  ## frameshit_variant
    "SO_0001587",  ## stop_gained
    "SO_0001574",  ## splice_acceptor_variant
    "SO_0001575",  ## splice_donor_variant
    "SO_0002012",  ## start_lost
    "SO_0001578",  ## stop_lost
    "SO_0001893",  ## transcript_ablation
    # "SO:0001889", ## transcript_amplification ## the Only HIGH impact that increase protein.
]

gof = ["SO_0002053"]
lof = ["SO_0002054"]


## Building Sequence Ontology
so_path = (
    "/Users/juanr/Desktop/Target_Engine/data_download/sequenceOntology_20221118.csv"
)
so_ontology = spark.read.csv(so_path, header=True)
building = (
    so_ontology.select(F.col("Accession"), F.col("Parents"))
    .withColumn("Parentalind", F.split(F.col("Parents"), ","))
    .withColumn("Parentalind", F.explode_outer("Parentalind"))
    .groupBy("Parentalind")
    .agg(F.collect_list(F.col("Accession")).alias("childrens"))
    .join(so_ontology, F.col("Parentalind") == so_ontology.Accession, "right")
)


### Load evidence datasources downloaded in January 2023:

otgenetics_evidence_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=ot_genetics_portal"
otgenetics = spark.read.parquet(otgenetics_evidence_path)
gene_burden_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=gene_burden"
gene_burden = spark.read.parquet(gene_burden_path)
eva_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=eva"
eva_germline = spark.read.parquet(eva_path)
eva_somatic_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=eva_somatic"
eva_somatic = spark.read.parquet(eva_somatic_path)
gel_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=genomics_england"
gel = spark.read.parquet(gel_path)
g2p_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=gene2phenotype"
g2p = spark.read.parquet(g2p_path)
uniprot_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=uniprot_literature"
uniprot = spark.read.parquet(uniprot_path)
uniprotvar_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=uniprot_variants"
uniprotvar = spark.read.parquet(uniprotvar_path)
orphanet_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=orphanet"
orphanet = spark.read.parquet(orphanet_path)
clingen_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=clingen"
clingen = spark.read.parquet(clingen_path)
cgc_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=cancer_gene_census"
cgc = spark.read.parquet(cgc_path)
intogen_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=intogen"
intogen = spark.read.parquet(intogen_path)
impc_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=impc"
impc = spark.read.parquet(impc_path)
chembl_evidences = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/evidence/sourceId=chembl/"
chembl = spark.read.parquet(chembl_evidences)


## others
target_path = (
    "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/targets/"
)
target = spark.read.parquet(target_path)
disease_path = (
    "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/diseases/"
)
diseases = spark.read.parquet(disease_path)
dis_name = diseases.select("id", "name")
indication_path = (
    "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/indication/"
)
indication = spark.read.parquet(indication_path)
## drug_path="/Users/juanr/Desktop/Target_Engine/downloadedEvidencesJanuary/molecule/"
## drug=spark.read.parquet(drug_path)
mecact_path = "/Users/juanr/Desktop/Target_Engine/DownloadFebruary_Release23.02/mechanismOfAction/"
mecact = spark.read.parquet(mecact_path)

#### GENE BURDEN

### We manually annotated those studies using LoF or PTV variants

burden_lof_path = (
    "/Users/juanr/Desktop/Target_Engine/Conteo_estudios_geneBurden_20230117.csv"
)
burden_lof = spark.read.csv(burden_lof_path, header=True)
burden_lof = burden_lof.withColumnRenamed("statisticalMethodOverview", "stMethod")

#### Para gene burden la funcion no tiene que hacer un filtrado de variantes

### EVA/ClinVar

##- Manually annotate which are the clinicalSignificances meaningfull: pathogenic, risk factor, protective

clinSign_germline_path = "/Users/juanr/Desktop/Target_Engine/eva_clinSig_20230117.csv"
clinSign_somatic_path = (
    "/Users/juanr/Desktop/Target_Engine/eva_somatic_clinSig_20230117.csv"
)

clinSign_germline = spark.read.csv(clinSign_germline_path, header=True)
clinSign_germline = clinSign_germline.withColumnRenamed(
    "clinicalSignificances", "significances"
)
clinSign_somatic = spark.read.csv(clinSign_somatic_path, header=True)
clinSign_somatic = clinSign_somatic.withColumnRenamed(
    "clinicalSignificances", "significances"
)

##-  Transform array of clinicalSignificances into Strings to check them.

eva_somatic_toAsses = eva_somatic.withColumn(
    "clinicalSignificances", F.concat_ws(",", F.col("clinicalSignificances"))
)

eva_germline_toAsses = eva_germline.withColumn(
    "clinicalSignificances", F.concat_ws(",", F.col("clinicalSignificances"))
)


## annotate TSG/oncogene/bivalent using 'hallmarks.attributes'

oncotsg_list = [
    "TSG",
    "oncogene",
    "Oncogene",
    "oncogene",
    "oncogene,TSG",
    "TSG,oncogene",
    "fusion,oncogene",
    "oncogene,fusion",
]

oncolabel = (
    target.select(
        "id", "approvedSymbol", F.explode_outer(F.col("hallmarks.attributes"))
    )
    .select("id", "approvedSymbol", "col.description")
    .filter(F.col("description").isin(oncotsg_list))
    .groupBy("id", "approvedSymbol")
    .agg(F.collect_set("description").alias("description"))
    .withColumn("description_splited", F.concat_ws(",", F.col("description")))
    .withColumn(
        "TSorOncogene",
        F.when(
            (
                F.col("description_splited").rlike("ncogene")
                & F.col("description_splited").rlike("TSG")
            ),
            F.lit("bivalent"),
        )
        .when(F.col("description_splited").rlike("ncogene(\s|$)"), F.lit("oncogene"))
        .when(F.col("description_splited").rlike("TSG(\s|$)"), F.lit("TSG"))
        .otherwise(F.lit("noEvaluable")),
    )
    .withColumnRenamed("id", "target_id")
)

#### rlike('('+Keywords+')(\s|$)'


### Hacer el join del actionType con el chembl para sacar los mecanismos de accion.
inhibitors = [
    "RNAI INHIBITOR",
    "NEGATIVE MODULATOR",
    "NEGATIVE ALLOSTERIC MODULATOR",
    "ANTAGONIST",
    "ANTISENSE INHIBITOR",
    "BLOCKER",
    "INHIBITOR",
    "DEGRADER",
    "INVERSE AGONIST",
    "ALLOSTERIC ANTAGONIST",
]

activators = [
    "PARTIAL AGONIST",
    "ACTIVATOR",
    "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR",
    "AGONIST",
    "SEQUESTERING AGENT",
]

columnas = ["activator", "inhibitor"]
both = activators + inhibitors

actiontype2 = (
    mecact.select("chemblIds", "actionType", "mechanismOfAction", "targets")
    .select(
        F.explode_outer("chemblIds").alias("drugId2"),
        "actionType",
        "mechanismOfAction",
        "targets",
    )
    .select(
        F.explode_outer("targets").alias("targetId2"),
        "drugId2",
        "actionType",
        "mechanismOfAction",
    )
    .dropDuplicates()
)


chembl1 = chembl.select(
    "targetId", "drugId", "diseaseId", "clinicalPhase", "diseaseFromSourceId"
)
chembl2 = (
    chembl1.join(
        actiontype2,
        (actiontype2.drugId2 == F.col("drugId"))
        & (actiontype2.targetId2 == F.col("targetId")),
        "left",
    ).drop("targetId2", "drugId2")
    ###.dropDuplicates()
    .withColumn(
        "twoCategories_new",
        F.when(F.col("actionType").isin(inhibitors), F.lit("inhibitor"))
        .when(F.col("actionType").isin(activators), F.lit("activator"))
        .otherwise(F.lit("noEvaluable")),
    )
)

chembl3 = (
    chembl2.filter(F.col("twoCategories_new") != "noEvaluable")
    .groupBy("targetId", "diseaseId")
    .pivot("twoCategories_new")
    .agg(F.count("targetId"))
)

chembl4 = chembl3.select(
    "targetId",
    "diseaseId",
    ##'clinicalPhase',
    *(F.col(c).cast("int").alias(c) for c in columnas)
).withColumn(
    "coherency",
    F.when(
        ((F.col("activator").isNotNull()) & (F.col("inhibitor").isNotNull())),
        F.when(
            (F.col("activator")) - (F.col("inhibitor")) != (F.col("activator")),
            F.lit("dispar"),
        ),
    ),
)

### Join all datasets

dfs = [
    otgenetics,
    gene_burden,
    eva_germline,
    eva_somatic,
    g2p,
    orphanet,
    cgc,
    intogen,
    impc,
    chembl,
]

all = dfs[0]
for df in dfs[1:]:
    all = all.unionByName(df, allowMissingColumns=True)
all.count()

#### 20230203 ###
prueba_assessment = (
    all.withColumn(
        "beta", F.col("beta").cast("float")
    )  ## from ot genetics & gene burden
    .withColumn(
        "OddsRatio", F.col("OddsRatio").cast("float")
    )  ## from ot genetics & gene burden
    .withColumn(
        "clinicalSignificances", F.concat_ws(",", F.col("clinicalSignificances"))
    )  ### from eva
    .withColumn(
        "exploded", F.explode_outer(F.col("mutatedSamples"))
    )  ### para cgc e intogen
    .withColumn(
        "variantConsequence", F.col("exploded.functionalConsequenceId")
    )  ### para cgc e intogen
    ### .withColumn('numberSamplesSameMutationType',F.col('exploded.numberSamplesWithMutationType'))### para cgc e intogen
    .withColumn(
        "mutatedSamplesVariantInfo",
        F.coalesce(F.col("mutatedSamples.functionalConsequenceId"), F.array()),
    )  ### para cgc e intogen
    .join(oncolabel, oncolabel.target_id == F.col("targetId"), "left")  ### para cgc
    .join(
        burden_lof, burden_lof.stMethod == F.col("statisticalMethodOverview"), "left"
    )  ### para gene burden
    .join(
        actiontype2,  ## para chembl
        (actiontype2.drugId2 == F.col("drugId"))
        & (actiontype2.targetId2 == F.col("targetId")),
        "left",
    )
    ##.drop('targetId2','drugId2')
    ###.dropDuplicates()
    .withColumn(
        "Assessment",
        #### Ot_genetics Portal ### updated to include the coloc+gwas analysis
        F.when(
            F.col("datasourceId") == "ot_genetics_portal",
            F.when(  ### label 14 evidences that are contradictory
                (
                    (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002315")
                    & (F.col("variantFunctionalConsequenceId").isin(var_filter_lof))
                ),
                F.lit("dispar"),
            )
            ### evidences with gwas+coloc increased expression without +var_lof
            .when(
                (
                    (F.col("beta").isNull())
                    & (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002315")
                    & (
                        F.col("variantFunctionalConsequenceId").isin(var_filter_lof)
                        == False
                    )
                ),
                F.when((F.col("OddsRatio") > 1), F.lit("GoF_risk")).when(
                    (F.col("OddsRatio") < 1), F.lit("GoF_protect")
                ),
            ).when(
                (
                    (F.col("oddsRatio").isNull())
                    & (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002315")
                    & (
                        F.col("variantFunctionalConsequenceId").isin(var_filter_lof)
                        == False
                    )
                ),
                F.when((F.col("beta") < 0), F.lit("GoF_protect")).when(
                    (F.col("beta") > 0), F.lit("GoF_risk")
                ),
            )
            ### evidences with coherent Gwas-coloc + var_lof
            .when(
                (
                    (F.col("beta").isNull())
                    & (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002316")
                    & (F.col("variantFunctionalConsequenceId").isin(var_filter_lof))
                ),
                F.when((F.col("OddsRatio") > 1), F.lit("LoF_risk")).when(
                    (F.col("OddsRatio") < 1), F.lit("LoF_protect")
                ),
            ).when(
                (
                    (F.col("oddsRatio").isNull())
                    & (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002316")
                    & (F.col("variantFunctionalConsequenceId").isin(var_filter_lof))
                ),
                F.when((F.col("beta") < 0), F.lit("LoF_protect")).when(
                    (F.col("beta") > 0), F.lit("LoF_risk")
                ),
            )
            ### evidences with colo+Gwas data but not variants
            .when(
                (
                    (F.col("beta").isNull())
                    & (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002316")
                    & (
                        F.col("variantFunctionalConsequenceId").isin(var_filter_lof)
                        == False
                    )
                ),
                F.when((F.col("OddsRatio") > 1), F.lit("LoF_risk")).when(
                    (F.col("OddsRatio") < 1), F.lit("LoF_protect")
                ),
            ).when(
                (
                    (F.col("oddsRatio").isNull())
                    & (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002316")
                    & (
                        F.col("variantFunctionalConsequenceId").isin(var_filter_lof)
                        == False
                    )
                ),
                F.when((F.col("beta") < 0), F.lit("LoF_protect")).when(
                    (F.col("beta") > 0), F.lit("LoF_risk")
                ),
            )
            ### evidences with coherent non/inconclusive gwas+coloc + var_lof
            .when(
                (
                    (F.col("beta").isNull())
                    & (
                        (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002314")
                        | (F.col("variantFunctionalConsequenceFromQtlId").isNull())
                    )
                    & (F.col("variantFunctionalConsequenceId").isin(var_filter_lof))
                ),
                F.when((F.col("OddsRatio") > 1), F.lit("LoF_risk")).when(
                    (F.col("OddsRatio") < 1), F.lit("LoF_protect")
                ),
            )
            .when(
                (
                    (F.col("oddsRatio").isNull())
                    & (
                        (F.col("variantFunctionalConsequenceFromQtlId") == "SO_0002314")
                        | (F.col("variantFunctionalConsequenceFromQtlId").isNull())
                    )
                    & (F.col("variantFunctionalConsequenceId").isin(var_filter_lof))
                ),
                F.when((F.col("beta") < 0), F.lit("LoF_protect")).when(
                    (F.col("beta") > 0), F.lit("LoF_risk")
                ),
            )
            .otherwise(F.lit("noEvaluable")),
        )  ### son tambien no data las que tiene riesgo pero no tienen LoF
        #### Gene burden
        .when(
            F.col("datasourceId") == "gene_burden",
            ### .filter(F.col('variantType').isin(var_filter))
            F.when(
                (
                    (F.col("whatToDo") == "get")
                    & (F.col("beta").isNull())
                    & (F.col("OddsRatio") > 1)
                ),
                F.lit("LoF_risk"),
            )
            .when(
                (
                    (F.col("whatToDo") == "get")
                    & (F.col("beta").isNull())
                    & (F.col("OddsRatio") < 1)
                ),
                F.lit("LoF_protect"),
            )
            .when(
                (
                    (F.col("whatToDo") == "get")
                    & (F.col("OddsRatio").isNull())
                    & (F.col("beta") > 0)
                ),
                F.lit("LoF_risk"),
            )
            .when(
                (
                    (F.col("whatToDo") == "get")
                    & (F.col("OddsRatio").isNull())
                    & (F.col("beta") < 0)
                ),
                F.lit("LoF_protect"),
            )
            .otherwise(
                F.lit("noEvaluable")
            ),  ### son tambien no data las que tiene riesgo pero no se ensayan LoF o PT
        )
        #### Eva_germline
        .when(
            F.col("datasourceId") == "eva",
            #### .filter(F.col('variantFunctionalConsequenceId').isin(var_filter_lof))
            F.when(
                (
                    ## (F.col('clinicalSignificances')!='likely pathogenic') &
                    (F.col("variantFunctionalConsequenceId").isin(var_filter_lof))
                    & F.col("clinicalSignificances").rlike("(pathogenic)$")
                ),
                F.lit("LoF_risk"),
            )
            .when(
                (
                    F.col("clinicalSignificances").contains("protective")
                    & F.col("variantFunctionalConsequenceId").isin(var_filter_lof)
                ),
                F.lit("LoF_protect"),
            )
            .otherwise(
                F.lit("noEvaluable")
            ),  ### Son todas aquellas que tenen info pero no son patogenicas/protective  + LoF
        )
        #### Eva_somatic
        .when(
            F.col("datasourceId") == "eva_somatic",
            #### .filter(F.col('variantFunctionalConsequenceId').isin(var_filter_lof))
            F.when(
                (
                    ##(F.col('clinicalSignificances')!='likely pathogenic') &
                    (F.col("variantFunctionalConsequenceId").isin(var_filter_lof))
                    & F.col("clinicalSignificances").rlike("(pathogenic)$")
                ),
                F.lit("LoF_risk"),
            )
            .when(
                (
                    F.col("clinicalSignificances").contains("protective")
                    & F.col("variantFunctionalConsequenceId").isin(var_filter_lof)
                ),
                F.lit("LoF_protect"),
            )
            .otherwise(
                F.lit("noEvaluable")
            ),  ### Son todas aquellas que tenen info pero no son patogenicas/protective  + LoF
        )
        #### G2P
        .when(
            F.col("datasourceId") == "gene2phenotype",
            F.when(
                F.col("variantFunctionalConsequenceId") == "SO_0002317",
                F.lit("LoF_risk"),
            )  ### absent gene product
            .when(
                F.col("variantFunctionalConsequenceId") == "SO_0002315",
                F.lit("GoF_risk"),
            )  ### increased gene product level
            .otherwise(F.lit("noEvaluable")),
        )
        #### Orphanet
        .when(
            F.col("datasourceId") == "orphanet",
            F.when(
                F.col("variantFunctionalConsequenceId") == "SO_0002054",
                F.lit("LoF_risk"),
            )  ### Loss of Function Variant
            .when(
                F.col("variantFunctionalConsequenceId") == "SO_0002053",
                F.lit("GoF_risk"),
            )  ### Gain_of_Function Variant
            .otherwise(F.lit("noEvaluable")),
        )
        #### CGC
        .when(
            F.col("datasourceId") == "cancer_gene_census",
            F.when(F.col("TSorOncogene") == "oncogene", F.lit("GoF_risk"))
            .when(F.col("TSorOncogene") == "TSG", F.lit("LoF_risk"))
            .when(F.col("TSorOncogene") == "bivalent", F.lit("bivalent_risk"))
            .otherwise(
                F.when(
                    F.arrays_overlap(
                        F.col("mutatedSamples.functionalConsequenceId"),
                        F.array([F.lit(i) for i in (var_filter_lof)]),
                    ),
                    F.lit("LoF_risk"),
                ).otherwise(F.lit("noEvaluable"))
            ),
        )  #### Aqui asumimos que todo lo que esta incluido da riesgo, pero solo podemos dar LoF porque ya no tienen dato de TSG/oncogen
        #### intogen
        .when(
            F.col("datasourceId") == "intogen",
            F.when(
                F.arrays_overlap(
                    F.col("mutatedSamples.functionalConsequenceId"),
                    F.array([F.lit(i) for i in (gof)]),
                ),
                F.lit("GoF_risk"),
            )
            .when(
                F.arrays_overlap(
                    F.col("mutatedSamples.functionalConsequenceId"),
                    F.array([F.lit(i) for i in (lof)]),
                ),
                F.lit("LoF_risk"),
            )
            .otherwise(F.lit("noEvaluable")),
        )
        #### impc
        .when(
            F.col("datasourceId") == "impc",
            F.when(F.col("diseaseId").isNotNull(), F.lit("KO_risk")).otherwise(
                F.lit("noEvaluable")
            ),
        )
        ### chembl
        .when(
            F.col("datasourceId") == "chembl",
            F.when(F.col("actionType").isin(inhibitors), F.lit("LoF_protect"))
            .when(F.col("actionType").isin(activators), F.lit("GoF_protect"))
            .otherwise(F.lit("noEvaluable")),
        ),
    )
    ### Homogenizar para contar todos los datos juntos:
    .withColumn(
        "homogenized",
        F.when(F.col("Assessment") == "KO_risk", F.lit("LoF_risk")).otherwise(
            F.col("Assessment")
        ),
    )
    .withColumn(
        "tendency",
        F.when(F.col("homogenized").contains("risk"), F.lit("Risk"))
        .when(F.col("homogenized").contains("protect"), F.lit("Protect"))
        .otherwise(F.lit("noEvaluable")),
    )
    .withColumn(
        "variation",
        F.when(F.col("homogenized").contains("LoF"), F.lit("LoF"))
        .when(F.col("homogenized").contains("GoF"), F.lit("GoF"))
        .otherwise(F.lit("noEvaluable")),
    )
)

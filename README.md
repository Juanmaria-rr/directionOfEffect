# Direction of Effect 

## Introduction

The direction of effect (DoE) in genetics refers to the observed relationship between a genetic variant and a specific trait or disease. It indicates whether a particular genetic variant increases (risk) or decreases (protect) the likelihood or severity of the trait or disease.

In this repository, we describe how to systematically address the direction of effect of targets on traits, using multiple sources of Open Target Evidences.

## Methodology

To assess the DoE we define two different categories with two outputs: 

    - ** variantEffect **: depending on the datasource, the variant/ level of expression/ drug mechanism is categorised in "loss of function (LoF)" or "gain of function (GoF)". At genetic level (i.e. OT genetics), only Loss of function variants are taken into account. For gene burden studies we manually annotated the valuable studies (including loss of function, protein truncating and/or damaging variants). At expression levels, annotated decreased/increased levels are translated to Loss /Gain of function, respectively. For drugs, inhibitors and activators are LoF and GoF, respectively.

    - ** directionOnTrait **: the direction is evaluated per datasource, translating the information into "risk" or "protect". Datasources with GWAS/Coloc information or collapsing analysis such as OT genetics and Gene burden provides the statistics Odds ratio and Beta values for evaluating "risk/protective". For clinical genetics datasources as clinVar, we only focus on pathogenic, protective or likely pathogenic variants. Cancer datasources or IMPC are translated to "risk" and ChEMBL drugs are all protective.

The final outputs are two columns, one per category. The combination of this content will give arise into the four posibilities of direction of effect. 


## Relevance in Drug Discovery and Development

Understanding the direction of effect enables the following:

### 1. Target Identification

Identifying potential drug targets is facilitated by considering the direction of effect. If a genetic variant with a given effect is associated with a specific disease, targeting the corresponding gene or protein becomes a promising approach for developing therapies to mitigate the disease's effects.

### 2. Drug Target Validation

The direction of effect aids in validating potential drug targets through the following steps:

- Perturb the target gene or protein in preclinical models or clinical trials.
- If mimicking the protective effect of a genetic variant by manipulating the target leads to positive outcomes, it suggests that developing drugs modulating the target could be therapeutically beneficial.

### 3. Personalized Medicine

The direction of effect guides the development of personalized medicine approaches by considering individual genetic profiles:

- Identify whether genes are associated with protective or risk effects for specific diseases or drug responses.
- Guide treatment decisions, including drug selection and dosage adjustments, tailored to each individual's genetic profile.

### 4. Biomarker Discovery

Discovering biomarkers is facilitated by considering the direction of effect:

- Examine genetic variants associated with risk or protective effects.
- Identify biomarkers predicting disease risk or treatment response, enabling more targeted and effective treatment strategies.


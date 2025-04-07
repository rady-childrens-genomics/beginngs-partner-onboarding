# BeginNGS® Partner Onboarding

# Overview

BeginNGS® (Begin Newborn Genome Sequencing) is a novel health care delivery system designed to screen newborns for genetic diseases — and connect their doctors with effective interventions (https://radygenomics.org/begin-ngs-newborn-sequencing/). BeginNGS is being undertaken by an open public-private consortium. The scope of BeginNGS is ambitious - to screen the genome sequences of all healthy or sick newborns or infants worldwide for all (~2,000) severe, childhood-onset genetic diseases (SCGD) that have effective therapies, provide confirmatory testing for screen positives, and provide guidance that equips physicians to initiate specific therapeutic interventions at or before symptom onset.

BeginNGS has quite a few activities. This GitHub repository relates to one of the more ambitious of these - federated queries of many large genomic sequencing datasets with large sets of variants that are believed causal for each BeginNGS SCGD and subsequent modeling of the genetic architecture of each SCGD across the lifespan and in most genetic ancestry groups (Kingsmore et al. Prequalification of genome-based newborn screening for severe childhood genetic diseases through federated training based on purifying hyperselection. Am J Hum Genet. In Press).

The narrow rationale for BeginNGS query federation is to overcome a major impediment to genome-based newborn screening (gNBS) - namely imprecision due to variants classified as pathogenic (P) or likely pathogenic (LP) that are not SCGD causal. Federated training predicated on purifying hyperselection provides a general framework to attain high precision in population screening. Federated training across many biobanks and clinical trials can provide a privacy-preserving mechanism for qualification of gNBS in diverse genetic ancestries.

The broader goal is to refine understanding of the genetic basis of the natural history of each SCGD in each genetic ancestry including prevalence, penetrance, expressivity, locus heterogeneity, variant (diplotype) heterogeneity, and untreated and treated outcomes across the lifespan. Examples of such analyses of the UK Biobank cohort are described in Kingsmore et al. Am J Hum Genet. In Press.

# Input Overview

- [Variants of interest](data/variants_of_interest_20231108.csv) (currently BeginNGS v2, 53,855 P and LP variants that map to 342 genes, 412 SCGD, and 1,603 SCGD therapeutic interventions) is normally encapsulated as a fixed resource inthe UDF, but can be implemented as an parameter. This is pre-annotated with consequence and population frequency information, but only chr-pos-ref-alt is used for the query itself
- [Blocklist](data/blocklist_20240329.csv) - entries classified as NSDCC (non-severe disease causing in childhood)
- [MOI](data/moi_20240805.txt) - mode of inheritance information
 
## Federated Query

![figure-run-federated-queries--UDF](https://github.com/user-attachments/assets/e38954d1-7e27-41c9-a461-8737e75462ce)


## Implementation for non-TileDB Users

The high level approach is:
1. Apply blocklist to variants of interest (if using blocklist)
2. Obtain vcf genotypes for variants of interest, merge on chr,pos,ref,alt
3. Classify compound hets, looking for co-occurring hets in sample/gene groups 
4. Merge sample metadata to get sex, use only consented subjects
5. Compute positive_genotypes using MOI rules
6. For each gene/subject grouping, compose a concatenated string of observed hit variants (diplotypes)
   
![BeginNGS-implementatio-for-non-TileDB-users](https://github.com/user-attachments/assets/45ccb8e4-ad59-4a72-aa63-265d3e749779)

### Source Code

Python source code for the BeginNGS query (both require a TileDB variant store)

[src/genotype_count_dispatch.py](src/genotype_count_dispatch.py) contains the mechanics of the BeginNGS query used in the federated implementation

[src/nbs_query_tiledb_repo.ipynb](src/nbs_query_tiledb_repo.ipynb) contains a more readable notebook-based implementation the BeginNGS query

## 1000 Genomes Validation Test Output

The following genotype and diplotype output was obtained by running the query on 3,202 high-coverage (Illumina 30X) DRAGEN 1000 Genomes Project samples (PMID:36055201, DOI: 10.1016/j.cell.2022.08.004) aligned to GRCh38 with Illumina DRAGEN v.4.2.7 (PMID: 39455800 DOI: 10.1038/s41587-024-02382-1. This can be used to validate external implementations.

The relevant 4.2.7 DRAGEN 1000 Genomes VCFs can be found in this Amazon S3 bucket:
`aws s3 ls --no-sign-request s3://1000genomes-dragen-v4-2-7/data/individuals/hg38_alt_masked_graph_v3/`

Query on 1kg with no blocklist

- [benchmarks/onekg_noblocklist_diplotypes.csv](benchmarks/onekg_noblocklist_diplotypes.csv)
- [benchmarks/onekg_noblocklist_genotypes.csv](benchmarks/onekg_noblocklist_genotypes.csv)

Query on 1kg with blocklist applied

- [benchmarks/onekg_diplotypes.csv](benchmarks/onekg_diplotypes.csv)
- [benchmarks/onekg_genotypes.csv](benchmarks/onekg_genotypes.csv)

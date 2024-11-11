# BeginNGS Partner Onboarding 

# Input Overview

- Variants of interest (currently 53,855 variants) is normally encapsulated as a fixed resource inthe UDF, but can be implemented as an parameter. This is pre-annotated with consequence and population frequency information, but only chr-pos-ref-alt is used for the query itself. 

## Federated Query

<img width="1589" alt="federated" src="https://github.com/user-attachments/assets/88c87f9f-b0a3-4754-8932-d6733adb5490">


## Implementation for non-TileDB Users
<img width="6930" alt="beginngs_flow" src="beginngs_flow.png">

### Source Code

Python source code for the BeginNGS query (both require a TileDB variant store)

[src/genotype_count_dispatch.py](src/genotype_count_dispatch.py) contains the mechanics of the BeginNGS query used in the federated implementation

[src/nbs_query_tiledb_repo.ipynb](src/nbs_query_tiledb_repo.ipynb) contains a more readable notebook-based implementation the BeginNGS query


## 1000 Genomes Validation Test Output

The following genotype and diplotype output was obtained by running the query on 3202 high-coverage DRAGEN 1000 Genomes samples. This can be used to validate external implementations.

- [benchmarks/onekg_diplotypes.csv](benchmarks/onekg_diplotypes.csv)
- [benchmarks/onekg_genotypes.csv](benchmarks/onekg_genotypes.csv)

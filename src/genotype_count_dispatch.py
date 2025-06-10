#!/usr/bin/env python

import logging
import math
import numpy as np
import pyarrow as pa
import pandas as pd
import tiledb
import tiledb.cloud
import tiledb.cloud.vcf
import tiledbvcf
import time

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import List, Optional
from tiledb.cloud.utilities import get_logger
from tiledb.cloud.vcf.vcf_toolbox import df_transform
from tiledb.cloud.vcf.vcf_toolbox.annotate import _annotate
from tiledb.cloud.vcf.vcf_toolbox.annotate import _annotate
from vcf_federated.filestore import read_lines
from tiledb.cloud.utilities import run_dag


def zygosity(gt: np.ndarray) -> str:
    """
    Convert a genotype array to a zygosity string classification.

    The zygosity is determined based on the number and values of alleles:
    - MISSING: Empty genotype or all alleles are -1 (missing)
    - HEMI: Single non-reference allele (hemizygous, typically X/Y chromosomes)
    - HOM_REF: All alleles are reference (0)
    - HOM_ALT: Multiple identical alternate alleles
    - HET: Mix of reference and alternate alleles (heterozygous)

    :param gt: Numpy array containing genotype allele indices where:
              0 = reference allele
              1+ = alternate allele
              -1 = missing allele
    :return: String indicating the zygosity classification
    """
    gt = list(gt)

    # All genotypes are missing
    if len(gt) == 0 or all(allele == -1 for allele in gt):
        return "MISSING"

    # One allele
    if len(gt) == 1:
        if gt[0] == 0:
            return "HOM_REF"
        return "HEMI"

    # More than one allele
    if all(allele == 0 for allele in gt):
        return "HOM_REF"
    if all(allele == gt[0] for allele in gt):
        return "HOM_ALT"
    return "HET"


def allele_string(ref: str, alt: str, gt: np.ndarray) -> str:
    """
    Convert a genotype array to a string representation of the allele combination.

    Creates a slash-delimited string showing the actual alleles present based on
    the reference allele, alternate allele, and genotype array. For example:
    - A/A for homozygous reference
    - A/T for heterozygous
    - T/T for homozygous alternate
    
    Only handles diploid (2 alleles) or haploid (1 allele) cases. For ploidy > 2,
    returns "NA".

    :param ref: Reference allele string (e.g. "A")
    :param alt: Alternate allele string (e.g. "T") 
    :param gt: Numpy array containing genotype allele indices where:
              0 = reference allele
              1 = alternate allele
    :return: String representation of the allele combination
    """

    def allele(i):
        if i == 0:
            return ref
        return alt

    gt = list(gt)

    # ploidy > 2 is not supported. We need more than one ref and one alt.
    if len(gt) > 2:
        return "NA"

    if len(gt) == 1:
        return allele(gt[0])

    if gt == [0, 0]:
        return f"{ref}/{ref}"
    if gt == [1, 1]:
        return f"{alt}/{alt}"

    return f"{ref}/{alt}"


@df_transform
def transform_vcf(
    vcf_df: pd.DataFrame,
    *,
    split_multiallelic: bool = True,
    add_zygosity: bool = True,
    add_allele_string: bool = True,
    drop_duplicates: bool = True,
    variant_loci: List[str] = None,
    verbose: bool = False,
    drop_samples: bool = False,
) -> pd.DataFrame:
    """
    Transform a VCF DataFrame by performing common preprocessing operations.

    This function applies several transformations to VCF data:
    1. Optionally drops duplicate variants
    2. Splits alleles into reference and alternate columns
    3. Handles multiallelic variants (split or join)
    4. Filters for specific variant loci if provided
    5. Adds zygosity classification
    6. Adds allele string representation
    7. Optionally drops sample information

    :param vcf_df: Input DataFrame containing VCF data
    :param split_multiallelic: If True, splits multiallelic variants into separate rows
    :param add_zygosity: If True, adds zygosity classification column
    :param add_allele_string: If True, adds allele string representation column
    :param drop_duplicates: If True, removes duplicate variant entries
    :param variant_loci: Optional list of specific variants to keep in "chrom-pos-ref-alt" format
    :param verbose: If True, enables debug logging
    :param drop_samples: If True, removes sample name column
    :return: Transformed DataFrame
    """

    level = logging.DEBUG if verbose else logging.INFO
    logger = get_logger(level)
    logger.propagate = False

    def log_event(message, t_prev=None):
        if t_prev:
            message += f" {time.time() - t_prev:.3f} sec"
        logger.debug(message)
        return time.time()

    t_prev = log_event("start transform")

    # Drop duplicates caused by multiple regions intersecting the same variant
    if drop_duplicates:
        vcf_df["str_alleles"] = vcf_df["alleles"].astype(str)
        vcf_df["str_gt"] = vcf_df["fmt_GT"].astype(str)
        vcf_df.drop_duplicates(
            subset=["contig", "pos_start", "sample_name", "str_alleles", "str_gt"],
            inplace=True,
        )
        vcf_df.drop(["str_alleles", "str_gt"], axis=1, inplace=True)
        t_prev = log_event("drop duplicates", t_prev)

    # Split alleles into ref and alt
    vcf_df["ref"] = vcf_df["alleles"].str[0]
    vcf_df["alt"] = vcf_df["alleles"].str[1:]
    vcf_df.drop("alleles", axis=1, inplace=True)
    t_prev = log_event("split ref/alt", t_prev)

    # Create an af column with the ALT IAF values
    alt_af = "info_TILEDB_ALT_IAF"
    if "info_TILEDB_IAF" in vcf_df:
        vcf_df[alt_af] = vcf_df["info_TILEDB_IAF"].str[1:]

    if split_multiallelic:
        # Split multiallelic variants
        explode_cols = ["alt"]

        if "af" in vcf_df:
            explode_cols.append(alt_af)

        vcf_df = vcf_df.explode(explode_cols)
        t_prev = log_event("split multiallelic", t_prev)
    else:
        # Convert alt to comma-separated string
        vcf_df["alt"] = vcf_df["alt"].str.join(",")
        t_prev = log_event("join multiallelic", t_prev)

    # Select records that match the requested chrom-pos-ref-alt
    if variant_loci:
        split_loci = [variant.split("-") for variant in variant_loci]

        loci_df = pd.DataFrame(
            {
                "contig": [x[0] for x in split_loci],
                "pos_start": [int(x[1]) for x in split_loci],
                "ref": [x[2] for x in split_loci],
                "alt": [x[3] for x in split_loci],
            }
        )

        vcf_df = pd.merge(
            loci_df, vcf_df, on=["contig", "pos_start", "ref", "alt"], how="inner"
        )
        t_prev = log_event("match chrom-pos-ref-alt records", t_prev)

        # Return if empty
        if len(vcf_df) == 0:
            return vcf_df

    # Add zygosity
    if add_zygosity:
        vcf_df["zygosity"] = vcf_df["fmt_GT"].apply(zygosity)
        t_prev = log_event("add zygosity", t_prev)

    # Add allele_string
    if add_allele_string:
        vcf_df["allele_string"] = vcf_df[["ref", "alt", "fmt_GT"]].apply(
            lambda row: allele_string(*row), axis=1
        )
        t_prev = log_event("add allele_string", t_prev)

    # Drop columns
    if drop_samples:
        vcf_df.drop(["sample_name"], axis=1, inplace=True)
    vcf_df.drop("fmt_GT", axis=1, inplace=True)
    return vcf_df

def aggregate_tables(tables: pa.Table) -> pa.Table:
    """Aggregates multiple PyArrow tables into a single table with counts.

    Args:
        tables: List of PyArrow tables to aggregate

    Returns:
        A PyArrow table containing the aggregated data with counts.
        Returns an empty table if no variants were found.
    """
    # Return if no variants were found
    tables = [x for x in tables if x is not None and x.num_rows > 0]
    if len(tables) == 0:
        return pa.table({})

    # Concatenate the tables
    table = pa.concat_tables(tables)

    df = table.to_pandas()

    # Aggregate and count
    df = (
        df.groupby(df.columns.tolist())
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )

    # Sort
    df.sort_values(["contig", "pos_start", "count"], inplace=True)

    return pa.Table.from_pandas(df)


def read_samples(vcf_uri: str, tiledb_cfg: Optional[dict] = None) -> List[str]:
    """Reads sample IDs from a TileDB-VCF dataset.

    Args:
        vcf_uri: URI of the TileDB-VCF dataset
        tiledb_cfg: Optional TileDB config dictionary

    Returns:
        List of sample IDs as strings
    """
    def fetch_samples(vcf_uri, tiledb_cfg):
        import pyarrow as pa
        import tiledbvcf

        sample_list = tiledbvcf.Dataset(vcf_uri, tiledb_config=tiledb_cfg).samples()
        # return as arrow
        return pa.array(sample_list)

    sample_dag = tiledb.cloud.dag.DAG(name="VCF-Federated-read-samples")
    samples_node = sample_dag.submit(fetch_samples, vcf_uri, tiledb_cfg)
    run_dag(sample_dag)
    return samples_node.result().to_pylist()


def genotype_query_by_variants(
    vcf_uri: str,
    *,
    variants: List[str],
    allele_string_mode: bool = False,
    zygosity_mode: bool = True,
    max_sample_batch_size: int = 1000,
    use_large_node: bool = False,
    verbose: bool = False,
    drop_samples: bool = True,
    drop_duplicates: bool = True,
    num_region_partitions: int = 1,
    threads: int = 2,
    tiledb_cfg: Optional[dict] = None,
    batch_mode: bool = False,

) -> pa.Table:
    """Queries genotypes for specific variants from a TileDB-VCF dataset.

    Args:
        vcf_uri: URI of the TileDB-VCF dataset
        variants: List of variant strings in format "chrom-pos-ref-alt"
        allele_string_mode: Whether to include allele strings in output
        zygosity_mode: Whether to include zygosity in output
        max_sample_batch_size: Maximum number of samples per batch
        use_large_node: Whether to use a large compute node
        verbose: Whether to output verbose logs
        drop_samples: Whether to drop sample IDs from output
        drop_duplicates: Whether to drop duplicate sample/variant combinations
        num_region_partitions: Number of region partitions for parallel processing
        threads: max_workers
        tiledb_cfg: cfg dict (mostly for x-payer)
        batch_mode: Whether to use batch task graphs instead of realtime

    Returns:
        PyArrow table containing the query results
    """
    import tiledb.cloud.vcf

    # Create a list of regions for each chromosome
    regions_by_chrom = defaultdict(list)
    for variant in variants:
        if not variant.startswith("#"):
            chrom, pos, _, _ = variant.split("-")
            regions_by_chrom[chrom].append(f"{chrom}:{pos}-{pos}")
        if not variant.startswith("#"):
            chrom, pos, _, _ = variant.split("-")
            regions_by_chrom[chrom].append(f"{chrom}:{pos}-{pos}")

    log = get_logger()
    log.propagate = False
    log.info(
        f"genotype_query_by_variants for {len(variants)} variants over {len(regions_by_chrom)} chrom"
    )

    # Read list of all samples
    samples = read_samples(vcf_uri, tiledb_cfg)

    def read_per_chrom(regions: List[str]) -> pa.Table:
        max_workers = math.ceil(len(samples) / max_sample_batch_size)

        # Submit the distributed query
        # Check if batch_mode is a valid parameter for tiledb.cloud.vcf.read
        import inspect
        read_params = inspect.signature(tiledb.cloud.vcf.read).parameters
        read_args = {
            "dataset_uri": vcf_uri,
            "regions": regions,
            "samples": samples,
            "transform_result": transform_vcf(
                split_multiallelic=True,
                add_zygosity=zygosity_mode,
                add_allele_string=allele_string_mode,
                variant_loci=variants,
                verbose=verbose,
                drop_samples=drop_samples,
            ),
            "max_sample_batch_size": max_sample_batch_size,
            "max_workers": max_workers,
            "num_region_partitions": num_region_partitions,
            "resource_class": "large" if use_large_node else None,
            "verbose": verbose,
            "config": tiledb_cfg,
        }
        
        # Check if batch_mode is a valid parameter for tiledb.cloud.vcf.read
        if "batch_mode" in read_params:
            read_args["batch_mode"] = batch_mode
            
        table = tiledb.cloud.vcf.read(**read_args)

        return table

    with ThreadPoolExecutor(max_workers=threads) as executor:
        tables = list(executor.map(read_per_chrom, regions_by_chrom.values()))
    agg_tables = aggregate_tables(tables)
    if drop_duplicates:
        # Keep only distinct sample/variant combinations
        agg_tables = pa.Table.from_pandas(agg_tables.to_pandas().drop_duplicates(subset=['sample_name', 'contig', 'pos_start', 'ref', 'alt']))
    return agg_tables

def genotype_query_by_bed_regions(
    vcf_uri: str,
    *,
    bed_regions: List[str],
    allele_string_mode: bool = False,
    max_sample_batch_size: int = 1000,
    use_large_node: bool = False,
    verbose: bool = False,
    drop_samples: bool = True,
    threads=2,
    tiledb_cfg: Optional[dict] = None,
) -> pa.Table:
    import tiledb.cloud.vcf

    # Create a list of regions for each chromosome
    regions_by_chrom = defaultdict(list)
    for region in bed_regions:
        # echo region being like chr1:0-999
        # BED is 0-based non-inclusive.
        # Covert to 1-based inclusive for VCF.
        # split on : or any whitespace
        # split on : or any whitespace
        chrom, start_pos, end_pos = (
            region.split(":")[0],
            int(region.split(":")[1].split("-")[0]) + 1,
            int(region.split(":")[1].split("-")[1]),
        )
        regions_by_chrom[chrom].append(f"{chrom}:{start_pos}-{end_pos}")

    # Read list of all samples
    samples = read_samples(vcf_uri, tiledb_cfg)

    def read_per_chrom(regions: List[str]) -> pa.Table:
        max_workers = math.ceil(len(samples) / max_sample_batch_size)

        # Submit the distributed query
        table = tiledb.cloud.vcf.read(
            vcf_uri,
            regions=regions,
            samples=samples,
            transform_result=transform_vcf(
                split_multiallelic=True,
                add_zygosity=not allele_string_mode,
                add_allele_string=allele_string_mode,
                verbose=verbose,
                drop_samples=drop_samples,
            ),
            max_sample_batch_size=max_sample_batch_size,
            max_workers=max_workers,
            num_region_partitions=1,
            resource_class="large" if use_large_node else None,
            verbose=verbose,
            config=tiledb_cfg,
        )

        return table

    # Read each chromosome in parallel
    # threads = len(regions_by_chrom.values())
    # threads = len(regions_by_chrom.values())
    with ThreadPoolExecutor(max_workers=threads) as executor:
        tables = list(executor.map(read_per_chrom, regions_by_chrom.values()))

    return aggregate_tables(tables)


def compute_compound_heterozygosity(vcf_df: pa.Table) -> pa.Table:
    """
    This takes an annotated vcf dataframe with a GENE column and zygosity column to determine cmpd hets
    """
    log = get_logger()
    vcf_df = vcf_df.to_pandas()
    hets = vcf_df[vcf_df["zygosity"] == "HET"]
    
    homs = vcf_df[vcf_df["zygosity"] == "HOM_ALT"]
    #just unique sample_name,GENE,zygosity
    homs = homs[["sample_name", "GENE"]].drop_duplicates()
    homs['passenger_het'] = True

    if len(hets) > 0:
        grouped_hets = (
            hets.groupby(["sample_name", "GENE"]).size().reset_index(name="het_count")
        )
        compound_hets = grouped_hets[grouped_hets["het_count"] > 1].copy()
        if len(compound_hets) > 0:
            # avoid the A value is trying to be set on a copy of a slice from a DataFrame error
            compound_hets.loc[:, "compound_event"] = "CMPD_HET"
            compound_hets.loc[:, "zygosity"] = "HET"  # for the join
            log.info(f"{len(compound_hets)} compound hets found")
            vcf_df = vcf_df.merge(
                compound_hets, on=["sample_name", "zygosity", "GENE"], how="left"
            )
        else:
            # columns should match up
            vcf_df.loc[:, "het_count"] = 0
            vcf_df.loc[:, "compound_event"] = ""
    else:
        vcf_df.loc[:, "het_count"] = 0
        vcf_df.loc[:, "compound_event"] = ""

    if len(homs) > 0:
    #set all vcf_df hets to compound_event HET_PASS if they are in homs sample, GENE
        log.info(f"found {len(homs)} homs")
        #import pdb; pdb.set_trace()
        vcf_df = vcf_df.merge(homs, on=["sample_name", "GENE"], how="left")
        vcf_df["passenger_het"] = vcf_df["passenger_het"].fillna(False)
        #if passenger_het is True, set compound_event to HET_PASS just using pandas not numpy
        vcf_df.loc[(vcf_df["passenger_het"] == True) & (vcf_df["zygosity"] == "HET"), "compound_event"] = "HET_PASS"
        #remove passenger_het
        #vcf_df = vcf_df.drop(columns=["passenger_het"])
    else:
        vcf_df.loc[:, "passenger_het"] = False
    #if compound_event is HET_PASS, set zygosity to HET        

    vcf_df["het_count"] = vcf_df["het_count"].fillna(0)
    vcf_df["het_count"] = vcf_df["het_count"].astype("int64")
    return pa.Table.from_pandas(vcf_df.reset_index())


def compute_gene_diplotype(vcf_df: pa.Table, gene: str) -> pa.Table:
    """
    This takes an annotated vcf dataframe, with a GENE column and a GENE to determine all diplotypes that exist in that gene
    """
    log = get_logger()
    vcf_df = vcf_df.to_pandas()
    vcf_df = vcf_df[vcf_df['GENE'] == gene]

    grouped_by_sample = (
            vcf_df.groupby(["sample_name"])
    )
    #for each sample display all the allele strings
    diplotypes = []
    for sample, sample_df in grouped_by_sample:
        allele_strings = sample_df["allele_string"].tolist()
        diplotypes.append({
            "sample_name": sample,
            "diplotype": "/".join(sorted(allele_strings))
        })
    
    diplotype_df = pd.DataFrame(diplotypes)
    vcf_df = vcf_df.merge(diplotype_df, on="sample_name", how="left")
    hets = vcf_df[vcf_df["zygosity"] == "HET"]
    if len(hets) > 0:
        grouped_hets = (
            hets.groupby(["sample_name", "GENE"]).size().reset_index(name="het_count")
        )
        compound_hets = grouped_hets[grouped_hets["het_count"] > 1].copy()
        if len(compound_hets) > 0:
            # avoid the A value is trying to be set on a copy of a slice from a DataFrame error
            compound_hets.loc[:, "compound_event"] = "CMPD_HET"
            compound_hets.loc[:, "zygosity"] = "HET"  # for the join
            log.info(f"{len(compound_hets)} compound hets found")
            vcf_df = vcf_df.merge(
                compound_hets, on=["sample_name", "zygosity", "GENE"], how="left"
            )
        else:
            # columns should match up
            vcf_df.loc[:, "het_count"] = 0
            vcf_df.loc[:, "compound_event"] = ""
    else:
        vcf_df.loc[:, "het_count"] = 0
        vcf_df.loc[:, "compound_event"] = ""

    vcf_df["het_count"] = vcf_df["het_count"].fillna(0)
    vcf_df["het_count"] = vcf_df["het_count"].astype("int64")
    return pa.Table.from_pandas(vcf_df.reset_index())


# Endpoints fot the dispatch function


def gc_variants(
    vcf_uri: str,
    *,
    variants_uri: str,
    variants: List[str],
    allele_string_mode: bool = False,
    max_sample_batch_size: int = 1000,
    use_large_node: bool = False,
    tiledb_cfg=None,
    verbose: bool = False,
) -> pa.Table:

    # Read variants
    if variants_uri is not None:
        variants = [
            variant
            for variant in read_lines(variants_uri)
            if not variant.startswith("#")
        ]

    return genotype_query_by_variants(
        vcf_uri,
        variants=variants,
        allele_string_mode=allele_string_mode,
        max_sample_batch_size=max_sample_batch_size,
        use_large_node=use_large_node,
        verbose=verbose,
        tiledb_cfg=tiledb_cfg,
    )


def gc_bed(
    vcf_uri: str,
    *,
    bed_lines: List[str],
    allele_string_mode: bool = False,
    max_sample_batch_size: int = 1000,
    use_large_node: bool = False,
    tiledb_cfg: Optional[dict] = None,
    verbose: bool = False,
) -> pa.Table:
    # Convert BED lines to BED regions
    bed_regions = [
        f"{line.split()[0]}:{line.split()[1]}-{line.split()[2]}" for line in bed_lines
    ]

    return genotype_query_by_bed_regions(
        vcf_uri,
        bed_regions=bed_regions,
        allele_string_mode=allele_string_mode,
        max_sample_batch_size=max_sample_batch_size,
        use_large_node=use_large_node,
        tiledb_cfg=tiledb_cfg,
        verbose=verbose,
    )


from vcf_federated.allele_count_dispatch import query_ensembl_by_gene_name
from vcf_federated.allele_count_dispatch import to_regions


def gc_gene(
    vcf_uri: str,
    *,
    ens_uri: str,
    gene_name: str,
    allele_string_mode: bool = False,
    max_sample_batch_size: int = 1000,
    use_large_node: bool = False,
    tiledb_cfg: Optional[dict] = None,
    verbose: bool = False,
):
    # Get BED regions for the gene
    kwargs = {"ensembl_uri": ens_uri}
    table = query_ensembl_by_gene_name(
        gene_name=gene_name, tiledb_cfg=tiledb_cfg, **kwargs
    )
    bed_regions = to_regions(table, bed=True)

    return genotype_query_by_bed_regions(
        vcf_uri,
        bed_regions=bed_regions,
        allele_string_mode=allele_string_mode,
        max_sample_batch_size=max_sample_batch_size,
        use_large_node=use_large_node,
        tiledb_cfg=tiledb_cfg,
        verbose=verbose,
    )


def fetch_approved_sample_metadata(
    sample_metadata_uri, sample_use_uri=None, tiledb_cfg=None, population_source: str = "rady", metadata_attrs: dict = None
):
    """
    Fetch sample metadata with proper sample use
    """
    log = get_logger()

    def query_sample_metadata(
        sample_metadata_uri: str,
        sample_stems: list,
        family_roles: list,
        genome_id_list: pa.lib.StringArray,
        population_source: str = "rady",
        metadata_attrs: dict = None,
    ) -> pa.Table:
        """
        This function queries the metadata to lookup sample names from genome ids
        Arguments:
            @param genome_id_list: the genome_ids of interest
        Returns:
            Sample data subset
        """
        
        sample_metadata_attrs = metadata_attrs["sample_metadata_attrs"][population_source]
        with tiledb.scope_ctx(tiledb_cfg):
            with tiledb.open(sample_metadata_uri) as A:
                df = A.query(attrs=sample_metadata_attrs).df[:]
                if sample_stems is not None and len(sample_stems) > 0:
                    df = df[df.sample_name.str.contains("|".join(sample_stems))]
                if family_roles is not None and len(family_roles) > 0:
                    df = df[df.family_relationship.isin(family_roles)]
                if genome_id_list is not None and len(genome_id_list) > 0:
                    genome_ids = genome_id_list.to_pylist()
                    df = df[df.genome_id.isin(genome_ids)]
                # remove -1 from end of sample names, just for rady
                if population_source == 'rady':
                    df["sample_name"] = df["sample_name"].str.replace(r"-1$", "", regex=True)
            return pa.Table.from_pandas(df)

    sample_graph = tiledb.cloud.dag.DAG(
        max_workers=1, name="FEDERATED QUERY: fetch samples nbs "
    )
    filtered_metadata_node = sample_graph.submit(
        query_sample_metadata,
        sample_metadata_uri=sample_metadata_uri,
        sample_stems=None,
        family_roles=None,
        genome_id_list=None,
        population_source=population_source,
        name=f"FEDERATED QUERY: Sample Metadata Filter",
        resource_class="standard",
        metadata_attrs=metadata_attrs,
    )
    sample_graph.compute()
    sample_graph.wait()
    filtered_metadata_pre_use = filtered_metadata_node.result().to_pandas()
    if sample_use_uri is not None:
        with tiledb.scope_ctx(tiledb_cfg):
            with tiledb.open(sample_use_uri, mode="r") as su:
                sample_use = su.query().df[:]
                filtered_metadata_use = filtered_metadata_pre_use.merge(
                    sample_use,
                    how="inner",
                    left_on=["individual_id", "family_relationship", "case_id"],
                    right_on=["ind_id", "family_relationship", "case_id"],
                )
                filtered_metadata = filtered_metadata_use[
                    filtered_metadata_use["additional_use"] == "Yes"
                ]
    else:
        filtered_metadata = filtered_metadata_pre_use
    return filtered_metadata


def apply_blocklist(variant_df, blocklist_uri, tiledb_cfg):
    log = get_logger()
    with tiledb.scope_ctx(tiledb_cfg):
        with tiledb.open(blocklist_uri, mode="r") as bl:
            blocklist = bl.query().df[:]
        # rename blocklist columns to CHROM, POS, REF, ALT
        blocklist = blocklist.rename(
            columns={"pos_start": "POS", "ref": "REF", "alt": "ALT"}
        )
        # drop the columns we don't need
        blocklist = blocklist[["CHROM", "POS", "REF", "ALT"]]
        blocked_entries = variant_df.merge(
            blocklist, how="inner", on=["CHROM", "POS", "REF", "ALT"]
        )
        blocked_entries_subset = blocked_entries[
            ["CHROM", "POS", "REF", "ALT", "GENE", "CLINVAR_ID"]
        ].drop_duplicates()
        block_all = variant_df.merge(
            blocked_entries_subset,
            on=["CHROM", "POS", "REF", "ALT", "GENE", "CLINVAR_ID"],
            how="left",
            indicator=True,
        )
        # remove the entries that are in the blocked_entries
        variant_df = block_all[block_all["_merge"] == "left_only"]
        return variant_df


def positive_genotype(row, sex_field, sex_lookup, accept_passenger_hets=False):
    log = get_logger()
    if row["moi"] == "Pattern unknown":
        return "Unknown"
    elif row["moi"] == "AD":
        return "Yes"
    elif row["moi"] == "AR":
        if row["zygosity"] == "HOM_ALT" or row["compound_event"] == "CMPD_HET" or (accept_passenger_hets and row["compound_event"] == "HET_PASS"):
            return "Yes"
        else:
            return "No"
    elif row["moi"] == "XR":
        if sex_field not in row:
            return "No"
        else:
            if row[sex_field] == sex_lookup["Female"]:
                if row["zygosity"] == "HOM_ALT" or row["compound_event"] == "CMPD_HET" or (accept_passenger_hets and row["compound_event"] == "HET_PASS"):
                    return "Yes"
            else:
                if row["zygosity"] == "HEMI":
                    return "Yes"
            return "No"
    elif row["moi"] == "XD":
        if sex_field not in row:
            return "Yes"
        else:
            if row[sex_field] == sex_lookup["Male"]:
                if row["zygosity"] == "HEMI":
                    return "Yes"
            else:
                return "Yes"


def beginngs_query(
    vcf_uri: str,
    moaavt_uri: str,
    moi_uri: str,
    blocklist_uri: str,
    bed_lines: List[str],
    sample_metadata_uri: str,
    sample_use_uri: str,
    use_blocklist: bool = True,
    tiledb_cfg=None,
    batch_mode: bool = False,
    max_sample_batch_size: int = 2000,
    threads: int = 4,
    aggregate_variants=True,
    verbose: bool = False,
    population_source: str = "rady",
    format_style: str = 'alexion',
    metadata_attrs: dict = None,
    positive_only: bool = True,
):
    """
    Query variants and genotypes from a mother-of-all-annotated variants table (MOAAVT).

    This function reads a MOAAVT containing annotated variants, queries the variants and their genotypes,
    and returns aggregated genotype and diplotype counts. It can optionally filter variants using a
    blocklist and/or BED regions.

    Args:
        vcf_uri: URI of the VCF array
        moaavt_uri: URI of the MOAAVT array containing annotated variants
        moi_uri: URI of the mode of inheritance array
        blocklist_uri: URI of the variant blocklist array
        bed_lines: Optional list of BED regions to restrict variant queries
        sample_metadata_uri: URI of the sample metadata array
        sample_use_uri: URI of the sample usage permissions array
        use_blocklist: Whether to filter variants using the blocklist (default: True)
        tiledb_cfg: Optional TileDB config dictionary
        batch_mode: Whether to process samples in batches (default: False)
        max_sample_batch_size: Maximum samples per batch when batch_mode=True (default: 2000)
        threads: Number of threads to use for processing (default: 4)
        aggregate_variants: Whether to aggregate variant counts (default: True)
        verbose: Enable verbose logging (default: False)
        population_source: Source population for metadata lookup (default: "rady")

    Returns:
        BeginNGSResult containing genotype and diplotype counts as pyarrow Tables
    """
    log = get_logger()
    log.propagate = False
    import pathlib

    import json
    import os
    import pkg_resources
    import yaml
    import sys

    if metadata_attrs is None:
        # this won't work within the UDF context, it must be passed from the package
        metadata_attrs_path = pkg_resources.resource_filename('vcf_federated', 'config/metadata_attrs.json')
        with open(metadata_attrs_path) as f:
            metadata_attrs = json.load(f)

    def load_variant_data(moaavt_uri, tiledb_cfg):
        with tiledb.scope_ctx(tiledb_cfg):
            with tiledb.open(moaavt_uri, mode="r") as mu:
                variant_df = mu.df[:]
        variant_df = variant_df.reset_index()
        variant_df = variant_df.rename(columns={"CLINVAR_VARIANT_ID": "CLINVAR_ID"})
        return variant_df

    def filter_variants_by_bed(variant_df, bed_lines):
        if not bed_lines:
            log.info("no bed regions found")
            return variant_df
            
        bed_regions = [
            f"{line.split()[0]}:{line.split()[1]}-{line.split()[2]}"
            for line in bed_lines
        ]
        log.info(f"there are {len(bed_regions)} bed regions")
        
        variants_in_bed_region = []
        for region in bed_regions:
            chrom, start_pos, end_pos = (
                region.split(":")[0],
                int(region.split(":")[1].split("-")[0]) + 1,
                int(region.split(":")[1].split("-")[1]),
            )
            variants_in_bed_region += [
                variant_df[
                    (variant_df["CHROM"] == chrom)
                    & (variant_df["POS"] >= start_pos)
                    & (variant_df["POS"] <= end_pos)
                ]
            ]
        return pd.concat(variants_in_bed_region)

    def get_genotype_data(vcf_uri, variant_df, batch_mode, max_sample_batch_size, threads, verbose, tiledb_cfg):
        variant_loci = variant_df.apply(lambda row: f"{row['CHROM']}-{row['POS']}-{row['REF']}-{row['ALT']}", axis=1).tolist()
        
        log.info(f"Querying genotype_query_by_variants for {len(variant_loci)} variants {variant_loci[:5]}...")
        log.info(f"genotype_query_by_variants")
        log.info(f"batch_mode: {batch_mode}")
        t_start = time.time()

        table = genotype_query_by_variants(
            vcf_uri,
            variants=variant_loci,
            allele_string_mode=True,
            zygosity_mode=True,
            batch_mode=batch_mode,
            max_sample_batch_size=max_sample_batch_size,
            threads=threads,
            use_large_node=False,
            verbose=verbose,
            drop_samples=False,
            num_region_partitions=1,
            tiledb_cfg=tiledb_cfg,
        ).to_pandas()

        log.info(f"genotype_query_by_variants finished")
        log.info(f"Done in {time.time() - t_start:.3f} seconds")

        return table

    def process_compound_heterozygosity(table):
        genes = table["GENE"].unique()
        cmpd_tables = []
        for gene in genes:
            gene_table = table[table["GENE"] == gene]
            log.info(f"{len(gene_table)} variants found for gene {gene}")
            cmpd_tables.append(
                compute_compound_heterozygosity(pa.Table.from_pandas(gene_table))
            )
        return pa.concat_tables(cmpd_tables)

    def merge_sample_metadata(vcf_df, sample_metadata_uri, sample_use_uri, population_source, metadata_attrs):
        if sample_metadata_uri is None:
            return vcf_df
            
        filtered_metadata = fetch_approved_sample_metadata(
            sample_metadata_uri, sample_use_uri, population_source=population_source, metadata_attrs=metadata_attrs
        )
        
        if population_source == "rady":
            log.info("removing sample name suffix for merge")
            vcf_df["sample_name"] = vcf_df["sample_name"].str.replace(r"-1$", "", regex=True)
            filtered_metadata["sample_name"] = filtered_metadata["sample_name"].str.replace(r"-1$", "", regex=True)
        elif population_source == "ukbb-alexion":
            log.info("renaming sample column to sample_name for ukbb-alexion")
            filtered_metadata = filtered_metadata.rename(columns={"cgr_sequence_id": "sample_name"})
            
        filtered_metadata["sample_name"] = filtered_metadata["sample_name"].astype(str)
        
        sm_samples = list(filtered_metadata["sample_name"])
        log.info(f"{len(sm_samples)} approved samples found in sample metadata")
        log.info(f"{len(vcf_df[['sample_name','contig','pos_start','ref','alt']].drop_duplicates())} tuples before sample filter")
        log.info(f"merging on sample_name {vcf_df['sample_name'][:5]} vs {filtered_metadata['sample_name'][:5]} (ensure strings)")
        
        vcf_df = pd.merge(vcf_df, filtered_metadata, on=["sample_name"], how="inner")
        log.info(f"{len(vcf_df[['sample_name','contig','pos_start','ref','alt']].drop_duplicates())} tuples after sample filter")
        
        return vcf_df

    def aggregate_variants_data(vcf_df,positive_only=True):
        log.info("Group by variant")
        log.info("Set zygosity for compound_event")
        vcf_df.loc[(vcf_df['compound_event'] == 'CMPD_HET') | (vcf_df['compound_event'] == 'HET_PASS'), 'zygosity'] = 'CMPD_HET'

        if vcf_df.duplicated(subset=['sample_name', 'chrposrefalt']).any():
            duplicates = vcf_df[vcf_df.duplicated(subset=['sample_name', 'chrposrefalt'], keep=False)]
            log.info(f"Found {len(duplicates)} duplicate rows:")
            log.info(duplicates)
            raise ValueError("Found duplicate sample_name/variant combinations - cannot pivot. Please deduplicate first.")

        if positive_only:
            log.info(f"There are {len(vcf_df)} variants before positive only filter")
            vcf_df = vcf_df[vcf_df['positive_genotype'] == 'Yes']
            log.info(f"There are {len(vcf_df)} variants after positive only filter")

        unique_genes = sorted(vcf_df['GENE'].unique())
        log.info(f"Found {len(unique_genes)} unique genes")

        all_gene_pivot_counts = []
        for gene in unique_genes:
            #if gene == 'GALT':
            #    import pdb; pdb.set_trace()
            gene_df = vcf_df[vcf_df['GENE'] == gene]
            #make sure this is ordered by pos_start
            gene_df = gene_df.sort_values('pos_start')
            #get the MOI for this gene
            moi = vcf_df[vcf_df['GENE'] == gene]['moi'].unique()[0]
            #assert this MOI is the same for all variants in this gene
            assert len(vcf_df[vcf_df['GENE'] == gene]['moi'].unique()) == 1
            log.info(f"Processing gene {gene} with {len(gene_df)} variants")
            
            gene_pivot = gene_df.pivot(index='sample_name',
                                     columns='chrposrefalt',
                                     values='zygosity')
            gene_pivot = gene_pivot.fillna('')
            gene_pivot_orig = gene_pivot.copy()
            def sort_by_zygosity(items):
                # Define zygosity order (highest priority first)
                zygosity_order = {
                    'HOM_ALT': 0,
                    'CMPD_HET': 1,
                    'HEMI': 2,
                    'HET': 3
                }
                # Sort items by zygosity order, then by position
                return sorted(items, key=lambda x: (zygosity_order.get(x[1], 3), x[0]))

            gene_pivot['DIPLOTYPE_RAW'] = gene_pivot_orig.apply(
                lambda row: ';'.join(
                    f"[{pos}]" if (zyg == 'CMPD_HET' or zyg == 'HEMI') 
                    else f"[{pos}];[{pos}]" if zyg == 'HOM_ALT' 
                    else f"[{pos}];+" 
                    for pos, zyg in sort_by_zygosity(row.items()) if zyg != ''
                ), 
                axis=1
            )

            

            def combine_diplotypes(diplotype):
                # Arrange diplotypes to segregate homozygous variants in separate brackets
                # Keep heterozygous variants in the second bracket if they coexist with homozygous variants
                # Compound het variants that do not coexist with homozygous variants should be segregated with
                # the lowest position variant occupying the first bracket and the rest in the second bracket

                #CFTR
                #gene_pivot.loc['HG00113']['DIPLOTYPE'] [chr7-117548758-G-T];[chr7-117590400-G-C;chr7-117592169-C-T]
                #gene_pivot.loc['HG01075']['DIPLOTYPE'] '[chr6-32039081-C-A;chr6-32040110-G-T];[chr6-32039081-C-A;chr6-32040110-G-T]'

                #CYP21A2
                #gene_pivot.loc['HG01708']['DIPLOTYPE'] '[chr6-32039081-C-A];[chr6-32039081-C-A;chr6-32038514-C-T]'

                # Split into individual variants
                variants = [v.strip('[]') for v in diplotype.split(';') if v]
                
                # If we don't have more than 2 variants, return original
                if len(variants) <= 2:
                    return diplotype
                
                # Count occurrences of each variant to determine homozygosity
                variant_counts = {}
                for variant in variants:
                    variant_counts[variant] = variant_counts.get(variant, 0) + 1
                
                # Get unique variants while preserving order
                unique_variants = list(dict.fromkeys(variants))
                
                # Separate homozygous and heterozygous variants
                hom_variants = [v for v in unique_variants if variant_counts[v] > 1]
                het_variants = [v for v in unique_variants if variant_counts[v] == 1]
                
                # Build diplotype string
                result = []
                
                if len(unique_variants) >= 1:
                    # First bracket: homozygous variants first, then first heterozygous variant if available
                    first_bracket = []
                    if hom_variants:
                        first_bracket.extend(hom_variants)
                    else:
                        first_bracket.append(het_variants[0])
                    result.append(f"[{';'.join(first_bracket)}]")
                    #[chr3-122284869-C-T];+	HETEROZYGOUS	CASR	AD	2
                    if not hom_variants and het_variants:  # Only heterozygous variants and at least one exists
                        first_bracket.append("+")
                    
                    # Second bracket: homozygous variants first, then remaining heterozygous variants
                    second_bracket = []
                    second_bracket.extend(hom_variants)  # Include homozygous variants again
                    if hom_variants and het_variants:
                        second_bracket.append(het_variants[0])
                    elif len(het_variants) > 1:
                        # compound het, one in each bracket
                        second_bracket.extend(het_variants[1:])
                    if second_bracket:
                        #sort by pos in chr9-34649445-A-G ascending
                        second_bracket.sort(key=lambda x: int(x.split('-')[1]))
                        result.append(f"[{';'.join(second_bracket)}]")
                
                return ';'.join(result)
                
            gene_pivot['DIPLOTYPE'] = gene_pivot['DIPLOTYPE_RAW'].apply(combine_diplotypes)

            gene_pivot['ZYGOSITY'] = gene_pivot_orig.apply(lambda row: 'HOMOZYGOUS' if 'HOM_ALT' in row.values else 'COMPOUND HETEROZYGOUS' if 'CMPD_HET' in row.values else  'HEMIZYGOUS' if 'HEMI' in row.values else 'HETEROZYGOUS' if any(v != '' for v in row.values) else '', axis=1)

             # For hemizygous cases, combine all variants into a single bracket
            def fix_hemizygous_diplotypes(row):
                if row['ZYGOSITY'] == 'HEMIZYGOUS':
                    # Extract all variants from the diplotype
                    variants = []
                    for bracket in row['DIPLOTYPE'].split(';'):
                        if bracket.startswith('[') and bracket.endswith(']'):
                            variants.extend(bracket[1:-1].split(';'))
                    # Sort variants by position
                    variants.sort(key=lambda x: int(x.split('-')[1]))
                    # Return single bracket with all variants
                    return f"[{';'.join(variants)}]"
                return row['DIPLOTYPE']
            

            gene_pivot['DIPLOTYPE'] = gene_pivot.apply(fix_hemizygous_diplotypes, axis=1)

            include_samples = False
            if include_samples:
                # Group by DIPLOTYPE and ZYGOSITY first
                grouped = gene_pivot.groupby(['DIPLOTYPE', 'ZYGOSITY'])
                # For each group, collect the sample names (which are in the index)
                gene_pivot_counts = grouped.apply(lambda x: ', '.join(x.index)).reset_index(name='samples')
            else:
                gene_pivot_counts = gene_pivot.groupby(['DIPLOTYPE','ZYGOSITY']).size().reset_index(name='count')
                gene_pivot_counts = gene_pivot_counts.sort_values('count', ascending=False)
            gene_pivot_counts['GENE'] = gene
            gene_pivot_counts['MOI'] = moi
            
            all_gene_pivot_counts.append(gene_pivot_counts)

        return pd.concat(all_gene_pivot_counts)

    def get_zygosity_counts(vcf_df):
        vcf_df_zygosity = (
            vcf_df.groupby(
                ["contig", "pos_start", "ref", "alt", "GENE", "CLINVAR_ID", "zygosity"]
            )["sample_name"]
            .nunique()
            .reset_index(name="sample_count")
        )

        vcf_df_positivity = (
            vcf_df.groupby(
                ["contig", "pos_start", "ref", "alt", "positive_genotype"]
            )["sample_name"]
            .nunique()
            .reset_index(name="sample_count")
        )

        log.info("pivoting by zygosity")
        vcf_df_zygosity_pivoted = vcf_df_zygosity.pivot_table(
            index=["contig", "pos_start", "ref", "alt", "GENE", "CLINVAR_ID"],
            columns="zygosity",
            values="sample_count",
            fill_value=0
        ).reset_index()

        log.info("pivoting by positivity") 
        vcf_df_positivity_pivoted = vcf_df_positivity.pivot_table(
            index=["contig", "pos_start", "ref", "alt"],
            columns="positive_genotype",
            values="sample_count",
            fill_value=0
        ).reset_index()
        vcf_df_positivity_pivoted.rename(columns={"Yes":"Positive"},inplace=True)
        vcf_df_positivity_pivoted.drop(columns=['No'],inplace=True)

        required_columns = ["contig", "pos_start", "ref", "alt"]
        log.info("merge zygosity and positivity")
        return pd.merge(
            vcf_df_zygosity_pivoted,
            vcf_df_positivity_pivoted,
            on=required_columns,
            how="left"
        )

    if verbose:
        log.info(f"beginngs query using vcf_uri: {vcf_uri} and moaavt:{moaavt_uri}")

    variant_df = load_variant_data(moaavt_uri, tiledb_cfg)
    log.info(f"{len(variant_df)} variants found in the moaavt")

    if use_blocklist and blocklist_uri is not None:
        variant_df = apply_blocklist(variant_df=variant_df, blocklist_uri=blocklist_uri, tiledb_cfg=tiledb_cfg)

    if bed_lines:
        variant_df = filter_variants_by_bed(variant_df, bed_lines)

    table = get_genotype_data(vcf_uri, variant_df, batch_mode, max_sample_batch_size, threads, verbose, tiledb_cfg)

    
    if len(table) == 0:
        log.info("No variants of interest found in the TileDB-VCF dataset")
        return pa.Table.from_pandas(table)

    table = table.drop("count", axis=1)
    variant_df_core = variant_df[["CHROM", "POS", "REF", "ALT", "GENE", "CLINVAR_ID"]]
    variant_df_core = variant_df_core.rename(
        columns={"CHROM": "contig", "POS": "pos_start", "REF": "ref", "ALT": "alt"}
    )
    table = table.merge(variant_df_core, on=["contig", "pos_start", "ref", "alt"], how="inner")

    table = process_compound_heterozygosity(table)
    vcf_df = table.to_pandas()
    vcf_df = merge_sample_metadata(vcf_df, sample_metadata_uri, sample_use_uri, population_source, metadata_attrs)


    if population_source not in metadata_attrs["sample_metadata_attrs"]:
        raise ValueError(f"Unknown sample metadata population source: {population_source}")
    if population_source not in metadata_attrs["sex_fields"]:
        raise ValueError(f"Unknown sex field population source: {population_source}")
    
    sample_metadata_attrs = metadata_attrs["sample_metadata_attrs"][population_source]
    sex_field = metadata_attrs["sex_fields"][population_source]["field"]
    sex_lookup = metadata_attrs["sex_fields"][population_source]["lookup"]

    if moi_uri is not None:
        with tiledb.open(moi_uri,config=tiledb_cfg) as moi_fh:
            moi_df = moi_fh.df[:]

        gene_poi_map = {gene.upper(): poi for gene, poi in zip(moi_df["GENE"], moi_df["MOI"])}
        vcf_df["moi"] = vcf_df["GENE"].map(gene_poi_map)
        vcf_df["positive_genotype"] = vcf_df.apply(positive_genotype, sex_field=sex_field, sex_lookup=sex_lookup, accept_passenger_hets=True, axis=1)
    else:
        vcf_df["positive_genotype"] = "N/A"
    
    vcf_df = vcf_df.drop_duplicates(subset=['sample_name', 'contig', 'pos_start', 'ref', 'alt'])
    if format_style == 'alexion':
        vcf_df['chrposrefalt'] = vcf_df.apply(lambda row: f"{row['contig']}-{row['pos_start']}-{row['ref']}-{row['alt']}", axis=1)
    else:
        vcf_df['chrposrefalt'] = vcf_df.apply(lambda row: f"{row['contig']}:{row['pos_start']}:{row['ref']}:{row['alt']}", axis=1)

    import pdb; pdb.set_trace()
    if aggregate_variants:
        pivot_df_counts = aggregate_variants_data(vcf_df,positive_only=positive_only)
        vcf_df = get_zygosity_counts(vcf_df)
    else:
        return vcf_df



    genotypes = pa.Table.from_pandas(vcf_df)
    diplotypes = pa.Table.from_pandas(pivot_df_counts)
    
    @dataclass
    class BeginNGSResult:
        genotypes: pa.Table
        diplotypes: pa.Table

    return BeginNGSResult(genotypes=genotypes, diplotypes=diplotypes)
    
def test_beginngs():


    
    df = beginngs_query(
        vcf_uri="tiledb://XXXXXXX",
        moaavt_uri="tiledb://XXXXXXX",
        moi_uri="tiledb://XXXXXXX",
        blocklist_uri="tiledb://XXXXXXX",
        sample_metadata_uri="tiledb://XXXXXXX",
        sample_use_uri="tiledb://XXXXXXX",
        use_blocklist=True,
        bed_lines=None,
        tiledb_cfg=None,
        batch_mode=False,
        max_sample_batch_size=2000,
        threads=4,
        namespace=None,
        remove_sample_name_suffix=True,
        verbose=False,
        population_source="rady"
    ).to_pandas()
    print(df)

if __name__ == "__main__":
    test_beginngs()


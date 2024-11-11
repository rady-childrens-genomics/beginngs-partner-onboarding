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
    Convert genotype to a zygosity string.

    :param gt: genotype
    :return: zygosity string
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
    Convert genotype to a allele_string string.

    :param ref: reference allele
    :param alt: alternate allele
    :param gt: genotype
    :return: allele_string string
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

) -> pa.Table:
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
        table = tiledb.cloud.vcf.read(
            vcf_uri,
            regions=regions,
            samples=samples,
            transform_result=transform_vcf(
                split_multiallelic=True,
                add_zygosity=zygosity_mode,
                add_allele_string=allele_string_mode,
                variant_loci=variants,
                verbose=verbose,
                drop_samples=drop_samples,
            ),
            max_sample_batch_size=max_sample_batch_size,
            max_workers=max_workers,
            num_region_partitions=num_region_partitions,
            resource_class="large" if use_large_node else None,
            verbose=verbose,
            config=tiledb_cfg,
        )

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
    sample_metadata_uri, sample_use_uri=None, tiledb_cfg=None, population_source: str = "rady"
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
    ) -> pa.Table:
        """
        This function queries the metadata to lookup sample names from genome ids
        Arguments:
            @param genome_id_list: the genome_ids of interest
        Returns:
            Sample data subset
        """
        metadata_attrs = {}
        metadata_attrs['rady'] = [
            "site_id",
            "study_id",
            "case_id",
            "individual_id",
            "biospecimen",
            "aliquot_id",
            "test_type",
            "family_relationship",
            "version",
            "sample_name",
            "gender",
            "fabric_report_id",
        ]
        metadata_attrs['ukbb-alexion'] = [
            "sample_name",
            "sex_encoded",
            "genetic_sex_encoded"
        ]
        metadata_attrs['onekg'] = [
            "case_id",
            "sample_name",
            "pid",
            "mid",
            "sex",
            "sexf",
            "pop",
            "reg",
            "population",
            "region",
            "phase3",
            "trio",
            "family_relationship",
            "individual_id"
        ]
        sample_metadata_attrs = metadata_attrs[population_source]
        with tiledb.scope_ctx(tiledb_cfg):
            with tiledb.open(sample_metadata_uri) as A:
                df = A.query(attrs=sample_metadata_attrs)[:]
                df = pd.DataFrame(df)
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


def positive_genotype(row, sex_field, sex_lookup, population_source):
    log = get_logger()
    if row["moi"] == "Pattern unknown":
        return "Unknown"
    elif row["moi"] == "AD":
        return "Yes"
    elif row["moi"] == "AR":
        if row["zygosity"] == "HOM_ALT" or row["compound_event"] == "CMPD_HET":
            return "Yes"
        else:
            return "No"
    elif row["moi"] == "XR":
        if sex_field not in row:
            return "No"
        else:
            if row[sex_field] == sex_lookup["Female"]:
                if row["zygosity"] == "HOM_ALT" or row["compound_event"] == "CMPD_HET":
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
    namespace=None,
    aggregate_variants=True,
    remove_sample_name_suffix=True,
    verbose: bool = False,
    population_source: str = "rady",
):
    """
    This reads in an Alexion-style mother-of-all-annotated variants table (moaavt), queries variants and return them annotated
    """
    log = get_logger()
    log.propagate = False
    import tomli
    import pathlib

    if verbose:
        log.info(f"beginngs query using vcf_uri: {vcf_uri} and moaavt:{moaavt_uri}")
    #TODO: we need to open the mooavt with the UDF
    with tiledb.scope_ctx(tiledb_cfg):
        with tiledb.open(moaavt_uri, mode="r") as mu:
            variant_df = mu.df[:]
        # reset this index
        # TODO: actually use an index in some intelligent manner
        variant_df = variant_df.reset_index()
        log.info(f"{len(variant_df)} variants found in the moaavt")



        variant_df = variant_df.rename(columns={"CLINVAR_VARIANT_ID": "CLINVAR_ID"})

        if use_blocklist == True and blocklist_uri is not None:
            variant_df = apply_blocklist(variant_df=variant_df, blocklist_uri=blocklist_uri,tiledb_cfg=tiledb_cfg)

        if bed_lines is not None and len(bed_lines) > 0:
            # Convert BED lines to BED regions
            bed_regions = [
                f"{line.split()[0]}:{line.split()[1]}-{line.split()[2]}"
                for line in bed_lines
            ]
            log.info(f"there are {len(bed_regions)} bed regions")
            variants_in_bed_region = []
            for region in bed_regions:
                # echo region being like chr1:0-999
                # BED is 0-based non-inclusive.
                # Don't touch the start, decrement the end
                chrom, start_pos, end_pos = (
                    region.split(":")[0],
                    int(region.split(":")[1].split("-")[0]) + 1,
                    int(region.split(":")[1].split("-")[1]),
                )
                # isolate variants from variant_df that are in the bed region

                variants_in_bed_region += [
                    variant_df[
                        (variant_df["CHROM"] == chrom)
                        & (variant_df["POS"] >= start_pos)
                        & (variant_df["POS"] <= end_pos)
                    ]
                ]
            variant_df = pd.concat(variants_in_bed_region)
        else:
            log.info(f"no bed regions found")

        # TODO: let's build this out of the CHROM POS REF ALT from moaavt
        variant_loci = variant_df.apply(lambda row: f"{row['CHROM']}-{row['POS']}-{row['REF']}-{row['ALT']}", axis=1).tolist()
        #variant_loci = list(variant_df["ID"])

        log.info(
            f"Querying genotype_query_by_variants for {len(variant_loci)} variants {variant_loci[:5]}..."
        )
        log.info(f"genotype_query_by_variants")
        t_start = time.time()

        #   contig  pos_start	ref	                                            alt	sample_name	zygosity	allele_string count
        # 0	chr1	1041678	CGCTCCGGCCAGTGCCAGGGTCGAGGTGAGCGGCTCCCCCGGGGGAGG	C	sample_1	HET	    T   1
        # 1	chr1	1041678	CGCTCCGGCCAGTGCCAGGGTCGAGGTGAGCGGCTCCCCCGGGGGAGG	C	sample_2	HET	    G/C 1
        # 2	chr1	1041678	CGCTCCGGCCAGTGCCAGGGTCGAGGTGAGCGGCTCCCCCGGGGGAGG	C	sample_3	HET	    T   1
        # 3	chr1	1041678	CGCTCCGGCCAGTGCCAGGGTCGAGGTGAGCGGCTCCCCCGGGGGAGG	C	sample_4	HET	    T/G 1
        # 4	chr1	1041678	CGCTCCGGCCAGTGCCAGGGTCGAGGTGAGCGGCTCCCCCGGGGGAGG	C	sample_5	HET     T/G 1
        # ...	...	...	...	...	...	...	...
        # 25995	chrX	154930903	C	T	sample_10	HEMI	1
        # 25996	chrX	154930950	G	C	sample_11	HET	1
        # 25997	chrX	154966103	C	T	sample_12	HEMI	1
        # 25998	chrX	154993141	T	G	sample_13	HET	1
        # 25999	chrX	154993141	T	G	sample_14	HET	1
        table = genotype_query_by_variants(
            vcf_uri,
            variants=variant_loci,
            allele_string_mode=True,
            zygosity_mode=True,
            max_sample_batch_size=2000,
            use_large_node=False,
            verbose=verbose,
            drop_samples=False,
            num_region_partitions=1,
            tiledb_cfg=tiledb_cfg,
        ).to_pandas()


        log.info(f"genotype_query_by_variants finished")


        if len(table) == 0:
            log.info(f"No variants of interest found in the TileDB-VCF dataset")
            return pa.Table.from_pandas(table)
        # drop count, not relevant if we have preserved the sample_name
        table = table.drop("count", axis=1)
        # contig  pos_start ref      alt zygosity  count
        variant_df_core = variant_df[
            ["CHROM", "POS", "REF", "ALT", "GENE", "CLINVAR_ID"]
        ]
        # rename to contig, pos_start, ref, alt
        variant_df_core = variant_df_core.rename(
            columns={"CHROM": "contig", "POS": "pos_start", "REF": "ref", "ALT": "alt"}
        )
        table = table.merge(
            variant_df_core,
            on=["contig", "pos_start", "ref", "alt"],
            how="inner",
        )
        # table = pa.Table.from_pandas(table)
        # split table by GENE and apply compute_comopund_heterozygosity to each group
        genes = table["GENE"].unique()

        cmpd_tables = []
        # TODO: in parallel either threaded or distributed
        for gene in genes:
            gene_table = table[table["GENE"] == gene]
            log.info(f"{len(gene_table)} variants found for gene {gene}")
            cmpd_tables.append(
                compute_compound_heterozygosity(pa.Table.from_pandas(gene_table))
            )

        table = pa.concat_tables(cmpd_tables)

        log.info(f"Done in {time.time() - t_start:.3f} seconds")

        vcf_df = table.to_pandas()
        
        if sample_metadata_uri is not None:
            filtered_metadata = fetch_approved_sample_metadata(
                sample_metadata_uri, sample_use_uri, population_source=population_source
            )
            log.info(f"filtered metadata: {len(filtered_metadata)} samples form filtered_metadata")
            # Convert sample_name to string if it's binary
            filtered_metadata['sample_name'] = filtered_metadata['sample_name'].astype(str)
            #remove -1 from sample_name
            if population_source == "rady":
                log.info("removing sample name suffix for merge")
                vcf_df["sample_name"] = vcf_df["sample_name"].str.replace(r"-1$", "", regex=True)
                filtered_metadata["sample_name"] = filtered_metadata["sample_name"].str.replace(r"-1$", "", regex=True)
            sm_samples = list(filtered_metadata["sample_name"])
            log.info(f"{len(sm_samples)} approved samples found in sample metadata")
            log.info(f"{len(vcf_df[['sample_name','contig','pos_start','ref','alt']].drop_duplicates())} tuples before sample filter")
            log.info(f"merging on sample_name {vcf_df['sample_name'][:5]} vs {filtered_metadata['sample_name'][:5]}")
            vcf_df = pd.merge(
                vcf_df, filtered_metadata, on=["sample_name"], how="inner"
            )
            log.info(f"{len(vcf_df[['sample_name','contig','pos_start','ref','alt']].drop_duplicates())} tuples after sample filter")



        # TODO: open this using the UDF auth

        if population_source == "rady":
            sex_field = "gender"
            sex_lookup = {"Male":"Male", "Female":"Female"}
        elif population_source == "onekg":
            sex_field = "sexf"
            sex_lookup = {"Male":1, "Female":2}
        elif population_source == "ukbb-alexion":
            sex_field = "sex_encoded"
            sex_lookup = {"Male":1, "Female":2}
        else:
            raise ValueError(f"Unknown population source: {population_source}")

        if moi_uri is not None:
            with tiledb.open(moi_uri) as moi_fh:
                moi_df = moi_fh.df[:]

            gene_poi_map = {}
            for gene, poi in zip(moi_df["GENE"], moi_df["MOI"]):
                gene_poi_map[gene.upper()] = poi

            vcf_df["moi"] = vcf_df["GENE"].map(gene_poi_map)
            vcf_df["positive_genotype"] = vcf_df.apply(positive_genotype, axis=1)
        else:
            vcf_df["positive_genotype"] = "N/A"
        
        #actually dedup on sample_name, contig, pos_start, ref, alt
        vcf_df = vcf_df.drop_duplicates(subset=['sample_name', 'contig', 'pos_start', 'ref', 'alt'])
        # count distinct sample names for each contig, pos_start, ref, alt

        vcf_df['chrposrefalt'] = vcf_df.apply(lambda row: f"{row['contig']}:{row['pos_start']}:{row['ref']}>{row['alt']}", axis=1)


        if aggregate_variants:
            log.info(f"Group by variant")
            log.info(f"Set zygosity for compound_event")
            # Set zygosity to CMPD_HET where compound_event is CMPD_HET
            vcf_df.loc[vcf_df['compound_event'] == 'CMPD_HET', 'zygosity'] = 'CMPD_HET'

            # Check for duplicate index entries before pivoting
            if vcf_df.duplicated(subset=['sample_name', 'chrposrefalt']).any():
                duplicates = vcf_df[vcf_df.duplicated(subset=['sample_name', 'chrposrefalt'], keep=False)]
                log.info(f"Found {len(duplicates)} duplicate rows:")
                log.info(duplicates)
                raise ValueError("Found duplicate sample_name/variant combinations - cannot pivot. Please deduplicate first.")
            
            # Get unique genes
            unique_genes = sorted(vcf_df['GENE'].unique())
            log.info(f"Found {len(unique_genes)} unique genes")

            # Create empty list to store gene-specific pivot tables
            all_gene_pivot_counts = []
            # Iterate over each gene
            for gene in unique_genes:
                gene_df = vcf_df[vcf_df['GENE'] == gene]
                log.info(f"Processing gene {gene} with {len(gene_df)} variants")
                
                # Create pivot table for this gene
                gene_pivot = gene_df.pivot(index='sample_name',
                                         columns='chrposrefalt',
                                         values='zygosity')
                gene_pivot = gene_pivot.fillna('')
                gene_pivot['variant_string'] = gene_pivot.apply(lambda row: ';'.join([f"{pos}({zyg})" if (zyg == 'HOM' or zyg == 'HEMI') else pos 
                                                        for pos, zyg in row.items() if zyg != '']), axis=1)
                log.info(f"{gene_pivot.head()}")
                
                gene_pivot_counts = gene_pivot.groupby('variant_string').size().reset_index(name='count')

                gene_pivot_counts = gene_pivot_counts.sort_values('count', ascending=False)
                gene_pivot_counts['GENE'] = gene
                all_gene_pivot_counts.append(gene_pivot_counts)

            pivot_df_counts = pd.concat(all_gene_pivot_counts)

            #TODO: we need to get an aggregation strategy to display positives in diplotypes
            #if moi_uri is not None:
                # Count positive genotypes per variant
            #    positive_counts = vcf_df.groupby('chrposrefalt')['positive_genotype'].apply(lambda x: (x == 'Yes').sum()).to_dict()
            #    pivot_as_df['positive_count'] = pivot_as_df.columns.map(positive_counts).fillna(0)
                
            
            vcf_df_zygosity = (
                vcf_df.groupby(
                    [
                        "contig", 
                        "pos_start",
                        "ref",
                        "alt",
                        "GENE",
                        "CLINVAR_ID",
                        "zygosity"
                    ]
                )["sample_name"]
                .nunique()
                .reset_index(name="sample_count")
            )

            vcf_df_positivity = (
                vcf_df.groupby(
                    [
                        "contig",
                        "pos_start",
                        "ref",
                        "alt",
                        "positive_genotype"
                    ]
                )["sample_name"]
                .nunique()
                .reset_index(name="sample_count")
            )


            log.info(f"pivoting by zygosity")
            vcf_df_zygosity_pivoted = vcf_df_zygosity.pivot_table(
                index=["contig", "pos_start", "ref", "alt", "GENE", "CLINVAR_ID"],
                columns="zygosity",
                values="sample_count",
                fill_value=0
            ).reset_index()
            


            log.info(f"pivoting by positivity")
            vcf_df_positivity_pivoted = vcf_df_positivity.pivot_table(
                index=["contig", "pos_start", "ref", "alt"],
                columns="positive_genotype",
                values="sample_count",
                fill_value=0
            ).reset_index()
            vcf_df_positivity_pivoted.rename(columns={"Yes":"Positive"},inplace=True)
            vcf_df_positivity_pivoted.drop(columns=['No'],inplace=True)

            # Check for missing columns
            required_columns = ["contig", "pos_start", "ref", "alt"]
            for col in required_columns:
                if col not in vcf_df_zygosity_pivoted.columns:
                    print(f"Missing column in vcf_df_zygosity_pivoted: {col}")
                if col not in vcf_df_positivity_pivoted.columns:
                    print(f"Missing column in vcf_df_positivity_pivoted: {col}")

            # Perform the merge
            log.info(f"merge zygosity and positivity")
            vcf_df = pd.merge(
                vcf_df_zygosity_pivoted,
                vcf_df_positivity_pivoted,
                on=required_columns,
                how="left"
            )

        @dataclass
        class BeginNGSResult:
            genotypes: pa.Table
            diplotypes: pa.Table
        genotypes = pa.Table.from_pandas(vcf_df)
        diplotypes = pa.Table.from_pandas(pivot_df_counts)
        return BeginNGSResult(genotypes=genotypes, diplotypes=diplotypes)
    
def test_beginngs():


    
    df = beginngs_query(
        vcf_uri="tiledb://XXXXXXXXXXXXXXXX",
        moaavt_uri="tiledb://XXXXXXXXXXXXXXXX",
        moi_uri="tiledb://XXXXXXXXXXXXXXXX",
        blocklist_uri="tiledb://XXXXXXXXXXXXXXXX",
        sample_metadata_uri="tiledb://XXXXXXXXXXXXXXXX",
        sample_use_uri="tiledb://XXXXXXXXXXXXXXXX",
        use_blocklist=True,
        bed_lines=None,
        tiledb_cfg=None
    ).to_pandas()
    print(df)

if __name__ == "__main__":
    test_beginngs()


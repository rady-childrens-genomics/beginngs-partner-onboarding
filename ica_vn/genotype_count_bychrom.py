import io
import os
import gzip
import csv
import argparse
import pandas as pd
from collections import Counter

from warnings import simplefilter

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

def get_base_filename(file_path):
    filename = os.path.basename(file_path)
    chr_name = filename.replace('_filtered.vcf', '')
    return chr_name

def get_vcf_names(vcf_path):
    with open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x.rstrip() for x in line.split('\t')]
                  break
    ifile.close()
    return vcf_names

def zygosity(gt: str) -> str:
    gt = str(gt)
    if "/" in gt:
        gt = gt.split("/")

    # All genotypes are missing
    if len(gt) == 0 or all(allele == "." for allele in gt):
        return "MISSING"

    # One allele
    if len(gt) == 1:
        if gt[0] == "0":
            return "HOM_REF"
        return "HEMI"

    # More than one allele
    if all(allele == "0" for allele in gt):
        return "HOM_REF"
    if all(allele == gt[0] for allele in gt):
        return "HOM_ALT"
    return "HET"

def categorize_variant(zygosity, moi, sex=0):
    if moi == "":
        return "Unknown"
    elif moi == "AD":
        if zygosity != "HOM_REF":
            return "Positive"
    elif moi == "AR":
        if zygosity == "HOM_ALT" or zygosity == "CMPD_HET":
            return "Positive"
        else:
            return "Negative"
    elif moi == "XR":
        if sex == 2:
            if zygosity == "HOM_ALT" or zygosity == "CMPD_HET":
                return "Positive"
        else:
            if zygosity == "HEMI":
                return "Positive"
        return "Negative"
    elif moi == "XD":
        if sex == 1:
            if zygosity == "HEMI":
                return "Positive"
        else:
            if zygosity != "HOM_REF":
                return "Positive"

def process_variant(gene_df, variants_per_sample, variant_count, moi_value):
    for sample in variants_per_sample.keys():
        sample_het = []
        sample_cmpd_het = []
        sample_homalt = []
        sample_hemi = []
        sample_sex = sex_dict.get(sample)
        genotype_per_sample = ""
        for _, row in gene_df.iterrows():
            genotype = row[sample]
            #getting non ref variants for the sample
            if genotype == 'HOM_ALT':
                sample_homalt.append(row['ID'])
            elif genotype == 'CMPD_HET':
                sample_cmpd_het.append(row['ID'])
            elif genotype == 'HET':
                sample_het.append(row['ID'])
            elif genotype == 'HEMI':
                sample_hemi.append(row['ID'])

        if (moi_value == 'AD' or (moi_value == 'XD' and sample_sex == 2)):
            if sample_het: #there should be just one het variant per sample in the gene, otherwise it would be a cmpd het
                genotype_per_sample = f"[{sample_het[0]}];+,HETEROZYGOUS,{row['GENE']},{moi_value}"
                #print(genotype_per_sample)
                variants_per_sample[sample].append(genotype_per_sample)
                variant_count[genotype_per_sample] += 1
        if moi_value == 'AR' or moi_value == 'AD' or ((moi_value == 'XR' or moi_value == 'XD') and sample_sex == 2): #not mutually exclusive
            if not sample_homalt and sample_cmpd_het:
                variants = ';'.join([f"[{item}]" for item in sample_cmpd_het])
                #print(variants)
                genotype_per_sample = f"{variants},COMPOUND HETEROZYGOUS,{row['GENE']},{moi_value}"
                variants_per_sample[sample].append(genotype_per_sample)
                variant_count[genotype_per_sample] += 1 
            elif sample_homalt:
                variants = ';'.join(sample_homalt)
                #print(variants)
                if sample_cmpd_het:
                    variants_cmpd = ';'.join(sample_cmpd_het)
                    genotype_per_sample = f"[{variants}];[{variants};{variants_cmpd}],HOMOZYGOUS,{row['GENE']},{moi_value}"
                else:
                    genotype_per_sample = f"[{variants}];[{variants}],HOMOZYGOUS,{row['GENE']},{moi_value}"
                variants_per_sample[sample].append(genotype_per_sample)
                variant_count[genotype_per_sample] += 1
        elif (moi_value == 'XR' or moi_value == 'XD') and sample_sex == 1:
                if sample_hemi:
                    variants = ';'.join([f"[{item}]" for item in sample_hemi])
                    genotype_per_sample = f"{variants},HEMIZYGOUS,{row['GENE']},{moi_value}"
                    variants_per_sample[sample].append(genotype_per_sample)
                    variant_count[genotype_per_sample] += 1
            
                    

def compute_diplotype (df_pos_only_perchrom):
    columns_to_drop = [col for col in df_pos_only_perchrom.columns if col not in ["CHROM","POS","REF","ALT","ID","QUAL","FILTER", \
        "INFO", "FORMAT", "CLINVAR_VARIANT_ID", "GENE", "cmpd_het_count_by_variant", "het_count_by_variant", "hom_alt_count_by_variant", "hemi_count_by_variant", "MOI", "Positive"] and (df_pos_only_perchrom[col] == 'HOM_REF').all()]

    # Drop those columns
    df_pos_only_perchrom = df_pos_only_perchrom.drop(columns=columns_to_drop)

    #Empty dictionary to store the variants per sample
    variants_per_sample = {col: [] for col in df_pos_only_perchrom.columns if col not in ["CHROM","POS","REF","ALT","ID","QUAL","FILTER", \
            "INFO", "FORMAT", "CLINVAR_VARIANT_ID", "GENE", "cmpd_het_count_by_variant", "het_count_by_variant", "hom_alt_count_by_variant", "hemi_count_by_variant", "MOI", "Positive"]}

    #Counter to count occurrences of variants
    variant_count = Counter()


    #Process samples by gene
    for gene, gene_df in df_pos_only_perchrom.groupby('GENE'):
        #cmpd_het_variants = []
        moi_value = gene_df['MOI'].iloc[0]

        if moi_value == 'AR' or moi_value == 'AD' or moi_value == 'XR' or moi_value == 'XD':   ##remove if
            process_variant(gene_df, variants_per_sample, variant_count, moi_value)
        
    return(variant_count)

def compute_compound_heterozygous(vcf_df, sample_columns):
        #Dictionary to track counts of HET variants per sample per gene
        het_counts = {}

        # Iterate over each gene and sample -> time consuming step
        for sample in vcf_df[sample_columns]:
            for gene in vcf_df['GENE'].unique():
                # Count the number of HET variants for the specific sample and gene
                count = len(vcf_df[(vcf_df['GENE'] == gene) & (vcf_df[sample] == 'HET')])
                count_hom = len(vcf_df[(vcf_df['GENE'] == gene) & (vcf_df[sample] == 'HOM_ALT')])
                #print("count is ", count)
                if count > 1 or (count >0 and count_hom>0):   ##checking if hom alt and het calls are present in the same gene:
                    if sample not in het_counts:
                        het_counts[sample] = {}
                    het_counts[sample][gene] = [count, count_hom]
        return(het_counts)

def convert_maleHET_to_HEMI(vcf_df, sample_columns): #to fix issue where male HET samples in X should be HEMI
    for sample in vcf_df[sample_columns]:
        sample_sex = sex_dict.get(sample)
        if sample_sex == 1:
            vcf_df.loc[(vcf_df[sample] == 'HET'), sample] = 'HEMI'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filtered_vcf_file_input", type=str, required=True,
        help="Path to the vcf input file")
    
    parser.add_argument(
        "--voi", type=str, required=True,
        help="Path to the voi input file")
    
    parser.add_argument(
        "--output", type=str, required=True,
        help="Output path")
    
    parser.add_argument(
        "--moi", type=str, required=True,
        help="Mode of inheritance file")
    
    parser.add_argument(
        "--sex", type=str, required=True,
        help="Metadata file")
    
    parser.add_argument(
        "--voi_annotated", type=str, required=True,
        help="Annotated voi")
    
    args = parser.parse_args()
    vcf_file_path = args.filtered_vcf_file_input
    voi = args.voi
    output_path = args.output
    input_moi = args.moi
    input_sex = args.sex
    voi_annotated = args.voi_annotated

    chr = get_base_filename(vcf_file_path)

    #get column names
    names = get_vcf_names(vcf_file_path)
    vcf = pd.read_csv(vcf_file_path, comment='#', sep="\t", header=None, names=names)

    variants = pd.read_csv(voi,sep = '\t', names=["#CHROM", "POS", "REF", "ALT"])
    variants_annotated = pd.read_csv(voi_annotated,sep = ',', usecols=["CHROM", "POS", "CLINVAR_VARIANT_ID", "REF", "ALT", "GENE"])

    # Merge to select only variants with selected alt
    vcf_df = pd.merge(variants, vcf, on=["#CHROM", "POS", "REF", "ALT"], how="inner")

    # Calculate zygosity of each variant from genotype
    for column in vcf_df:
        if column not in (["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]):
            vcf_df[column] = vcf_df[column].apply(zygosity)

    vcf_df = vcf_df.rename(columns={"#CHROM": "CHROM"})

    # Add gene annotation
    vcf_df = vcf_df.merge(variants_annotated, on=["CHROM", "POS", "REF", "ALT"], how="left")

    # For postproc vcfs, get genotype columns (samples):
    gt_columns = [col for col in vcf_df.columns if col not in (["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "CLINVAR_VARIANT_ID", "GENE"])]

    het_counts = compute_compound_heterozygous(vcf_df, gt_columns)

    # Modify the original DataFrame based on the het_counts information
    for sample in vcf_df[gt_columns]:
        for gene, counts in het_counts.get(sample, {}).items():
            if counts[0] > 1 or (counts[0] >0 and counts[1]>0):
                # Update the rows where the sample is "HET" for the current gene to "CMPD_HET"
                vcf_df.loc[(vcf_df['GENE'] == gene) & (vcf_df[sample] == 'HET'), sample] = 'CMPD_HET'

    # Merge with mode of inheritance
    moi = pd.read_csv(input_moi,sep = '\t', names=["GENE", "MOI"], skiprows=1)

    vcf_df = vcf_df.merge(moi, on=["GENE"], how="left")

    # Get gender information from metadata file
    sex = pd.read_csv(input_sex, sep = ' ', usecols=["SampleID", "Sex"])

    global sex_dict
    sex_dict = sex.set_index('SampleID')['Sex'].to_dict()
    
    # Need to convert the X HET variants in male samples into HEMI variants before counting
    for _, gene_df in vcf_df.groupby('GENE'):
        moi_value = gene_df['MOI'].iloc[0]
        if moi_value == 'XD':
            convert_maleHET_to_HEMI(vcf_df, gt_columns)

    # Count types of occourring genotypes for each variant
    vcf_df["cmpd_het_count_by_variant"] = vcf_df[gt_columns].apply(lambda x: x.str.contains("CMPD_HET")).sum( axis=1)

    vcf_df["het_count_by_variant"] = vcf_df[gt_columns].apply(lambda x: x.str.contains("^HET")).sum( axis=1)

    vcf_df["hom_alt_count_by_variant"] = vcf_df[gt_columns].apply(lambda x: x.str.contains("HOM_ALT")).sum( axis=1)

    vcf_df["hemi_count_by_variant"] = vcf_df[gt_columns].apply(lambda x: x.str.contains("^HEMI")).sum( axis=1)    

    # Create empty positives table
    positives_table = pd.DataFrame(index=vcf_df.index, columns=gt_columns)

    # Populate positives table
    for variant in vcf_df.index:
        inheritance_mode = vcf_df.loc[variant, 'MOI']  # Get inheritance mode for current variant
        for sample in gt_columns:
            sample_sex = sex_dict.get(sample)

            positives_table.loc[variant, sample] = categorize_variant(vcf_df.loc[variant, sample], inheritance_mode , sample_sex)

    # Count positives for each variant
    positives_table["Positive"] = positives_table[gt_columns].apply(lambda x: x.str.contains("Positive")).sum( axis=1)

    # Save genotypes into csv output table

    output_genotypes_1 = vcf_df[['CHROM', 'POS', 'REF', 'ALT', 'GENE', 'CLINVAR_VARIANT_ID', 'cmpd_het_count_by_variant', 'hemi_count_by_variant' ,'het_count_by_variant', 'hom_alt_count_by_variant']]

    output_genotypes_2 = positives_table[["Positive"]]

    output_genotypes = pd.concat([output_genotypes_1,output_genotypes_2], axis=1)

    output_genotypes.to_csv(f'{output_path}/{chr}_genotypes.csv', index=False)

    vcf_df["Positive"] = positives_table["Positive"]

    # Filter by only positive variants
    df_filtered = vcf_df[vcf_df.iloc[:,-1] >0]

    diplotypes = compute_diplotype(df_filtered)

    
    with open (f'{output_path}/{chr}_diplotypes.csv', mode='w', newline='') as file:
        header = ["DYPLOTYPE", "ZYGOSITY", "GENE", "MOI", "COUNT"]
        writer = csv.writer(file)
        writer.writerow(header)
        for key, count in diplotypes.items():
            values = key.split(',')
            row = values + [count]
            writer.writerow(row)


if __name__ == "__main__":
    main()
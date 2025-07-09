#!/bin/bash

path="../../1000Genomes/ms-vcf/post-proc-msvcf/"
moi="input_files/moi_20240805.txt"
sex="input_files/20130606_g1k_3202_samples_ped_population.txt"
voi_annotated="input_files/variants_of_interest_20231108.csv"
blocklist="input_files/blocklist_20241114.csv"

voi="input_files/voi_simple.tsv"
tail -n +2 "$voi_annotated" | csvcut -c 3,4,6,7 | csvformat -T > "${voi}"

if [ -n "$blocklist" ]; then
  grep -v -f <(awk -F'\t' '{print "^" $2 "\t" $3 "\t" $4 "\t" $5 "$"}' "${blocklist}") "${voi}" > "input_files/filtered_voi.tsv"
  voi="input_files/filtered_voi.tsv"
fi

mkdir -p filtered_postproc/
mkdir -p output_genotypes/

###Launch bcftools and python script on all input files in the path###

for folder in "${path}"*; do
    chrom_folder=$(basename "${folder}")
        for file in "$folder"/postproc-vcf/*
        do if [[ $file == *dragen.vcf.gz ]]
        then
            filtered_vcf="filtered_postproc/${chrom_folder}_filtered.vcf"
            (
            bcftools view -O v -R "$voi" "$file" > "$filtered_vcf"
            python3 genotype_count_bychrom.py --filtered_vcf_file_input "$filtered_vcf" --voi "$voi" \
            --output output_genotypes --moi "$moi" --sex "$sex" --voi_annotated "$voi_annotated"
                
            ) &
            pids+=($!)
        fi
        done
done


for pid in "${pids[@]}"; do
    wait "$pid"
done

echo "done"


filenum=0
for file in output_genotypes/*genotypes.csv
    do if [[ $filenum == 0 ]]
    then
            cat "$file" > genotypes.csv
    else
    awk 'NR>1' "$file" >> genotypes.csv
    fi
    filenum=$((filenum + 1))
done



filenum=0
for file in output_genotypes/*diplotypes.csv
    do if [[ $filenum == 0 ]]
    then
            cat "$file" > diplotypes.csv
    else
    awk 'NR>1' "$file" >> diplotypes.csv
    fi
    filenum=$((filenum + 1))
done

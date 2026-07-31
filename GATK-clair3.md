**GATK-"Only retain biallelic SNPs."

```
bcftools view -m2 -M2 -v snps run1_ARS_UCD_v2_0.genomeWide.filteredpass.SNPs.vcf.gz -o run1_ARS_biallelic_SNPs.vcf.gz
```

**"Split the 20 samples."

```
samples=(
  3DWF22T0036 3DWF22T0070 3DWF22T0240 3DWF22T0415 3DWF22T0549
  3DWFN0018 3DWFP0355 3DWFP0398 3DWFP0413 3DWFQ0177
  3DWFQ0404 3DWFQ0458 3DWFQ0459 3DWFR0017 3DWFR0281
  3DWFS0069 3DWFS0123 3DWFS0137 3DWFS0202 WGWFL0149
)

for sample in "${samples[@]}"; do
  echo "sample processing: $sample"
  
##Extract samples and compress the output.
  bcftools view -s "$sample" run1_ARS_biallelic_SNPs.vcf.gz -o "${sample}.vcf.gz"

##Create index
  bcftools index "${sample}.vcf.gz"
done

echo "done！"
```

**Filtering. Strict filtering criteria:** DP between 20 and 40, GQ > 30; **relaxed filtering criteria:** DP between 10 and 40, GQ > 20.

```
while read sample; do
  if [ -f "${sample}.vcf.gz" ]; then
    echo "Processing ${sample}.vcf.gz..."
    bcftools filter -i 'FORMAT/DP>=20 && FORMAT/DP<=40 && FORMAT/GQ>30' \
      "${sample}.vcf.gz" -o "${sample}.filtered-DP20-40-GQ30.vcf.gz" -O z
    tabix "${sample}.filtered-DP20-40-GQ30.vcf.gz"
  fi
done < sample.list
while read sample; do
  if [ -f "${sample}.vcf.gz" ]; then
    echo "Processing ${sample}.vcf.gz..."
    bcftools filter -i 'FORMAT/DP>=10 && FORMAT/DP<=40 && FORMAT/GQ>20' \
      "${sample}.vcf.gz" -o "${sample}.filtered-DP10-40-GQ20.vcf.gz" -O z
    tabix "${sample}.filtered-DP10-40-GQ20.vcf.gz"
  fi
done < sample.list
```

**Using gtcheck for comparison**

```
while read s1 s2; do 
    echo "Comparing: $s1 vs $s2"
   bcftools gtcheck -g /storage/public/home/2024050455/08.compare/Extract/${s1}.filtered-DP10-40-GQ20.vcf.gz /storage/public/home/2024050455/08.compare/Sample/${s2}_DP6-20_GQ10.vcf.gz | awk -v p="${s1}_vs_${s2}" '{print p "\t" $0}'
done < /storage/public/home/2024050455/08.compare/Compare/sample_pair.list > all_concordance-gDP10-40-GQ20-cDP6-20-GQ10.txt
while read s1 s2; do
    echo "Comparing: $s1 vs $s2"
    bcftools gtcheck -g /storage/public/home/2024050455/08.compare/Extract/${s1}.filtered-DP10-40-GQ20.vcf.gz /storage/public/home/2024050455/08.compare/Sample/${s2}_DP10-20_GQ20.vcf.gz | awk -v p="${s1}_vs_${s2}" '{print p "\t" $0}'
done < /storage/public/home/2024050455/08.compare/Compare/sample_pair.list > all_concordance-gDP10-40-GQ20-cDP10-20-GQ20.txt
while read s1 s2; do
    echo "Comparing: $s1 vs $s2"
bcftools gtcheck -g /storage/public/home/2024050455/08.compare/Extract/${s1}.filtered-DP20-40-GQ30.vcf.gz /storage/public/home/2024050455/08.compare/Sample/${s2}_DP6-20_GQ10.vcf.gz| awk -v p="${s1}_vs_${s2}" '{print p "\t" $0}'
done < /storage/public/home/2024050455/08.compare/Compare/sample_pair.list > all_concordance-gDP20-40-GQ30-cDP6-20-GQ10.txt
while read s1 s2; do
    echo "Comparing: $s1 vs $s2"
bcftools gtcheck -g /storage/public/home/2024050455/08.compare/Extract/${s1}.filtered-DP20-40-GQ30.vcf.gz /storage/public/home/2024050455/08.compare/Sample/${s2}_DP10-20_GQ20.vcf.gz| awk -v p="${s1}_vs_${s2}" '{print p "\t" $0}'
done < /storage/public/home/2024050455/08.compare/Compare/sample_pair.list > all_concordance-gDP20-40-GQ30-cDP10-20-GQ20.txt
```

**Output:**  
Number of sites compared, discordance, inconsistency rate (discordant / sites compared), consistency rate (1 − inconsistency), number of sites skipped due to no match, and overlap fraction (sites compared / (sites compared + skipped)).

```
#!/bin/bash

files=(
all_concordance-gDP10-40-GQ20-cDP10-20-GQ20.txt
all_concordance-gDP10-40-GQ20-cDP6-20-GQ10.txt
all_concordance-gDP20-40-GQ30-cDP10-20-GQ20.txt
all_concordance-gDP20-40-GQ30-cDP6-20-GQ10.txt
)

for f in "${files[@]}"; do
    out=${f/all_concordance/gtcheck_summary}
    out=${out/.txt/.tsv}

    echo "Processing $f -> $out"

    awk '
    {
        pair = $1

        # Sites_Compared / Sites_Skipped
        if ($3 == "sites-compared") sc[pair] = int($4)
        if ($3 == "sites-skipped-no-match") ss[pair] = int($4)

        # Discordance：DC 行，第5列
        if ($2 == "DC") dc[pair] = int($5)
    }
    END {
        print "Sample_Pair\tSites_Compared\tSites_Skipped\tDiscordance\tInconsistency_Rate\tConsistency_Rate\tOverlap_Fraction"
        for (p in sc) {
            sc_i = sc[p] + 0
            ss_i = ss[p] + 0
            dc_i = (p in dc) ? dc[p] : 0

            inc = (sc_i > 0) ? dc_i / sc_i : 0
            con = 1 - inc
            ov  = (sc_i + ss_i > 0) ? sc_i / (sc_i + ss_i) : 0

            printf "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\n",
                   p, sc_i, ss_i, dc_i, inc, con, ov
        }
    }
    ' "$f" > "$out"

done

echo "All files processed."

```

**calculate unique

```
#!/bin/bash
set -euo pipefail

# =========================

PAIR_LIST="/storage/public/home/2024050455/08.compare/Compare/sample_pair.list"
GATK_DIR="/storage/public/home/2024050455/08.compare/Extract"
CLAIR3_DIR="/storage/public/home/2024050455/08.compare/Sample"
OUT_DIR="/storage/public/home/2024050455/08.compare/Compare/unique_snp"
mkdir -p "$OUT_DIR"

# 
declare -A COMPARISONS=(
  ["gDP10-40-GQ20-cDP6-20-GQ10"]="filtered-DP10-40-GQ20.vcf.gz _DP6-20_GQ10.vcf.gz"
  ["gDP10-40-GQ20-cDP10-20-GQ20"]="filtered-DP10-40-GQ20.vcf.gz _DP10-20_GQ20.vcf.gz"
  ["gDP20-40-GQ30-cDP6-20-GQ10"]="filtered-DP20-40-GQ30.vcf.gz _DP6-20_GQ10.vcf.gz"
  ["gDP20-40-GQ30-cDP10-20-GQ20"]="filtered-DP20-40-GQ30.vcf.gz _DP10-20_GQ20.vcf.gz"
)

# =========================

for tag in "${!COMPARISONS[@]}"; do
    echo "=== Processing $tag ==="

    read GATK_SUFFIX CLAIR3_SUFFIX <<< "${COMPARISONS[$tag]}"

    OUT_FILE="$OUT_DIR/unique_snp_${tag}.tsv"
    echo -e "Sample_Pair\tSNP-unique-GATK\tSNP-unique-Clair3" > "$OUT_FILE"

    while read -r s1 s2; do
        pair="${s1}_vs_${s2}"

        GATK_VCF="$GATK_DIR/${s1}.${GATK_SUFFIX}"
        CLAIR3_VCF="$CLAIR3_DIR/${s2}${CLAIR3_SUFFIX}"

        if [[ ! -f "$GATK_VCF" || ! -f "$CLAIR3_VCF" ]]; then
            echo "  [SKIP] $pair (VCF missing)"
            continue
        fi

        tmp=$(mktemp -d)

        # Extract SNP sites（CHROM + POS）
        bcftools view -v snps "$GATK_VCF" \
          | bcftools query -f '%CHROM\t%POS\n' \
          | sort -u > "$tmp/gatk.snps"

        bcftools view -v snps "$CLAIR3_VCF" \
          | bcftools query -f '%CHROM\t%POS\n' \
          | sort -u > "$tmp/clair3.snps"

        # calculate unique
        gatk_unique=$(comm -23 "$tmp/gatk.snps" "$tmp/clair3.snps" | wc -l)
        clair3_unique=$(comm -13 "$tmp/gatk.snps" "$tmp/clair3.snps" | wc -l)

        echo -e "$pair\t$gatk_unique\t$clair3_unique" >> "$OUT_FILE"

        rm -rf "$tmp"
        echo "  $pair done: GATK=$gatk_unique Clair3=$clair3_unique"

    done < "$PAIR_LIST"

    echo "Result saved to $OUT_FILE"
    echo
done

```

We obtained the file "ConcordanceAnalysisClair3vsGATK_20260116". From the table, we found that the filtering strategy "gDP20-40-GQ30 vs cDP10-20-GQ20" performed best.

**
```
#!/bin/bash
set -euo pipefail

# ============================================================================
# Script Name: GATK_vs_Clair3_Concordance_Analysis.sh
# Description: Compare variant calling results between GATK and Clair3
#             1. Extract common sites and compare genotypes
#             2. Generate concordance matrix (0=concordant, 1=discordant)
#             3. Generate detailed matrix (with genotype info for discordant sites)
#             4. Calculate discordance rate per site
#             5. Classify sites by quality (Good/Bad/Low_support)
#             6. Extract flanking sequences for high-quality SNPs
#             7. Visualization analysis (embedded Python script)
# Author: Xuewei
# ============================================================================

# ===============================
# Configuration Parameters (modify as needed)
# ===============================
GATK_PATH="/storage/public/home/2024050455/08.compare/Extract"
CLAIR3_PATH="/storage/public/home/2024050455/08.compare/Sample"
REF_FASTA="ARS-UCD2.0.fa"          # Reference genome file path
OUTDIR="concordance-analysis"       # Output directory
TMPDIR="temp_files"                 # Temporary file directory

# Sample pair mapping (GATK sample name <-> Clair3 sample name)
# Index positions must correspond one-to-one
samples_gatk=(
    3DWF22T0036 3DWF22T0070 3DWF22T0240 3DWF22T0415 3DWF22T0549
    3DWFN0018 3DWFP0355 3DWFP0398 3DWFP0413 3DWFQ0177
    3DWFQ0404 3DWFQ0458 3DWFQ0459 3DWFR0017 3DWFR0281
    3DWFS0069 3DWFS0123 3DWFS0137 3DWFS0202 WGWFL0149
)

samples_clair=(
    sample15 sample06 sample12 sample11 sample16
    sample03 sample07 sample20 sample04 sample05
    sample01 sample13 sample02 sample10 sample17
    sample18 sample14 sample09 sample19 sample08
)

# ===============================
# Create output directories
# ===============================
mkdir -p "${OUTDIR}" "${TMPDIR}"

# ============================================================================
# Step 1: Extract common sites for all sample pairs and compare genotypes
# Input: GATK and Clair3 VCF files (per sample pair)
# Output: ${OUTDIR}/all_sites_concordance.txt
#         Contains detailed information and concordance status for each site
# ============================================================================
echo "=========================================="
echo "Step 1: Extracting common sites and comparing genotypes"
echo "=========================================="

OUTFILE="${OUTDIR}/all_sites_concordance.txt"

# Write header
echo -e "CHROM\tPOS\tREF\tALT\tGATK_Sample\tGATK_GT\tGATK_GQ\tGATK_DP\tClair3_Sample\tClair3_GT\tClair3_GQ\tClair3_DP\tDiscordant" \
    > "${OUTFILE}"

# Iterate through each sample pair
for i in "${!samples_gatk[@]}"; do
    gatk_sample="${samples_gatk[$i]}"
    clair_sample="${samples_clair[$i]}"
    
    gatk_vcf="${GATK_PATH}/${gatk_sample}.filtered-DP20-40-GQ30.vcf.gz"
    clair_vcf="${CLAIR3_PATH}/${clair_sample}_DP10-20_GQ20.vcf.gz"
    
    echo "Processing sample pair: ${gatk_sample} vs ${clair_sample}"
    
    # 1.1 Extract common sites from both VCFs (-n=2 requires both files have the site)
    bcftools isec -n=2 -w1 "${gatk_vcf}" "${clair_vcf}" -Oz \
        -o "${TMPDIR}/${gatk_sample}.common.vcf.gz"
    bcftools isec -n=2 -w2 "${gatk_vcf}" "${clair_vcf}" -Oz \
        -o "${TMPDIR}/${clair_sample}.common.vcf.gz"
    
    # 1.2 Create indexes
    bcftools index -f "${TMPDIR}/${gatk_sample}.common.vcf.gz"
    bcftools index -f "${TMPDIR}/${clair_sample}.common.vcf.gz"
    
    # 1.3 Extract key fields: CHROM, POS, REF, ALT, GT, GQ, DP
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ\t%DP]\n' \
        "${TMPDIR}/${gatk_sample}.common.vcf.gz" > "${TMPDIR}/${gatk_sample}.info"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ\t%DP]\n' \
        "${TMPDIR}/${clair_sample}.common.vcf.gz" > "${TMPDIR}/${clair_sample}.info"
    
    # 1.4 Merge information from both files and compare genotype concordance
    # Note: norm_gt function normalizes genotypes (e.g., 0/1 and 1/0 both become 0/1)
    paste "${TMPDIR}/${gatk_sample}.info" "${TMPDIR}/${clair_sample}.info" | \
    awk -v gs="${gatk_sample}" -v cs="${clair_sample}" '
    BEGIN { OFS="\t" }
    
    # Normalize genotype: replace | with / and sort alleles
    function norm_gt(gt, a, b) {
        gsub(/\|/, "/", gt)
        split(gt, arr, "/")
        a = arr[1]; b = arr[2]
        if (a > b) return b "/" a
        else return a "/" b
    }
    
    {
        # GATK fields
        gatk_gt = norm_gt($5)
        gatk_gq = $6
        gatk_dp = $7
        
        # Clair3 fields
        clair_gt = norm_gt($12)
        clair_gq = $13
        clair_dp = $14
        
        # Check concordance (0=concordant, 1=discordant)
        discord = (gatk_gt == clair_gt) ? 0 : 1
        
        print $1, $2, $3, $4,
              gs, gatk_gt, gatk_gq, gatk_dp,
              cs, clair_gt, clair_gq, clair_dp,
              discord
    }' >> "${OUTFILE}"
    
done

echo "Step 1 completed! Output file: ${OUTFILE}"
echo ""

# ============================================================================
# Step 2: Generate concordance matrix (wide format, values 0/1)
# Input: ${OUTDIR}/all_sites_concordance.txt
# Output: ${OUTDIR}/site_level_concordance_matrix.tsv
#         Rows=SNP sites, Columns=sample pairs, Values=0(concordant)/1(discordant)/NA(absent)
# ============================================================================
echo "=========================================="
echo "Step 2: Generating concordance matrix (0/1 format)"
echo "=========================================="

INPUT="${OUTDIR}/all_sites_concordance.txt"
MATRIX_OUT="${OUTDIR}/site_level_concordance_matrix.tsv"
TMP="${OUTDIR}/.body.tmp"

# 2.1 Generate header (all sample pair names)
awk -v OFS="\t" '
NR==1 { next }
{
    pairs[$5 "_vs_" $9] = 1
}
END {
    n = asorti(pairs, pair_list)
    printf "CHROM\tPOS\tREF\tALT"
    for (i = 1; i <= n; i++)
        printf "\t%s", pair_list[i]
    printf "\n"
}
' "${INPUT}" > "${MATRIX_OUT}"

# 2.2 Generate data body (concordance status for each site across sample pairs)
awk -v OFS="\t" '
NR==1 { next }

{
    site = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
    pair = $5 "_vs_" $9
    discord = $13     # 0=concordant, 1=discordant
    
    sites[site] = 1
    pairs[pair] = 1
    data[site, pair] = discord
}

END {
    n = asorti(pairs, pair_list)
    
    for (s in sites) {
        split(s, a, SUBSEP)
        printf "%s\t%s\t%s\t%s", a[1], a[2], a[3], a[4]
        
        for (i = 1; i <= n; i++) {
            key = s SUBSEP pair_list[i]
            if (key in data)
                printf "\t%s", data[key]
            else
                printf "\tNA"
        }
        printf "\n"
    }
}
' "${INPUT}" | sort -k1,1 -k2,2n > "${TMP}"

cat "${TMP}" >> "${MATRIX_OUT}"
rm -f "${TMP}"

echo "Step 2 completed! Output file: ${MATRIX_OUT}"
echo ""

# ============================================================================
# Step 3: Generate detailed matrix (show specific genotypes for discordant sites)
# Input: ${OUTDIR}/all_sites_concordance.txt
# Output: ${OUTDIR}/site_level_concordance_detailed.tsv
#         Same as Step 2, but discordant sites show GATK(GT,GQ,DP)|CLAIR3(GT,GQ,DP)
# ============================================================================
echo "=========================================="
echo "Step 3: Generating detailed concordance matrix"
echo "=========================================="

DETAILED_OUT="${OUTDIR}/site_level_concordance_detailed.tsv"

# 3.1 Generate header (same as Step 2)
awk -v OFS="\t" '
NR==1 { next }
{
    pairs[$5 "_vs_" $9] = 1
}
END {
    n = asorti(pairs, pair_list)
    printf "CHROM\tPOS\tREF\tALT"
    for (i = 1; i <= n; i++)
        printf "\t%s", pair_list[i]
    printf "\n"
}
' "${INPUT}" > "${DETAILED_OUT}"

# 3.2 Generate data body (output detailed genotype info for discordant sites)
awk -v OFS="\t" '
NR==1 { next }

{
    site = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
    pair = $5 "_vs_" $9
    
    # Output 0 for concordant, detailed genotype info for discordant
    if ($13 == 0) {
        cell = 0
    } else {
        cell = "GATK(" $6 "," $7 "," $8 ")|CLAIR3(" $10 "," $11 "," $12 ")"
    }
    
    sites[site] = 1
    pairs[pair] = 1
    data[site, pair] = cell
}

END {
    n = asorti(pairs, pair_list)
    
    for (s in sites) {
        split(s, a, SUBSEP)
        printf "%s\t%s\t%s\t%s", a[1], a[2], a[3], a[4]
        
        for (i = 1; i <= n; i++) {
            key = s SUBSEP pair_list[i]
            if (key in data)
                printf "\t%s", data[key]
            else
                printf "\tNA"
        }
        printf "\n"
    }
}
' "${INPUT}" | sort -k1,1 -k2,2n > "${TMP}"

cat "${TMP}" >> "${DETAILED_OUT}"
rm -f "${TMP}"

echo "Step 3 completed! Output file: ${DETAILED_OUT}"
echo ""

# ============================================================================
# Step 4: Calculate discordance rate per site
# Input: ${OUTDIR}/site_level_concordance_matrix.tsv
# Output: ${OUTDIR}/site_discordance_summary.tsv
#         Contains discordant count, concordant count, total comparisons, and discordance rate
# ============================================================================
echo "=========================================="
echo "Step 4: Calculating discordance rate per site"
echo "=========================================="

SUMMARY_OUT="${OUTDIR}/site_discordance_summary.tsv"

awk -v OFS="\t" '
NR==1 {
    print "CHROM","POS","REF","ALT",
          "Discordant_Pair_Count",
          "Concordant_Pair_Count",
          "Total_Compared_Pairs",
          "Discordance_Rate"
    next
}

{
    discord = 0
    concord = 0
    total = 0
    
    for (i = 5; i <= NF; i++) {
        if ($i == "NA") continue
        total++
        if ($i == 0)
            concord++
        else
            discord++
    }
    
    if (total > 0)
        rate = discord / total
    else
        rate = "NA"
    
    print $1, $2, $3, $4,
          discord,
          concord,
          total,
          rate
}
' "${MATRIX_OUT}" > "${SUMMARY_OUT}"

echo "Step 4 completed! Output file: ${SUMMARY_OUT}"
echo ""

# ============================================================================
# Step 5: Classify sites by quality
# Classification criteria:
#   - Good: Total comparisons >= 16 AND Discordance rate <= 5%
#   - Bad:  Total comparisons >= 16 AND Discordance rate > 5%
#   - Low_support: Total comparisons < 16
# Input: ${OUTDIR}/site_discordance_summary.tsv
# Output: ${OUTDIR}/site_quality_classification.tsv
# ============================================================================
echo "=========================================="
echo "Step 5: Site quality classification"
echo "=========================================="

CLASS_OUT="${OUTDIR}/site_quality_classification.tsv"

awk -v OFS="\t" '
NR==1 {
    print $0, "SNP_Class"
    next
}

{
    total = $7
    rate  = $8
    
    if (total >= 16) {
        if (rate <= 0.05)
            class = "Good"
        else
            class = "Bad"
    } else {
        class = "Low_support"
    }
    
    print $0, class
}
' "${SUMMARY_OUT}" > "${CLASS_OUT}"

echo "Step 5 completed! Output file: ${CLASS_OUT}"
echo ""

# ============================================================================
# Step 6: Extract flanking sequences for high-quality SNPs (requires manual filtering first)
# Note: This step requires extracting "Good" sites from site_quality_classification.tsv
#       Save as concordance-analysis/good_snps.tsv (format: CHROM POS REF ALT)
# ============================================================================
echo "=========================================="
echo "Step 6: Extracting flanking sequences for high-quality SNPs"
echo "=========================================="

# Check if good_snps.tsv exists
if [ -f "${OUTDIR}/good_snps.tsv" ]; then
    # Generate 9bp region BED file (5bp upstream and 4bp downstream of SNP, total 9bp)
    awk '{
        start = $2 - 5
        if (start < 0) start = 0
        end   = $2 + 4
        print $1"\t"start"\t"end
    }' "${OUTDIR}/good_snps.tsv" > "${OUTDIR}/good_snps_9bp.bed"
    
    # Extract sequences using bedtools
    if command -v bedtools &> /dev/null; then
        bedtools getfasta \
            -fi "${REF_FASTA}" \
            -bed "${OUTDIR}/good_snps_9bp.bed" \
            -fo "${OUTDIR}/good_snps_9bp.fa"
        echo "Step 6 completed! Output file: ${OUTDIR}/good_snps_9bp.fa"
    else
        echo "Warning: bedtools not installed, skipping sequence extraction"
    fi
else
    echo "Skipping Step 6: ${OUTDIR}/good_snps.tsv not found"
    echo "Hint: Please manually extract Good sites from ${CLASS_OUT} to generate good_snps.tsv"
    echo "      Example command: awk '\$9==\"Good\" {print \$1,\$2,\$3,\$4}' ${CLASS_OUT} > ${OUTDIR}/good_snps.tsv"
fi
echo ""

# ============================================================================
# Step 7: Visualization analysis (embedded Python script)
# Input: ${OUTDIR}/site_level_concordance_matrix.tsv
# Output: ${OUTDIR}/snp_sharing_concordance.pdf and ${OUTDIR}/snp_stats.tsv
# ============================================================================
echo "=========================================="
echo "Step 7: Visualization analysis (Python)"
echo "=========================================="

# Create temporary Python script
PYTHON_SCRIPT=$(mktemp)

cat > "${PYTHON_SCRIPT}" << 'EOF'
#!/usr/bin/env python3
"""
SNP Matrix Analysis Script
==========================
Count SNPs appearing in N sample pairs and calculate genotype concordance
Input: site_level_concordance_matrix.tsv (wide table, columns=sample pairs, values=0/1/NA)
Output: Bar plot snp_sharing_concordance.pdf + statistics table snp_stats.tsv
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
import os

# ─────────────────────────────────────────────
# 1. Read file
# ─────────────────────────────────────────────
input_file = "concordance-analysis/site_level_concordance_matrix.tsv"

if not os.path.exists(input_file):
    print(f"[ERROR] File not found: {input_file}")
    print("Please ensure Step 2 has been run to generate this file")
    sys.exit(1)

print(f"[INFO] Reading file: {input_file}")
df = pd.read_csv(input_file, sep="\t", low_memory=False)
print(f"[INFO] Raw data dimensions: {df.shape[0]} rows × {df.shape[1]} columns")

# ─────────────────────────────────────────────
# 2. Identify metadata columns vs sample pair columns
#    Sample pair column pattern: contains "_vs_" (e.g., 3DWF22T0036_vs_sample15)
# ─────────────────────────────────────────────
meta_cols = ["CHROM", "POS", "REF", "ALT"]
sample_cols = [c for c in df.columns if c not in meta_cols]
n_samples = len(sample_cols)
print(f"[INFO] Detected sample pair columns: {n_samples}")
print(f"[INFO] Sample pair column examples: {sample_cols[:3]} ...")

# ─────────────────────────────────────────────
# 3. Convert sample pair columns to numeric (0=concordant, 1=discordant, NA=missing)
# ─────────────────────────────────────────────
sample_df = df[sample_cols].apply(pd.to_numeric, errors="coerce")

# presence matrix: value 0 indicates "concordant in this sample pair"
# Note: "concordant" means identical genotypes (Discordant=0)
presence = (~sample_df.isna()) & (sample_df == 0)

# ─────────────────────────────────────────────
# 4. Calculate how many sample pairs each SNP appears in (concordantly)
# ─────────────────────────────────────────────
snp_count = presence.sum(axis=1)   # Number of sample pairs with concordance for this row
df["n_samples_present"] = snp_count

# ─────────────────────────────────────────────
# 5. Calculate genotype concordance
#    Definition: Among sample pairs where this SNP appears, proportion of concordant calls
#    i.e., proportion of non-NA values that are 0 (concordant)
# ─────────────────────────────────────────────
def calc_concordance(row):
    vals = row.dropna()
    if len(vals) == 0:
        return np.nan
    # Concordance = proportion of values that are 0 (concordant)
    mode_count = (vals == 0).sum()
    return mode_count / len(vals)

print("[INFO] Calculating concordance...")
df["concordance"] = sample_df.apply(calc_concordance, axis=1)

# ─────────────────────────────────────────────
# 6. Group statistics by "number of sample pairs present"
# ─────────────────────────────────────────────
all_counts = range(1, n_samples + 1)

stats = []
for k in all_counts:
    subset = df[df["n_samples_present"] == k]
    n_snps = len(subset)
    mean_conc = subset["concordance"].mean()
    stats.append({
        "n_sample_pairs": k,
        "n_snps": n_snps,
        "mean_concordance": round(mean_conc, 4) if not np.isnan(mean_conc) else np.nan
    })

stats_df = pd.DataFrame(stats)

# Save statistics table to output directory
output_stats = "concordance-analysis/snp_stats.tsv"
stats_df.to_csv(output_stats, sep="\t", index=False)
print(f"[INFO] Statistics table saved: {output_stats}")
print(stats_df.to_string(index=False))

# ─────────────────────────────────────────────
# 7. Generate dual-axis bar plot
# ─────────────────────────────────────────────
fig, ax1 = plt.subplots(figsize=(14, 6))

x = stats_df["n_sample_pairs"].values
y_snps = stats_df["n_snps"].values
y_conc = stats_df["mean_concordance"].values

# ── Primary axis: SNP count bar plot
colors = plt.cm.Blues(np.linspace(0.35, 0.85, len(x)))
bars = ax1.bar(x, y_snps, color=colors, edgecolor="white",
               linewidth=0.6, zorder=2, label="SNP count")

ax1.set_xlabel("Number of Sample Pairs with Concordant Genotype", fontsize=13)
ax1.set_ylabel("Number of SNPs", fontsize=13, color="#2166ac")
ax1.tick_params(axis="y", labelcolor="#2166ac")
ax1.set_xticks(x)
ax1.set_xticklabels([str(i) for i in x], fontsize=10)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
ax1.set_xlim(0.3, n_samples + 0.7)
ax1.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax1.set_axisbelow(True)

# Add numbers on top of bars (only show if value > 0)
if max(y_snps) > 0:
    for bar, val in zip(bars, y_snps):
        if val > 0:
            ax1.text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + max(y_snps) * 0.01,
                     f"{val:,}", ha="center", va="bottom", fontsize=7.5,
                     color="#2166ac", fontweight="bold")

# ── Secondary axis: Concordance line plot
ax2 = ax1.twinx()
valid = ~np.isnan(y_conc)
if any(valid):
    ax2.plot(x[valid], y_conc[valid] * 100, color="#d6604d",
             marker="o", markersize=5, linewidth=2,
             label="Mean concordance", zorder=3)
ax2.set_ylabel("Mean Concordance (%)", fontsize=13, color="#d6604d")
ax2.tick_params(axis="y", labelcolor="#d6604d")
ax2.set_ylim(0, 110)
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))

# ── Combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2,
           loc="upper right", fontsize=10, framealpha=0.85)

plt.title("SNP Sharing Across Sample Pairs and Genotype Concordance",
          fontsize=14, fontweight="bold", pad=12)
plt.tight_layout()

output_figure = "concordance-analysis/snp_sharing_concordance.pdf"
plt.savefig(output_figure, bbox_inches="tight")
plt.close()
print(f"[INFO] Figure saved: {output_figure}")
print("[DONE] Analysis complete!")
EOF

# Run Python script
python3 "${PYTHON_SCRIPT}"

# Delete temporary Python script
rm -f "${PYTHON_SCRIPT}"

echo "Step 7 completed!"
echo ""

# ============================================================================
# Clean up temporary files (optional)
# ============================================================================
echo "=========================================="
echo "Cleaning up temporary files"
echo "=========================================="

read -p "Delete temporary file directory ${TMPDIR}? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    rm -rf "${TMPDIR}"
    echo "Temporary files deleted"
else
    echo "Keeping temporary files: ${TMPDIR}"
fi

echo ""
echo "=========================================="
echo "All analyses completed!"
echo "=========================================="
echo "Output files:"
echo "  1. ${OUTDIR}/all_sites_concordance.txt          - Detailed comparison for all sites"
echo "  2. ${OUTDIR}/site_level_concordance_matrix.tsv  - Concordance matrix (0/1)"
echo "  3. ${OUTDIR}/site_level_concordance_detailed.tsv- Detailed matrix"
echo "  4. ${OUTDIR}/site_discordance_summary.tsv       - Site discordance rate statistics"
echo "  5. ${OUTDIR}/site_quality_classification.tsv    - Site quality classification"
echo "  6. ${OUTDIR}/good_snps_9bp.fa                   - Flanking sequences for Good SNPs (requires manual filtering)"
echo "  7. ${OUTDIR}/snp_sharing_concordance.pdf        - Visualization plot"
echo "  8. ${OUTDIR}/snp_stats.tsv                      - Statistics table"
echo "=========================================="
```

**This step re-classifies each SNP by the number of sample pairs in which it had a genotype comparison (concordant or discordant), rather than only counting concordant comparisons as "present," which had previously caused mean concordance to be artificially inflated near 100%. Genotype concordance is now calculated independently for each site, so the resulting bar chart and concordance line reflect true, unbiased patterns of SNP sharing and agreement across sample pairs.

```
#!/usr/bin/env python3
"""
SNP Matrix Analysis Script (Fixed Version)
===========================================
For each SNP site, count how many sample pairs it was COMPARED IN
(regardless of whether the comparison was concordant or discordant),
and separately compute genotype concordance for that site.

Fix from the previous version:
Previously, "presence" was defined as value == 0, meaning only
concordant comparisons were counted as "present". Since sites were
then grouped by that same count, every group ended up almost entirely
concordant by construction, making mean concordance ~100% regardless
of the underlying data (a circular result).

Here, "presence" is defined as value is not NA, i.e. the site had a
comparison result in that sample pair (concordant=0 or discordant=1).
Concordance is computed independently of this presence count.

Input : site_level_concordance_matrix.tsv
        (wide table; rows = SNP sites, columns = sample pairs,
         values = 0/1/NA)
Output: snp_stats_fixed.tsv                (summary table)
        snp_sharing_concordance_fixed.pdf  (bar + line chart)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
import os

# -------------------------------------------------
# 1. Load input matrix
# -------------------------------------------------
input_file = sys.argv[1] if len(sys.argv) > 1 else "concordance-analysis/site_level_concordance_matrix.tsv"

if not os.path.exists(input_file):
    print(f"[ERROR] File not found: {input_file}")
    print("Usage: python analyze_snp_matrix_fixed.py site_level_concordance_matrix.tsv")
    sys.exit(1)

output_dir = os.path.dirname(os.path.abspath(input_file))
output_stats = os.path.join(output_dir, "snp_stats_fixed.tsv")
output_figure = os.path.join(output_dir, "snp_sharing_concordance_fixed.pdf")

df = pd.read_csv(input_file, sep="\t", low_memory=False)
total_sites = df.shape[0]
print(f"[INFO] Loaded {total_sites:,} SNP sites, {df.shape[1]} columns")

# -------------------------------------------------
# 2. Split metadata columns from sample-pair columns
# -------------------------------------------------
meta_cols = ["CHROM", "POS", "REF", "ALT"]
sample_cols = [c for c in df.columns if c not in meta_cols]
n_pairs = len(sample_cols)

sample_df = df[sample_cols].apply(pd.to_numeric, errors="coerce")

# -------------------------------------------------
# 3. Presence = compared in this sample pair (value is 0 or 1, not NA)
# -------------------------------------------------
presence = ~sample_df.isna()
df["n_samples_present"] = presence.sum(axis=1)

# -------------------------------------------------
# 4. Concordance = fraction of concordant (0) calls among comparisons
# -------------------------------------------------
def calc_concordance(row):
    vals = row.dropna()
    if len(vals) == 0:
        return np.nan
    return (vals == 0).sum() / len(vals)

df["concordance"] = sample_df.apply(calc_concordance, axis=1)

# -------------------------------------------------
# 5. Group sites by number of sample pairs they were compared in
# -------------------------------------------------
stats = []
for k in range(1, n_pairs + 1):
    subset = df[df["n_samples_present"] == k]
    mean_conc = subset["concordance"].mean()
    stats.append({
        "n_sample_pairs": k,
        "n_snps": len(subset),
        "mean_concordance": round(mean_conc, 4) if not np.isnan(mean_conc) else np.nan
    })
stats_df = pd.DataFrame(stats)

# Sanity check: total SNPs across bars should equal total rows in the matrix
bar_total = stats_df["n_snps"].sum()
print(f"[CHECK] Sum across bars = {bar_total:,} | Matrix rows = {total_sites:,} | "
      f"{'OK' if bar_total == total_sites else 'MISMATCH'}")

stats_df.to_csv(output_stats, sep="\t", index=False)
print(f"[INFO] Summary table saved: {output_stats}")

# -------------------------------------------------
# 6. Plot: SNP count (bars) + mean concordance (line)
# -------------------------------------------------
fig, ax1 = plt.subplots(figsize=(14, 6))
x = stats_df["n_sample_pairs"].values
y_snps = stats_df["n_snps"].values
y_conc = stats_df["mean_concordance"].values

colors = plt.cm.Blues(np.linspace(0.35, 0.85, len(x)))
bars = ax1.bar(x, y_snps, color=colors, edgecolor="white", linewidth=0.6, zorder=2, label="SNP count")
ax1.set_xlabel("Number of Sample Pairs Where SNP Was Compared", fontsize=13)
ax1.set_ylabel("Number of SNPs", fontsize=13, color="#2166ac")
ax1.tick_params(axis="y", labelcolor="#2166ac")
ax1.set_xticks(x)
ax1.set_xticklabels([str(i) for i in x], fontsize=10)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
ax1.set_xlim(0.3, n_pairs + 0.7)
ax1.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax1.set_axisbelow(True)

for bar, val in zip(bars, y_snps):
    if val > 0:
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(y_snps) * 0.01,
                  f"{val:,}", ha="center", va="bottom", fontsize=7.5, color="#2166ac", fontweight="bold")

ax2 = ax1.twinx()
valid = ~np.isnan(y_conc)
ax2.plot(x[valid], y_conc[valid] * 100, color="#d6604d", marker="o", markersize=5,
          linewidth=2, label="Mean concordance", zorder=3)
ax2.set_ylabel("Mean Concordance (%)", fontsize=13, color="#d6604d")
ax2.tick_params(axis="y", labelcolor="#d6604d")
ax2.set_ylim(0, 110)
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=10, framealpha=0.85)

plt.title("SNP Sharing Across Sample Pairs and Genotype Concordance (Fixed)",
          fontsize=14, fontweight="bold", pad=12)
plt.tight_layout()
plt.savefig(output_figure, bbox_inches="tight")
plt.close()
print(f"[INFO] Figure saved: {output_figure}")
```

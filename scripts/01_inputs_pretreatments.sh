#!/bin/bash

########################################################################

set -euo pipefail
trap 'echo "ERROR : Pipeline failed at line $LINENO. Exiting." >&2' ERR

echo -e "\n############ FASTA INDEXING ############\n"

if [ -z "$REF_FASTA" ]; then
    echo "ERROR : NO .fasta or .fa file provided."
    exit 1
fi

echo "Reference genome (FASTA) : $REF_FASTA"
echo "Indexing ..."
samtools faidx "$REF_FASTA" # -o "${OUTPUT_DIR}/$(basename "${REF_FASTA}").fai"

########################################################################

echo -e "\n############ VCF PRETREATMETS ############\n"

if [ -z "$VCF" ]; then
    echo "ERROR : NO .vcf or .vcf.gz file provided."
    exit 1
fi


# zip if uncompressed
if [[ "$VCF" == *.vcf ]]; then
    echo "Compressing VCF to $VCF_COMPRESSED ..."
    # check if already exists
    if [ ! -f "$VCF_COMPRESSED" ]; then
        bgzip -c "$VCF" > "$VCF_COMPRESSED"
    else
        echo "INFO : Compressed VCF already exists. Skipping bgzip."
    fi
fi




echo "Indexing VCF ..."
bcftools index "$VCF_COMPRESSED" -o "${OUTPUT_DIR}/$(basename "${VCF_COMPRESSED}").csi"

echo "Filtering VCF (QUAL>30, MAF>0.01, DP>10, F_MISSING<=$MISSING_FILTER ) ..."

# TD add -Ou for intermediate steps ? faster ?
# add depth if needed after normalisation : bcftools filter -Ou -S . -e 'FMT/DP<10'
# normalisation is needed for bcftools to correctly ignore all indels and correctly count the remaining filters 
bcftools view -e 'ALT~"\*"' -Ou "$VCF_COMPRESSED" \
  | bcftools norm -m -any \
  | bcftools view -v snps -i "QUAL>30 && MAF>0.01 && F_MISSING<=$MISSING_FILTER" -Oz -o "$VCF_FILTERED" \
  || { echo "ERROR : bcftools failed."; exit 1; }

# index the filtered VCF
echo "Indexing filtered VCF ..."
bcftools index "$VCF_FILTERED" -o "${VCF_FILTERED}.csi"

########################################################################

echo "############ GFF SORTING ############"

if [ -z "$GFF" ]; then
    echo "ERROR : NO .gff file in provided."
    exit 1
fi

echo "Found GFF : $GFF"

echo "GFF sorting ..."
awk -F'\t' '{
    priority = 3
    if ($3 == "mRNA") priority = 1
    else if ($3 == "CDS") priority = 2
    else if ($3 == "UTR") priority = 2
    # temp cols for sorting:
    # contig , start pos , end pos , priority , original GFF(line)
    printf "%s\t%d\t%d\t%d\t%s\n", $1, $4, $5, priority, $0
}' "$GFF" | sort -k1,1V -k2,2n -k3,3nr -k4,4n | cut -f5- > "$GFF_SORTED"
    echo "Sorted GFF in : $GFF_SORTED"


echo "FILTERING & SORTING DONE."
echo ""
echo "########################################################################"
echo ""
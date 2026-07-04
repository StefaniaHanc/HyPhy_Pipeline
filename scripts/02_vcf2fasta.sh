#!/bin/bash

set -euo pipefail
trap 'echo "ERROR : Pipeline failed at line $LINENO. Exiting." >&2' ERR

FASTA_FINAL_DIR="${OUTPUT_DIR}/FASTA_BY_GENE_CDS"

if [ -d "$FASTA_FINAL_DIR" ]; then
    echo "WARNING : Output dir $FASTA_FINAL_DIR exists. Overwriting ..."
    rm -rf "$FASTA_FINAL_DIR"
fi

if [ ! -f "$SCRIPT_DIR/scripts/vcf2fasta/vcf2fasta.py" ]; then
    echo "ERROR : vcf2fasta.py not found."
    exit 1
fi


########################################################################

# VCF2FASTA software was adapted to correctly output IUPAC ambiguitty codes, manage minus strands and indexing corrections
echo "Extracting gene sequences ... (may take a while)"

# command with options that are needed for HyPhy run
# this will project the SNPs on the reference genome and extract genes based on the annotation file
# -blend : # concatenates GFF entries of FEAT (CDS) into a single alignment
# --inframe : # forces the first codon of the sequence to be inframe
python3 "$SCRIPT_DIR/scripts/vcf2fasta/vcf2fasta.py" \
    --fasta "$REF_FASTA" \
    --vcf "$VCF_FILTERED" \
    --gff "$GFF_SORTED" \
    --feat CDS \
    --blend \  
    --inframe \  
    --out "$FASTA_FINAL_DIR"

echo "FASTA files of genes in : $FASTA_FINAL_DIR"

########################################################################

echo "########### SEPARATING MONOMORPHIC GENES ... "###########"
MONOMORPH_DIR="${OUTPUT_DIR}/MONOMORPH_GENES"
mkdir -p "$MONOMORPH_DIR"

SUMMARY_FILE="${FASTA_FINAL_DIR}/variant_counts_summary.tsv"

if [ -f "$SUMMARY_FILE" ]; then
    echo "Moving monomorphic genes ..."
    # skip header
    tail -n +2 "$SUMMARY_FILE" | while IFS=$'\t' read -r GENE COUNT; do
        if [ "$COUNT" -eq 0 ]; then
            if [ -f "${FASTA_FINAL_DIR}/${GENE}.fa" ]; then
                mv "${FASTA_FINAL_DIR}/${GENE}.fa" "$MONOMORPH_DIR/"
            fi
        fi
    done
    
    MOVED_COUNT=$(find "$MONOMORPH_DIR" -maxdepth 1 -name "*.fa" | wc -l)
    echo "INFO : contains $MOVED_COUNT monomorphic genes."

else
    echo "ERROR : $SUMMARY_FILE not found."
    exit 1
fi

echo "VCF2FASTA DONE."

echo -e "\n########################################################################\n"
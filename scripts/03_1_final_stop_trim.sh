#!/bin/bash

set -euo pipefail
trap 'echo "ERROR : Pipeline failed at line $LINENO. Exiting." >&2' ERR

FASTA_FINAL_DIR=${FASTA_FINAL_DIR:-"."} 
GENE=${GENE:-""} 
GLOBAL_OUT=${OUTPUT_DIR:-"."}

# global output folder
HYPHY_OUT_DIR="${GLOBAL_OUT}/SEQS_FOR_HYPHY"
SUMMARY_FILE="${GLOBAL_OUT}/final_stops_variants_info.tsv"

mkdir -p "$HYPHY_OUT_DIR"

echo -e "Gene_Name\tCodon\tNb_of_indivs" > "$SUMMARY_FILE"
echo "Processing final stop codons ..."

# handle cases with no matches
shopt -s nullglob

for file_path in "${FASTA_FINAL_DIR}/${GENE}"*.fa "${FASTA_FINAL_DIR}/${GENE}"*.fasta; do
    filename=$(basename "$file_path")
    gene_name="${filename%.*}"
    out_file="$HYPHY_OUT_DIR/${gene_name}.fa"
    
    # process final stops
    awk -v gene="$gene_name" -v outfile="$out_file" -v report="$SUMMARY_FILE" '
        BEGIN {
            # allowed stop codons including IUPAC that only translates into STOP
            allowed["TAA"] = 1;
            allowed["TGA"] = 1;
            allowed["TAG"] = 1;
            allowed["TRA"] = 1;
            allowed["TAR"] = 1;
        }
        
        # FASTA headers
        /^>/ {
            if (seq != "") process_seq();
            header = $0;
            seq = "";
            next;
        }
        
        # sequence lines, remove enter
        {
            gsub(/\r/, "");
            seq = seq $0;
        }
        
        # final sequence; write reports
        END {
            if (seq != "") process_seq();
            
            # non-standard codons (STOP codon non-synonymous mutations) => summary file
            for (codon in non_standard) {
                print gene "\t" codon "\t" non_standard[codon] >> report;
            }
        }
        


        # equalize lengths of all sequences


        function process_seq() {
            len = length(seq);
            if (len >= 3) {
                # last 3 bases
                last3 = toupper(substr(seq, len-2, 3));
                if (!(last3 in allowed)) {
                    non_standard[last3]++;
                }
                trimmed = substr(seq, 1, len-3);
                print header > outfile;
                
                # 80-char blocks
                trimmed_len = length(trimmed);
                if (trimmed_len == 0) {
                    print "" > outfile;
                } else {
                    for (i = 1; i <= trimmed_len; i += 80) {
                        print substr(trimmed, i, 80) > outfile;
                    }
                }
            } else {
                print header > outfile;
                print seq > outfile;
            }
        }
    ' "$file_path"
done

echo ""
echo "DONE. Terminal-codon trimmed sequences in : $HYPHY_OUT_DIR"
echo "Report in $SUMMARY_FILE"

echo ""
echo "########################################################################"
echo ""


echo -e "Moving FASTAS to separate folders. \n This step is essential for HyPhy treatment of multiple genes in parallel.)"
for file in "$HYPHY_OUT_DIR"/*.fa; do
    genename="${file%.fa}"
    mkdir -p "$genename"
    mv "$file" "$genename/"
done

shopt -u nullglob
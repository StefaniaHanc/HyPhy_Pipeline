#!/bin/bash

set -euo pipefail
trap 'echo "ERROR : Pipeline failed at line $LINENO. Exiting." >&2' ERR

INPUT_DIR="${PRUNING_INPUT_DIR:-${FASTA_FINAL_DIR:-.}}"                           
TREE="${TREE:-.}"    
MAX_DROP="${MAX_DROP:-30}"                                    

# outputs
GLOBAL_OUT="${OUTPUT_DIR:-.}"
BASE_OUT="${GLOBAL_OUT}/PROCESSED"
DIR_PRUNED="$DIR_WITH_STOPS/PRUNED_SEQS_FOR_HYPHY"
DIR_TRIMMED="$DIR_WITH_STOPS/TRIMMED_SEQS_FOR_HYPHY"
DIR_CLEAN="$DIR_WITH_STOPS/CLEAN_SEQS" 

REPORT_TSV="$BASE_OUT/internal_stops_report.tsv"
TMP_STOPS="$BASE_OUT/.global_stops.tmp"

mkdir -p "$DIR_PRUNED" "$DIR_TRIMMED" "$DIR_CLEAN"
rm -f "$REPORT_TSV" "$TMP_STOPS"

# headers
echo -e "Gene_Name\tAction\tDropped_Count\t%_Gene_Kept\tTotal_Stops_In_Gene\tFirst_Stop_Pos\tDropped_Individuals_List" > "$REPORT_TSV"
echo -e "\n############################################\n"
echo "=============== Internal STOP codon processing ... =================="
echo "Pruning max $MAX_DROP individuals with stops at the 1st stop location."

########################################################################

find "$INPUT_DIR" -maxdepth 2 -type f \( -name "*.fa" -o -name "*.fasta" \) | while read -r fasta_file; do
    filename=$(basename "$fasta_file")
    gene="${filename%.*}"

    # prevent infinite loop if output dir is inside input dir
    if [[ "$fasta_file" == *"$BASE_OUT"* ]]; then
        continue
    fi

    # folder setup
    mkdir -p "$DIR_PRUNED/$gene" "$DIR_TRIMMED/$gene" "$DIR_CLEAN/$gene"


    #passing shell variables
    awk -v gene="$gene" \
        -v dir_pruned="$DIR_PRUNED" \
        -v dir_trimmed="$DIR_TRIMMED" \
        -v dir_clean="$DIR_CLEAN" \
        -v report="$REPORT_TSV" \
        -v max_drop="$MAX_DROP" \
        -v tmp_stops="$TMP_STOPS" '
    BEGIN { FS=" "; OFS="\t" }
    


    #######################################################

    # parse fasta
    /^>/ {
        # save to seqs array with header as key
        if (seq != "") seqs[header] = seq
        # store new header, reset seq var
        header = $0
        seq = ""
        next
    }
    {
        gsub(/[ \t\r\n]+/, "", $0)
        seq = seq toupper($0)
    }
    # after the file has been read
    END {
        # save last seq to the array
        if (seq != "") seqs[header] = seq
        # array of stop codons
        allowed["TAA"]=1; allowed["TAG"]=1; allowed["TGA"]=1; allowed["TAR"]=1; allowed["TRA"]=1;

        # init tracking
        global_first_stop = -1
        total_stops = 0
        max_len = 0
        
        # loop to find stops
        for (h in seqs) {
            s = seqs[h]
            len = length(s)
            # max seq length
            if (len > max_len) max_len = len
            first_stop[h] = -1           
            for (i = 1; i <= len - 2; i += 3) {
                codon = substr(s, i, 3)
                # increment if current codon in array
                if (codon in allowed) {
                    total_stops++
                    # record position
                    if (first_stop[h] == -1) first_stop[h] = i
                    # update global earliest stop codon position
                    if (global_first_stop == -1 || i < global_first_stop) global_first_stop = i
                }
            }
        }

        # track how many seqs share a stop codon at the absolute earliest position
        count_at_first = 0
        if (global_first_stop != -1) {
            for (h in seqs) {
                if (first_stop[h] == global_first_stop) {
                    count_at_first++
                    # mark to be dropped
                    drop_this_ind[h] = 1
                    dropped_taxa[count_at_first] = h
                    # get codons causing issue
                    codon_str = substr(seqs[h], global_first_stop, 3)
                    first_stop_type[codon_str]++
                }
            }
        }

        # default outputs
        action = "Clean"
        cut_length = max_len
        gene_kept_pct = 100
        dropped_str = "-"
        dropped_raw = ""
        num_dropped = 0

        # handle premature stops
        if (global_first_stop != -1) {
            # if < than threshold
            if (count_at_first <= max_drop) {
                num_dropped = count_at_first
                # format names
                for (i = 1; i <= num_dropped; i++) {
                    raw_name = dropped_taxa[i]
                    gsub(/^>/, "", raw_name) # strip >
                    dropped_str = (dropped_str == "-") ? raw_name : dropped_str ", " raw_name
                    dropped_raw = (dropped_raw == "") ? raw_name : dropped_raw "\n" raw_name
                }
                
                # check for + stops
                next_first_stop = -1
                for (h in seqs) {
                    if (h in drop_this_ind) continue  # skip taxa if already dropping
                    if (first_stop[h] != -1) {
                        if (next_first_stop == -1 || first_stop[h] < next_first_stop) {
                            next_first_stop = first_stop[h]
                        }
                    }
                }
                
                # if remaining seqs still have stops, prune the taxa & truncate rest of seq
                if (next_first_stop != -1) {
                    action = "Pruned+Shortened"
                    cut_length = next_first_stop - 1
                } else {
                    # oftherwise just drop
                    action = "Pruned"
                }

                # out dirs
                out_file = dir_pruned "/" gene "/" gene ".fa"
                print dropped_raw > (dir_pruned "/" gene "/" gene ".to_drop")
                
            } else {
                # if threshold exceeded, trim
                action = "Shortened_All"
                cut_length = global_first_stop - 1  # cut before stop
                out_file = dir_trimmed "/" gene "/" gene ".fa"
                dropped_str = "ALL_CUT_TO_POS_" cut_length
            }
            # % of seq remaining
            gene_kept_pct = (cut_length / max_len) * 100
        } else {
            # if no stop codons 
            action = "Clean"
            out_file = dir_clean "/" gene "/" gene ".fa"
        }

        # write processed seqs
        for (h in seqs) {
            if ((action == "Pruned" || action == "Pruned+Shortened") && (h in drop_this_ind)) continue
            clean_seq = substr(seqs[h], 1, cut_length)
            
            print h > out_file
            # cut to 80 chars outputs
            for (i = 1; i <= length(clean_seq); i += 80) {
                print substr(clean_seq, i, 80) > out_file
            }
        }
        # summary stats for gene
        first_pos_str = (global_first_stop != -1) ? global_first_stop : "-"
        print gene, action, num_dropped, sprintf("%.2f", gene_kept_pct), total_stops, first_pos_str, dropped_str >> report

        # counts of specific stop codon types to tmp
        if (global_first_stop != -1 && count_at_first > 0) {
            print (first_stop_type["TAA"]+0) "\t" (first_stop_type["TAG"]+0) "\t" (first_stop_type["TGA"]+0) "\t" (first_stop_type["TAR"]+0) "\t" (first_stop_type["TRA"]+0) >> tmp_stops
        }
    }
    ' "$fasta_file"

    # 3. BASH CLEANUP: Remove the empty directories awk didn't need
    find "$DIR_PRUNED/$gene" "$DIR_TRIMMED/$gene" "$DIR_CLEAN/$gene" -empty -type d -delete

########################################################################

    # TREE PRUNING 
    DROP_FILE="$DIR_PRUNED/${gene}/${gene}.to_drop"
    FASTA_FILE="$DIR_PRUNED/${gene}/${gene}.fa"
    OUT_TREE="$DIR_PRUNED/${gene}/${gene}.nwk"

    if [ -f "$DROP_FILE" ]; then
        if [ -f "$TREE" ]; then
            # Securely pass variables to Python
            export PY_TREE="$TREE"
            export PY_FASTA="$FASTA_FILE"
            export PY_OUT="$OUT_TREE"
            export PY_GENE="$gene"
            
            python3 -c "
            
import sys, os
try:
    from Bio import Phylo
except ImportError:
    print(f'WARNING : biopython required to prune tree. {os.environ.get(\"PY_GENE\")} tree skipped', file=sys.stderr)
    sys.exit(1)

# get file paths
try:
    tree_path = os.environ.get('PY_TREE')
    fasta_path = os.environ.get('PY_FASTA')
    out_tree = os.environ.get('PY_OUT')
    gene_name = os.environ.get('PY_GENE')

    tree = Phylo.read(tree_path, 'newick')
    
    # empty set to store names of the taxa
    fasta_taxa = set()
    #extract sequence names
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                fasta_taxa.add(line.strip()[1:])

    # list of the names of all terminal nodes            
    tree_taxa = [tip.name for tip in tree.get_terminals() if tip.name is not None]
    
    for tip in tree_taxa:
        if tip not in fasta_taxa:
            try:
                #remove the tip
                tree.prune(tip)
            except ValueError:
                pass
    # should not be less than 3 
    if len(tree.get_terminals()) >= 3:
        Phylo.write(tree, out_tree, 'newick')
    else:
        print(f'WARNING : tree for {gene_name} pruned to < 3 taxa.', file=sys.stderr)
        # delete corresponding fasta
        if os.path.exists(fasta_path):
            os.remove(fasta_path)

except Exception as e:
    print(f'Error processing tree for {os.environ.get(\"PY_GENE\")}: {e}', file=sys.stderr)
"
            rm -f "$DROP_FILE"
        else
            echo "WARNING : Global tree file $TREE not found - $gene." >&2
        fi
    fi

done

########################################################################

# SUMMARY
echo "First global STOP codon summary : "
if [ -f "$TMP_STOPS" ]; then
    awk '{taa+=$1; tag+=$2; tga+=$3; tar+=$4; tra+=$5} END {
        print " TAA: " taa+0 " | TAG: " tag+0 " | TGA: " tga+0 " | TAR: " tar+0 " | TRA: " tra+0
    }' "$TMP_STOPS"
    rm -f "$TMP_STOPS"
else
    echo " No internal stop codons found across any genes."
fi

echo -e "\n############################################\n"

echo "DONE internal STOP codon processing."
echo "TSV Report : $REPORT_TSV"
echo "Pruned (≤ 30 stops) : $DIR_PRUNED/"
echo "Shortened (> 30)    : $DIR_TRIMMED/"
echo "Cleaned Fastas      : $DIR_CLEAN/"
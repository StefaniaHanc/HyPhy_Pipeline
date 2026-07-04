#!/bin/bash

set -euo pipefail
trap 'echo "ERROR : Pipeline failed at line $LINENO. Exiting." >&2' ERR

if [ -z "$TREE" ] || [ ! -f "$TREE" ]; then
    echo "ERROR : Global tree not given or file not found ($TREE). Cannot run HyPhy."
    exit 1
fi

PRUNED_DIR="${OUTPUT_DIR}/PROCESSED/PRUNED_SEQS_FOR_HYPHY"
TRIMMED_DIR="${OUTPUT_DIR}/PROCESSED/TRIMMED_SEQS_FOR_HYPHY"
CLEAN_DIR="${OUTPUT_DIR}/PROCESSED/CLEAN_SEQS"

HYPHY_OUT_DIR="${OUTPUT_DIR}/HYPHY_OUTPUTS"
HYPHY_LOG_DIR="${HYPHY_OUT_DIR}/logs"
ERROR_LIST="${HYPHY_LOG_DIR}/failed_genes_HyphyLogError.txt"


mkdir -p "$HYPHY_OUT_DIR" "$HYPHY_LOG_DIR"
true > "$ERROR_LIST" 

########################################################################

# positional parameters -> represent the order of the arguments passed to the function when we call it via parallel
run_hyphy() {
    local file=$1
    local global_tree=$2
    local out_dir=$3
    local log_dir=$4
    local model=$5
    
    local filename=$(basename "$file")
    local gene_name="${filename%.*}"
    local file_dir=$(dirname "$file")
    
    # tree selection
    # .nwk file in the same directory as the FASTA (maxdepth 1)
    local local_tree=$(find "$file_dir" -maxdepth 1 -name "*.nwk" 2>/dev/null | head -n 1)
    local active_tree
    
    if [[ -n "$local_tree" ]]; then
        active_tree="$local_tree"
        echo "Local tree used : $gene_name: $active_tree" >> "$log_dir/${gene_name}.terminal.log"
    else
        active_tree="$global_tree"
        echo "Global tree used : $gene_name: $active_tree" >> "$log_dir/${gene_name}.terminal.log"
    fi

########################################################################
# HyPhy exec

    # CPU=1 limits each gene analysis to 1 core for quicker runs
    # &>> to redirect and append stdout and stderr to the log file

    hyphy CPU=1 ENV=TOLERATE_NUMERICAL_ERRORS=1 "$model" \
        --alignment "$file" \
        --tree "$active_tree" \
        --output "$out_dir/${gene_name}.${model}.json" &>> "$log_dir/${gene_name}.terminal.log"

    # if non-zero exit code
    if [ $? -ne 0 ]; then
        echo "$gene_name" >> "$out_dir/failed_genes.txt"
        echo "Failed : $gene_name" >> "$log_dir/failed_genes_HyphyError.txt"
        return 1
    fi
}

# EXPORT -> propagates the function and variables to subprocesses 
# makes them into environment vars so GNU Parallel can see them
export -f run_hyphy

########################################################################

# GNU PARALLEL

echo "Parallelisation for HyPhy $HYPHY_MODEL using $THREADS cores ..."

# detects all fasta files in 3 directories from the previous step
find "$PRUNED_DIR" "$TRIMMED_DIR" "$CLEAN_DIR" -mindepth 2 -name "*.fa" 2>/dev/null \
    | parallel --env run_hyphy --jobs "$THREADS" --progress run_hyphy {} "$TREE" "$HYPHY_OUT_DIR" "$HYPHY_LOG_DIR" "$HYPHY_MODEL"


echo -e "\nHyPhy run finished."
echo "Results saved in : $HYPHY_OUT_DIR"
echo "Terminal logs in : $HYPHY_LOG_DIR"

if [ -s "$HYPHY_OUT_DIR/failed_genes.txt" ]; then
    echo "WARNING : Some genes failed. Check $HYPHY_OUT_DIR/failed_genes.txt"
    echo "Error log in : $ERROR_LIST"
fi
echo "##############################################################"
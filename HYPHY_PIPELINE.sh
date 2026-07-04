#!/bin/bash

# in case of error, outputs the last non-zero status line of code executed (exits, unasigned var check)
set -euo pipefail
trap 'echo "ERROR : Pipeline failed at line $LINENO. Exiting." >&2' ERR

# description of parameters to be used by the pipeline
usage() {
    echo "Usage : path/HYPHY_PIPELINE.sh --ref <REF_FASTA> --vcf <VCF> --gff <GFF> --tree <TREE> -m <HYPHY_MODEL> [options]"
    echo ""
    echo "Mandatory : "
    echo "  --ref     Path to reference FASTA file."
    echo "  --vcf     Path to VCF file."
    echo "  --gff     Path to GFF annotation file."
    echo "  -m        HyPhy model (BUSTED, aBSREL)."
    echo "  --tree    Path to global phylogenetic Newick tree."
    echo ""
    echo "Optional : "
    echo "  -o        Output directory (Default : ./OUTPUTS)."
    echo "  -t        Number of cores (CPUs) to use (Default: max available)."
    echo "  -f        Maximum missing data allowed per site (0.0-1.0) (Default : 0.2)."
    echo "  -s        Path to species file."
    echo "  --gea     Path to file containing GEA genes (gene per line)."
    echo "  --gmt     Path to gene ontology GMT file (for step 7)"
    echo "  --start-at Step number to start from (Default: 1)."
    echo "  --end-at   Step number to end at (Default: 7)."
    exit 1
}



########################################################################

# making sure the required tools are pre-installed

check_dependencies() {
    local dependencies=("awk" "sort" "samtools" "bcftools" "bgzip" "python3" "Rscript" "hyphy" "parallel")
    local missing=0

    echo "Checking software dependencies ..."
    
    # bash tools
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "ERROR : Required tool '$cmd' is not installed or not in PATH." >&2
            missing=1
        fi
    done

    # Python modules 
    if command -v python3 &> /dev/null; then
        for mod in Bio pandas numpy matplotlib pysam art; do
            if ! python3 -c "import $mod" &> /dev/null; then  
                echo "ERROR : Python module '$mod' is missing. (pip install $mod)" >&2
                missing=1
            fi
        done
    fi

    # if command -v Rscript &> /dev/null; then
    #     for pkg in clusterProfiler enrichplot GO.db tidyverse; do
    #         if ! Rscript -e "library('$pkg')" &> /dev/null; then
    #             echo "ERROR : R package '$pkg' is missing. Please install it." >&2
    #             missing=1
    #         fi
    #     done
    # fi

    if [ "$missing" -eq 1 ]; then
        echo "======= EXITING due to missing dependencies. =======" >&2
        exit 1
    fi
    echo "All dependencies found."
}


########################################################################

# setting up vars
START_STEP=1
END_STEP=7
REF_FASTA=""
VCF=""
GFF=""
GMT_FILE=""
HYPHY_MODEL=""
TREE=""
# default is max available cores detected (you may want to use less) ; works for linux & iOS
THREADS=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)  
MISSING_FILTER=0.20 # by default allow ups to 20% missing data per site
GEA_FILE="no"
SPECIES_FILE=""
OUT_BASE="$PWD"


# user's inputs
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --ref) REF_FASTA="$2"; shift ;;
        --vcf) VCF="$2"; shift ;;
        --gff) GFF="$2"; shift ;;
        -m|--model) HYPHY_MODEL="$2"; shift ;;
        -o|--out-dir) OUT_BASE="$2"; shift ;;
        --tree) TREE="$2"; shift ;;
        --gmt) GMT_FILE="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;; 
        -f|--missing) MISSING_FILTER="$2"; shift ;;
        --gea) GEA_FILE="$2"; shift ;;
        -s|--species) SPECIES_FILE="$2"; shift ;;
        --start-at) START_STEP="$2"; shift ;;
        --end-at) END_STEP="$2"; shift ;;
        -h|--help) usage ;;
        *) 
            echo "Unknown parameter : $1"
            usage
            ;;
    esac
    shift
done


# calling the main function to check dependencies
check_dependencies

# checking mandatory parameters
if [ -z "$REF_FASTA" ] || [ -z "$VCF" ] || [ -z "$GFF" ] || [ -z "$HYPHY_MODEL" ] || [ -z "$TREE" ] ; then
    echo "ERROR : Missing mandatory parameters." >&2
    usage
fi


# verif if mandatory files provided
for file in "$REF_FASTA" "$VCF" "$GFF" "$TREE"; do
    if [ ! -f "$file" ]; then
        echo "ERROR : Required file not found: $file" >&2
        exit 1
    fi
done

########################################################################

# exporting the provided variables to be processed across scripts via 'export' ...
OUT_BASE=$(realpath "$OUT_BASE")
export OUTPUT_DIR="${OUT_BASE}/OUTPUTS"
mkdir -p "$OUTPUT_DIR"

# absolute path of the dir where the script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
export SCRIPT_DIR
export REF_FASTA=$(realpath "$REF_FASTA")
export VCF=$(realpath "$VCF")
export GFF=$(realpath "$GFF")
export TREE=$(realpath "$TREE")
export TREE
export MISSING_FILTER
export THREADS

HYPHY_MODEL=$(echo "$HYPHY_MODEL" | tr '[:upper:]' '[:lower:]')
VALID_MODELS=("busted" "absrel") 
if [[ ! " ${VALID_MODELS[*]} " =~ " ${HYPHY_MODEL} " ]]; then
    echo "ERROR : Invalid HyPhy model '$HYPHY_MODEL'. Available : ${VALID_MODELS[*]}" >&2
    exit 1
fi

export HYPHY_MODEL

# to generate the GMT file, you can use the "http://biit.cs.ut.ee/gmt-helper/#" website
if [ -n "$GMT_FILE" ] && [ -f "$GMT_FILE" ]; then
    export GMT_FILE=$(realpath "$GMT_FILE")
    echo "GMT file provided. Gene ontology enrichment will run."
else
    export GMT_FILE=""
fi

if [[ "$VCF" == *.vcf ]]; then
    export VCF_COMPRESSED="${OUTPUT_DIR}/$(basename "$VCF").gz"
else
    export VCF_COMPRESSED="$VCF"
fi

VCF_BASE=$(basename "$VCF" .gz)
VCF_BASE=${VCF_BASE%.vcf} 

# convert in case of relative paths 

# naming filtered VCF (to contain filters info) 
# to be changed in case the fixed filters are changed
export VCF_FILTERED="${OUTPUT_DIR}/${VCF_BASE}_SNP_MAFsup1percent_minQUAL30_minDP10_MISS${MISSING_FILTER}.vcf.gz"
export GFF_SORTED="${OUTPUT_DIR}/$(basename "$GFF" .gff)_sorted.gff"
export FASTA_FINAL_DIR="${OUTPUT_DIR}/FASTA_BY_GENE_CDS"

# gff output naming
if [ -n "$GFF" ]; then
    export GFF_SORTED="${OUTPUT_DIR}/$(basename "$GFF" .gff)_sorted.gff"
fi

export FASTA_FINAL_DIR="${OUTPUT_DIR}/FASTA_BY_GENE_CDS"



if [ -n "$SPECIES_FILE" ]; then
    export SPECIES_FILE=$(realpath "$SPECIES_FILE")
fi

if [ "$GEA_FILE" != "no" ] && [ -n "$GEA_FILE" ]; then
    export GEA_FILE=$(realpath "$GEA_FILE")
fi

########################################################################

# variables info 
echo -e "\n############################################\n"
echo "Outputs directory : $OUTPUT_DIR"
echo "HyPhy model       : $HYPHY_MODEL"
echo "Starting step     : $START_STEP"
echo "Ending step       : $END_STEP"
echo "Threads           : $THREADS"
echo -e "\n############################################\n"



# function allows to skip certain steps or stop before the pipeline finishes
# this might be usefil in the case of interrupted runs / specific needs
run_step() {
    local step_num="$1" desc="$2" commands="$3"
    
    local should_run=$(awk -v start="$START_STEP" -v end="$END_STEP" -v curr="$step_num" \
        'BEGIN { if (curr >= start && curr <= end) print 1; else print 0 }')

    if [ "$should_run" -eq 1 ]; then
        echo -e "\n$step_num. $desc"
        eval "$commands"  # bash -c $commands 
    else
        echo "Skipping step $step_num ($desc) ..."
    fi
}


########################################################################


LOG_FILE="${OUTPUT_DIR}/pipeline_run.log"
# redirect errors to terminal & log; stderr to stdout
exec > >(tee -a "$LOG_FILE") 2>&1
echo "Pipeline log in $LOG_FILE"

########################################################################

# pipeline steps execution according to the option input ( --start-at <1-7> and --end-at <1-7> )

# Pretreatments step includes indexing of FASTA and VCF, VCF compression, VCF filtering, GFF sorting
run_step 1 "=================== 1. Running pretreatments ===================" \
    "bash $SCRIPT_DIR/scripts/01_inputs_pretreatments.sh"


# VCF to FASTA conversion script runs the VCF2FASTA software & separates monomorphic genes to save analysis time
run_step 2 "=================== 2. VCF to FASTA conversion ===================" \
    "bash $SCRIPT_DIR/scripts/02_vcf2fasta.sh"

# Processing STOP codons is added to fight HyPhy's interruptions
run_step 3 "=================== 3. Processing STOP codons ===================" "
    # trimming the last codons from each sequence
    echo '=== 3.1 Processing terminal STOP codons ... ==='
    bash '$SCRIPT_DIR/scripts/03_1_final_stop_trim.sh'
    
    # strategy to increase the number of processed genes which treates genes where at least one premature STOP codon is introduced
    # by either deleting individuals (from tree and fasta) or trimming sequences to the last stop codon (not including)
    # statistics are generated for exploiting the stop codon occurance for each individual gene

    echo '=== 3.2 Processing premature STOP codons ... ==='
    export PRUNING_INPUT_DIR='${OUTPUT_DIR}/SEQS_FOR_HYPHY'
    bash '$SCRIPT_DIR/scripts/03_2_internal_stop_conditional_pruning_trimming_stats.sh'"

# Running HyPhy Analysis using the available CPUs (or < -t > option) & CPU restriction to limit computational times (data merging is causing computational times to expand) TD
run_step 4 "=================== 4. Running HyPhy Analysis with $HYPHY_MODEL model ===================" \
    "bash $SCRIPT_DIR/scripts/04_Hyphy.sh"



# extracting relevant information from the log files generated for each gene / branch
case "$HYPHY_MODEL" in
    "busted")
        run_step 5 "=================== 5. Assembling BUSTED gene results ===================" \
            "bash $SCRIPT_DIR/scripts/05_Hyphy_BUSTED_results.sh" ;;
    "absrel")
        run_step 5 "=================== 5. Assembling aBSREL branch results ===================" \
            "bash $SCRIPT_DIR/scripts/05_Hyphy_aBSREL_results.sh" ;;
    *) # td
esac


# only for busted:
if [ "$HYPHY_MODEL" == "busted" ]; then
    # Plotting significant results (p-val < 0.05) across contigs
    export HYPHY_RESULTS_FILE="${OUTPUT_DIR}/HYPHY_selection_results.tsv"
    run_step 6 "=================== 6. Plotting significant results ===================" \
        "python3 '$SCRIPT_DIR/scripts/06_Manhattan_positive_selection.py' --stats '$HYPHY_RESULTS_FILE' --gff '$GFF_SORTED' --gea '$GEA_FILE'"

    # GO term enrichment analysis
    if [ -n "$GMT_FILE" ]; then
        FDR_FILE="${OUTPUT_DIR}/genes_passing_FDR.tsv"
        if [ -f "$FDR_FILE" ] && [ "\$(wc -l < "$FDR_FILE")" -gt 1 ]; then
            run_step 7 "=================== 7. Gene Ontology Enrichment ===================" \
                "Rscript '$SCRIPT_DIR/scripts/07_GOenrichment.R' '$FDR_FILE' '$GMT_FILE'"
        else
            echo "WARNING : No significant genes found. Skipping GO enrichment."
        fi
    fi
else
    echo -e "\n=================== Skipping 6 & 7 ==================="
    echo "Manhattan plotting and GO enrichment are configured for the BUSTED model outputs."
fi

echo -e "\n#########################################################################################"
echo "======================    PIPELINE FINISHED SUCCESSFULLY.    ======================== "

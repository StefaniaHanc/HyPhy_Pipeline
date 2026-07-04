#!/bin/bash
set -euo pipefail
trap 'echo "ERROR : Step 5 failed at line $LINENO. Exiting." >&2' ERR

echo "============ Assembling and summarizing HyPhy results ... =============="

echo ""

########################################################################

LOG_DIR="${OUTPUT_DIR}/HYPHY_OUTPUTS/logs"
RESULTS_DIR="${OUTPUT_DIR}"
HYPHY_RESULTS_FILE="${RESULTS_DIR}/HYPHY_selection_results.tsv"
MONOMORPH_FILE="${OUTPUT_DIR}/FASTA_BY_GENE_CDS/variant_counts_summary.tsv"
OUTPUT_LOG="${OUTPUT_DIR}/Hyphy_summary_log.txt"

mkdir -p "$RESULTS_DIR"

if [[ ! -d "$LOG_DIR" ]]; then
    echo "ERROR : Log directory $LOG_DIR does not exist. Did HyPhy run and complete?" >&2
    exit 1
fi

########################################################################

echo "Parsing HyPhy logs into $HYPHY_RESULTS_FILE ..."

FORMAT="%-18s %-15s %-7s %-10s %-12s %-10s %-12s %-14s %-12s %-12s %-25s %-15s %-5s\n"
printf "$FORMAT" "GENE_NAME" "SELECTION" "SIGNIF" "P_VALUE" "OMEGA" "PROPORTION" "GLOBAL_OMEGA" "LOG_L" "AIC_C" "ERROR" "WARNING" "JSON_PVAL" "GEA" > "$HYPHY_RESULTS_FILE"

# nullglob skips if no .terminal.log files exist 
shopt -s nullglob

for filepath in "$LOG_DIR"/*.terminal.log; do
    gene_name=$(basename "$filepath" .terminal.log)
    dir_name=$(dirname "$filepath")
    json_file="$dir_name/${gene_name}.${HYPHY_MODEL}.json"
    
    # extract JSON p-value (full)
    if [[ -f "$json_file" ]]; then
        json_pval=$(grep -m 1 -oP '"p-value":\K[0-9.eE+-]+' "$json_file" || echo "NA")
    else
        json_pval="NA"
    fi
    [[ -z "$json_pval" ]] && json_pval="NA"

    # -w = whole word match, -F = string match
    gea_val=0
    if [ "$GEA_FILE" != "no" ] && [ -f "$GEA_FILE" ]; then
        if grep -q -w -F "$gene_name" "$GEA_FILE"; then gea_val=1; fi
    fi
########################################################################

    # process logs 
    awk -v gene="$gene_name" -v json_pval="$json_pval" -v gea_val="$gea_val" -v warning_phrase="derived from a single codon site" -v fmt="$FORMAT" '
    
    BEGIN {
        signif = "NA"; pval = "NA"; 
        warning = "NA"; error = "NA";
        section = ""; global_omega = "NA";
        omega_full = "NA"; prop_full = "NA";
        omega_constrained = "NA"; prop_constrained = "NA";
    }

    # SECTION

    /Obtaining the global omega/          { section = "global_omega" }
    /Performing the full/                 { section = "full_model" }
    /Performing the constrained/          { section = "constrained_model" }

########################################################################

    # SIGNIFICANCE AND P-VALUE

    /Likelihood ratio test.*p[ \t]*=[ \t]*[0-9.]+/ {
        match($0, /p[ \t]*=[ \t]*[0-9.]+/)
        pval_str = substr($0, RSTART, RLENGTH)
        split(pval_str, p_arr, "=")
        raw_pval = p_arr[2]; gsub(/^[ \t]+|[ \t]+$/, "", raw_pval)
        pval = sprintf("%.4f", raw_pval)
        signif = (raw_pval < 0.05) ? 1 : 0
    }

########################################################################

    # OMEGA & PROPORTIONS

    /Diversifying selection/ && section == "full_model" {
        split($0, arr, "|")
        if (arr[3] != "") {
            omega_full = arr[3]; gsub(/^[ \t]+|[ \t]+$/, "", omega_full)
            prop_full  = arr[4]; gsub(/^[ \t]+|[ \t]+$/, "", prop_full)
        }
    }

    /Neutral evolution/ && section == "constrained_model" {
        split($0, arr, "|")
        if (arr[3] != "") {
            omega_constrained = arr[3]; gsub(/^[ \t]+|[ \t]+$/, "", omega_constrained)
            prop_constrained  = arr[4]; gsub(/^[ \t]+|[ \t]+$/, "", prop_constrained)
        }
    }

########################################################################

    # GLOBAL OMEGA

    /\* non-synonymous\/synonymous rate ratio for \*test\*/ {
        global_omega = $NF
    }

    # LOG(L) & AIC-c

    /^\* Log\(L\)/ {
        match($0, /Log\(L\)[ \t]*=[ \t]*[-0-9.]+/)
        if (RSTART) {
            temp_logL = substr($0, RSTART, RLENGTH)
            split(temp_logL, l_arr, "=")
            val_logL = l_arr[2]; gsub(/[ \t]+/, "", val_logL)
            if (section == "full_model") logL_full = val_logL
            if (section == "constrained_model") logL_constrained = val_logL
        }
        
        match($0, /AIC-c[ \t]*=[ \t]*[-0-9.]+/)
        if (RSTART) {
            temp_aic = substr($0, RSTART, RLENGTH)
            split(temp_aic, a_arr, "=")
            val_aic = a_arr[2]; gsub(/[ \t]+/, "", val_aic)
            if (section == "full_model") aic_full = val_aic
            if (section == "constrained_model") aic_constrained = val_aic
        }
    }

########################################################################

    # WARNINGS & ERRORS

    $0 ~ warning_phrase { warning = "Single_codon_derived" }
    tolower($0) ~ /stop codons shown/ { error = "StopCodon" }
    tolower($0) ~ /sequences are not properly aligned/ { error = "Not_Properly_Aligned" }
    tolower($0) ~ /internal error in computebranchcache/ { error = "ComputeBranchCache" }

########################################################################

    # FINAL

    END {
        gsub(/[\r\n]/, "", error)
        
        if (signif == 1) {
            final_logL  = logL_full
            final_aic   = aic_full
            final_omega = omega_full
            final_prop  = (prop_full != "NA" && prop_full != "") ? prop_full "%" : "NA"
            selection   = "Diversifying"
        } else {
            final_logL  = logL_constrained
            final_aic   = aic_constrained

            if (omega_constrained != "NA" && omega_constrained != "") {
                final_omega = omega_constrained
                final_prop  = prop_constrained "%"
            } else {
                final_omega = global_omega
                final_prop  = "NA"
            }
            
            if (global_omega == "" || global_omega == "NA") {
                selection = "NA"
            } else if (global_omega < 0.95) {
                selection = "Purifying"
            } else if (global_omega > 1.05) {
                selection = "Diversifying_NS" 
            } else {
                selection = "Neutral"
            }
        }
        
        if (final_omega == "") final_omega = "NA"
        if (final_logL == "") final_logL = "NA"
        if (final_aic == "") final_aic = "NA"
        if (global_omega == "") global_omega = "NA"

        printf fmt, gene, selection, signif, pval, final_omega, final_prop, global_omega, final_logL, final_aic, error, warning, json_pval, gea_val
    }
    ' "$filepath" >> "$HYPHY_RESULTS_FILE"


done
shopt -u nullglob


########################################################################

# multiple testing
echo "Bonferroni and FDR (Benjamini-Hochberg) corrections ..."

CORRECTED_TSV_FDR="${RESULTS_DIR}/HYPHY_selection_corrected_FDR.tsv"
CORRECTED_TSV_BONF="${RESULTS_DIR}/HYPHY_selection_corrected_Bonferroni.tsv"
ALIGNED_FILE_FDR="${RESULTS_DIR}/genes_passing_FDR.tsv"
ALIGNED_FILE_BONF="${RESULTS_DIR}/genes_passing_Bonferroni.tsv"

# total tests (- header & 'NA')
TOTAL_TESTS=$(awk 'NR>1 && $4 != "NA"' "$HYPHY_RESULTS_FILE" | wc -l)

if [[ "$TOTAL_TESTS" -gt 0 ]]; then
    head -n 1 "$HYPHY_RESULTS_FILE" > "$CORRECTED_TSV_FDR"
    head -n 1 "$HYPHY_RESULTS_FILE" > "$CORRECTED_TSV_BONF"
    
    # sort pval by ascending (col 4)
    awk 'NR>1 && $4 != "NA"' "$HYPHY_RESULTS_FILE" | sort -k4,4n > "${RESULTS_DIR}/tmp_sorted.tsv"
    # calc of thresholds
    awk -v N="$TOTAL_TESTS" -v alpha=0.05 -v fdr_out="$CORRECTED_TSV_FDR" -v bonf_out="$CORRECTED_TSV_BONF" '
    BEGIN { 
        max_fdr_rank = 0 
        bonf_thresh = alpha / N
    }
    {
        pval = $4
        fdr_thresh = (NR / N) * alpha
        
        # FDR
        if ((pval+0) <= (fdr_thresh+0)) {
            max_fdr_rank = NR
        }
        
        # Bonferroni
        if (pval+0 <= bonf_thresh+0) {
            print $0 >> bonf_out
        }
        
        lines[NR] = $0
    }
    END {
        # print lines passing FDR
        for (i=1; i<=NR; i++) {
            if (i <= max_fdr_rank) {
                print lines[i] >> fdr_out
            }
        }
    }' "${RESULTS_DIR}/tmp_sorted.tsv"
    
    rm -f "${RESULTS_DIR}/tmp_sorted.tsv"

    # format
    awk '{$1=$1; print}' "$CORRECTED_TSV_FDR" | column -t > "$ALIGNED_FILE_FDR"
    awk '{$1=$1; print}' "$CORRECTED_TSV_BONF" | column -t > "$ALIGNED_FILE_BONF"
    
    echo "Genes passing FDR correction : $ALIGNED_FILE_FDR"
    echo "Genes passing Bonferroni correction : $ALIGNED_FILE_BONF"
else
    echo "No valid p-values found."
fi

########################################################################

# RESULTS SUMMARY CALCULS
echo "Calculating summary stats ..."

calc_pct() {
    local part=$1 total=$2
    if [[ "$total" -eq 0 || -z "$total" ]]; then
        echo "0.00%"
    else
        awk -v p="$part" -v t="$total" 'BEGIN { printf "%.2f%%", (p/t)*100 }'
    fi
}

analyze_hyphy_results() {
    local file_path="$1" monomorph_count="$2" 

    local total_lines=$(wc -l < "$file_path")
    local total_genes=$((total_lines - 1))
    
    local pos_sel=$(awk '$2=="Diversifying"' "$file_path" | wc -l)
    local pos_sel_ns=$(awk '$2=="Diversifying_NS"' "$file_path" | wc -l)
    local neutral_sel=$(awk '$2=="Neutral"' "$file_path" | wc -l)
    local purifying_sel=$(awk '$2=="Purifying"' "$file_path" | wc -l)

    local pos_sel_single_codon=$(awk '$2=="Diversifying" && $0 ~ /Single_codon_derived/' "$file_path" | wc -l)
    local not_pos_sel_single_codon=$(awk '$3=="0" && $0 ~ /Single_codon_derived/' "$file_path" | wc -l)
    local single_codon_total=$(awk '$0 ~ /Single_codon_derived/' "$file_path" | wc -l)
    
    local internal_stops=$(awk '$0 ~ /StopCodon/' "$file_path" | wc -l)
    local no_sel=$(awk '$2 ~ /NA/' "$file_path" | wc -l)
    local pos_sel_multi_codon=$((pos_sel - pos_sel_single_codon))

    echo -e "\n################################################################\n"
    echo " Results Summary"
    echo -e "\n################################################################\n"
    echo "File source : $file_path"
    echo "Total genes : $total_genes"
    echo "Genes under significant positive selection (total) : $pos_sel ($(calc_pct "$pos_sel" "$total_genes"))"
    echo "Genes under positive selection, but not significant  : $pos_sel_ns ($(calc_pct "$pos_sel_ns" "$total_genes"))"
    echo "Genes under neutral selection : $neutral_sel ($(calc_pct "$neutral_sel" "$total_genes"))"
    echo "Genes under purifying selection : $purifying_sel ($(calc_pct "$purifying_sel" "$total_genes"))"

    echo -e "\n################################################################\n"
    echo "WARNING : Single_codon_derived"
    echo -e "\n################################################################\n"
    echo "Genes under positive selection derived from single codon site : $pos_sel_single_codon "
    echo "Genes under purifying/neutral derived from single codon site : $not_pos_sel_single_codon "
    echo "Genes derived from single codon site (total of warnings) : $single_codon_total "
    
    echo -e "\n################################################################\n"
    echo "Internal stop codons : $internal_stops"
    echo "No results (should correspond to internal stops nb) : $no_sel"
    
    if [[ -n "$monomorph_count" && "$monomorph_count" -gt 0 ]]; then
        echo -e "\n################################################################\n"
        echo "Monomorph genes excluded prior to HyPhy : $monomorph_count"
        echo -e "\n################################################################\n"
    fi
}

# monomorph counts
monomorph_total=0
if [[ -f "$MONOMORPH_FILE" ]]; then
    monomorph_total=$(awk '$2=="0"' "$MONOMORPH_FILE" | wc -l)
fi

########################################################################

# output log
analyze_hyphy_results "$HYPHY_RESULTS_FILE" "$monomorph_total" > "$OUTPUT_LOG"

echo "Results generation complete. Results :"
echo "Details : $HYPHY_RESULTS_FILE"
echo "Corrected p-val : $ALIGNED_FILE_FDR"
echo "Corrected p-val : $ALIGNED_FILE_BONF"
echo "Summary file : $OUTPUT_LOG"
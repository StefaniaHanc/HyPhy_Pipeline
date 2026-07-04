#!/bin/bash

set -euo pipefail
trap 'echo "ERROR : aBSREL summary failed at line $LINENO. Exiting." >&2' ERR

echo "============ Assembling aBSREL results ... =============="

ABSREL_RESULTS="${OUTPUT_DIR}/HYPHY_aBSREL_selection_results.tsv"
LOGS="${OUTPUT_DIR}/HYPHY_OUTPUTS/logs"

if [ -z "${SPECIES_FILE:-}" ] || [ ! -f "$SPECIES_FILE" ]; then
    echo "ERROR: SPECIES_FILE is not defined or does not exist." >&2
    exit 1
fi

FORMAT="%-18s %-15s %-20s %-12s %-10s %-12s %-40s %-60s\n"
printf "$FORMAT" "GENE_NAME" "SELECTION" "NB_BRANCHES_POS_SEL" "AVG_P_VAL" "TESTED_NB" "ERRORS" "SPECIES_UNDER_SEL" "BRANCHES_POS_SEL" > "$ABSREL_RESULTS"

shopt -s nullglob


for gene in "$LOGS"/*.terminal.log; do
    gene_name=$(basename "$gene" .terminal.log)

    # aBSREL info
    awk -v gene="$gene_name" -v fmt="$FORMAT" '
    BEGIN {
        # vars init
        selection = "NA"; avg_p_val = "NA"; num_tested = "NA"; 
        nb_branches_pos_sel = 0; error = "NA"; 
        capture_branches = 0; sum_p_val = 0; count_p_val = 0;
    }
    
    # read species appartenance first
    NR==FNR {
        if ($1 != "" && $2 != "") {
            species_map[$1] = $2
        }
        next
    }

    # terminal log file
    
    # stop codon genes detection
    /[Ss]top codons shown in capital letters/ { error = "StopCodon" }
    
    # execution errors (outside StopCodon)
    /Check errors\.log for execution error details\./ { 
        if (error != "StopCodon") { error = "ExecError" }
    }
    
    # other error
    (/ERROR/ || /[Ee]rror/ || /Exception/ || /Segmentation fault/ || /Assertion/) { 
        if (error != "StopCodon" && error != "ExecError") { error = "Yes" }
    }
    
    # summary line 
    /found \*\*[0-9]+\*\* branches under selection/ {
        split($0, arr, "\\*\\*")
        num_sel = arr[2]
        num_tested = arr[4]

        if (num_sel > 0) {
            selection = "Diversifying"
            capture_branches = 1
        } else {
            selection = "Other"
            capture_branches = 0
        }
    }

    # branches under selection and calculate average pval
    /^\* .+, p-value = / {
        if (capture_branches == 1) {
            branch = $2; gsub(/,/, "", branch) # branch name
            pval = $5 + 0 # make numeric
            
            # lookup species or extract node number
            if (branch in species_map) {
                sp = substr(species_map[branch], 1, 3)
            } else if (branch ~ /Node/ ~ /^[0-9]+$/) {
                sp = branch
                gsub(/[^0-9]/, "", sp)
                if (sp == "") sp = branch
            } else {
                sp = "NA"
            }

            # tracking
            sum_p_val += pval
            count_p_val++
            nb_branches_pos_sel++
            
            # store data for branch sorting later
            b_pvals[branch] = pval
            
            # store data for species average calculation
            sp_sum[sp] += pval
            sp_count[sp]++
        }
    }

    END {
        # calculate overall gene average pval
        if (count_p_val > 0) {
            avg_p_val = sprintf("%.5f", sum_p_val / count_p_val)
        }

        # SORTING by pval

        # process & sort branches by pval - ascending
        b_idx = 1
        for (b in b_pvals) { b_array[b_idx++] = b }
        total_b = b_idx - 1
        
        # sort branches
        for (i = 1; i <= total_b; i++) {
            for (j = i + 1; j <= total_b; j++) {
                if (b_pvals[b_array[i]] > b_pvals[b_array[j]]) {
                    tmp = b_array[i]; b_array[i] = b_array[j]; b_array[j] = tmp
                }
            }
        }
        
        # branch:p_val
        branches_pos_sel = ""
        for (i = 1; i <= total_b; i++) {
            b = b_array[i]
            val = sprintf("%s:%g", b, b_pvals[b])
            branches_pos_sel = (branches_pos_sel == "") ? val : branches_pos_sel "," val
        }
        if (branches_pos_sel == "") branches_pos_sel = "NA"


        # process and sort species by average pval - ascending
        s_idx = 1
        for (sp in sp_sum) {
            sp_avg[sp] = sp_sum[sp] / sp_count[sp]
            s_array[s_idx++] = sp
        }
        total_s = s_idx - 1

        # sort species
        for (i = 1; i <= total_s; i++) {
            for (j = i + 1; j <= total_s; j++) {
                if (sp_avg[s_array[i]] > sp_avg[s_array[j]]) {
                    tmp = s_array[i]; s_array[i] = s_array[j]; s_array[j] = tmp
                }
            }
        }

        # species:avg_pval
        species_sel_str = ""
        for (i = 1; i <= total_s; i++) {
            sp = s_array[i]
            val = sprintf("%s:%.5f", sp, sp_avg[sp])
            species_sel_str = (species_sel_str == "") ? val : species_sel_str "," val
        }
        if (species_sel_str == "") species_sel_str = "NA"


        # line/gene
        printf fmt, gene, selection, nb_branches_pos_sel, avg_p_val, num_tested, error, species_sel_str, branches_pos_sel
    }
    ' "$SPECIES_FILE" "$gene" >> "$ABSREL_RESULTS"
done
shopt -u nullglob

echo "DONE . aBSREL branch/species results in : $ABSREL_RESULTS"
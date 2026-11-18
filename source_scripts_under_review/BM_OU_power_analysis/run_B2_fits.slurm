#!/bin/bash
#SBATCH --job-name=B2_fits
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --output=slurm-B2-fits-%j.out
#SBATCH --error=slurm-B2-fits-%j.err

set -uo pipefail

PROJECT_DIR="$HOME/projects/uosa/Carolin_Kosiol/manel_aitelhadj/BM_OU_analysis_all_species"
RB="/mnt/apps/users/maehadj/conda/envs/revbayes_env/bin/rb"

cd "$PROJECT_DIR" || exit 1
mkdir -p scripts_scratch/B2_fits scripts_scratch/logs

run_cell() {
    local n_tips="$1" label="$2" rep="$3"
    local sim_dir="scripts_scratch/B2_sim_data/N${n_tips}_${label}/rep${rep}"
    local tree_path="$sim_dir/tree.nwk"
    local data_path="$sim_dir/data.nex"

    for model in bm ou; do
        local result_root="scripts_scratch/B2_fits/N${n_tips}_${label}_rep${rep}_${model}"
        if [[ -f "${result_root}/marginal_likelihood.txt" ]]; then
            echo "SKIP N${n_tips} ${label} rep${rep} ${model}"
            continue
        fi
        "$RB" > "scripts_scratch/logs/fit_N${n_tips}_${label}_rep${rep}_${model}.log" 2>&1 <<EOF_RB
tree_path = "$tree_path"
data_path = "$data_path"
result_root = "$result_root"
source("scripts_scratch/${model}_simulation_fit.rev")
EOF_RB
        [[ $? -eq 0 ]] && echo "[$(date '+%T')] END N${n_tips} ${label} rep${rep} ${model}" \
                        || echo "[$(date '+%T')] FAILED N${n_tips} ${label} rep${rep} ${model}" >&2
    done
}
export -f run_cell
export RB

parallel --colsep ' ' --jobs "$SLURM_NTASKS" \
    --joblog "B2_fits_${SLURM_JOB_ID}.log" \
    run_cell {1} {2} {3} :::: scripts_scratch/B2_joblist.txt

echo "Finished B2 fits at $(date)"

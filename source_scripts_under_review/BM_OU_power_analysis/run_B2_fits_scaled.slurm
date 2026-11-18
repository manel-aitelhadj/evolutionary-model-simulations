#!/bin/bash
#SBATCH --job-name=B2_fits_scaled
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --output=slurm-B2-fits-scaled-%j.out
#SBATCH --error=slurm-B2-fits-scaled-%j.err
set -uo pipefail
PROJECT_DIR="$HOME/projects/uosa/Carolin_Kosiol/manel_aitelhadj/BM_OU_analysis_all_species"
RB="/mnt/apps/users/maehadj/conda/envs/revbayes_env/bin/rb"
cd "$PROJECT_DIR" || exit 1
mkdir -p scripts_scratch/B2_fits_scaled scripts_scratch/logs

run_cell() {
    local n_tips="$1" label="$2" rep="$3"
    local sim_dir="scripts_scratch/B2_sim_data_scaled/N${n_tips}_${label}/rep${rep}"
    for model in bm ou; do
        local result_root="scripts_scratch/B2_fits_scaled/N${n_tips}_${label}_rep${rep}_${model}"
        if [[ -f "${result_root}/marginal_likelihood.txt" ]]; then continue; fi
        "$RB" > "scripts_scratch/logs/fit_scaled_N${n_tips}_${label}_rep${rep}_${model}.log" 2>&1 <<EOF_RB
tree_path = "$sim_dir/tree.nwk"
data_path = "$sim_dir/data.nex"
result_root = "$result_root"
source("scripts_scratch/${model}_simulation_fit.rev")
EOF_RB
    done
}
export -f run_cell
export RB
parallel --colsep ' ' --jobs "$SLURM_NTASKS" --joblog "B2_fits_scaled_${SLURM_JOB_ID}.log" \
    run_cell {1} {2} {3} :::: scripts_scratch/B2_joblist_scaled.txt
echo "Finished B2 scaled fits at $(date)"

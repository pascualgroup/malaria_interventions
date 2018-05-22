# Organization
Experiments are organized in 4 levels:

1. **Parameter space**: This is the set of parameters such as biting rates, interventions, etc.
2. **Scenario**: Each parameter space can be run for the 3 scenarios.
3. **Experiment**: Experiments are variations on the parameter space. By definition '00' denotes the basic parameter space that is used to reach a steady-state and create a checkpoint. '01' is control. It has the *exact same parameters* as '00'. It loads the checkpoint created by 00 and continues the simulation. Any other number (e.g., 02, 03...) denotes a variation from 00 (e.g., to test interventions, different biting rates or different seasonal patterns).
4. **Runs**: Each combination of paramter space, scenario and experiment can be run multiple times. Importantly, the seed for the random number (`RANDOM_SEED` parameter) is consistent across all experiments for a given run in a given parameter space. For example, within parameter space PS03, any run *r* of experiments E00, E01, E02,... will have the same value for `RANDOM_SEED`. I did this because E01 (control), should be an exact continuation of E00 (checkpoint).

The study design is layed in [this google sheet](https://docs.google.com/spreadsheets/d/1AetmLv-3sxpv9blupDA04pF_Y0RYOod38sRxGu1SOuM/edit?usp=sharing).

# File names
**Parameter file names** are denoted as: `PSxx_y_Ezz_Rw.py`, where *xx* is the id of the parameter space, *y* is the scenario ('S', 'N', or 'G'), *zz* is the id of the experiment and *w* is the id of the run. Note that *w* does not have a leading 0 as it is an integer, while *xx* and *zz* can have as they are strings.

**job file names** are similarly denoted: `PSxxyEzz.sbatch` (no underscores so the job name will be shorter). The job file name does not have the run because the runs are run as arrays in slurm. The runs to run are specified in the sbatch: `#SBATCH --array='1-10'`.

**sqlite file names** are `PSxx_y_Ezz.sqlite`. If an sqlite is a checkpoint it will be: `PSxx_y_E00_CP.py` (00 is always the experiment name for creating a checkpoint).

# General workflow pipeline
The neutral scenarios are always counterpart to the immune selection scenario. Therefore, their paramter files can only be created AFTER the counterpart immune selection scenario is run. Within each scenario, experiments depend on E00 because the system needs to reach a steady state. Therefore, E00 should be run first and when it finished the rest can be run. Hence, a typical pipeline should look lije this:

1. Generate parameter files for all experiments in a given paramter space (e.g., 01) in S.
2. Run `PS01SE00.sbatch` for `PS01_S_E00.py`
3. Run the rest of the experiments that depend on E00. This can be done manually, wating for `PS01SE00.sbatch` to finish or with: `sbatch -d afterok:jobid PS01SE01.sbatch`, where jobid is the job id of `PS01SE00.sbatch`
4. Copy the resulting sqlite files to my local computer: `cp PS01_S*.sqlite shai@192.170.193.140:/home/shai/Documents/malaria_interventions_sqlite`.
5. Generate parameter files for the N and G scenarios based on `PS01_S_E01_Rw`. Note that only the control experiemt (E01) acn be used for that because the N and G are created by matching duration of infection of the counterpart S experiment. In the R code, the neutral scenario is matched by averaging the duration of infection. This averaging is done across all available runs for E01 (i.e., PS01_S_E01_R1, PS01_S_E01_R2, PS01_S_E01_R3, ...). In contrast, the G scenario is matched by fitting a curve and this is done per run. Therefore, PS01_G_Exx_R**1** will match PS01_S_E01_R**1**, PS01_G_Exx_R**2** will match PS01_S_E01_R**2**, etc (always match to E01, separately per run).
6. Run `PS01NE00.sbatch` for `PS01_N_E00.py` (and analogously for the G scenario).
7.  Run the rest of the experiments that depend on E00 (as in step 3) for N and G.

# Seasonality
Seasonality is implemented using an ODE model for the mosquito population.

# Interventions

## IRS

## MDA

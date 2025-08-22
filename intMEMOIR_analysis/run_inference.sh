#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --job-name=inference_adb_fixedTrees_condRoot_withRho
#SBATCH --array=1-5
#SBATCH --output=outs/%x_%a.out

SEED=$SLURM_ARRAY_TASK_ID
analysis=$SLURM_JOB_NAME

java -Dglass.platform=Monocle -Dmonocle.platform=Headless \
--module-path $HOME/javafx-sdk-17.0.6-linux-monocle/lib --add-modules=javafx.base,javafx.fxml \
-jar $HOME/beast_runs/ADB_06-08-25.jar \
-version_file $HOME/beast_runs/version_files/adb_version.xml \
-version_file $HOME/beast_runs/version_files/feast_version.xml \
-version_file $HOME/beast_runs/version_files/beast_version.xml \
-version_file $HOME/beast_runs/version_files/tidetree_version.xml \
-overwrite \
-seed $SEED \
-loglevel debug \
-statefile ${analysis}/chain${SEED}.state \
${analysis}.xml > ${analysis}/chain${SEED}.out

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --job-name=treesize_comparison
#SBATCH --array=1-40
#SBATCH --output=outs/%A_%a.out

i=$SLURM_ARRAY_TASK_ID
treeNr=$(sed "$((i+1))q;d" "tree_data.csv" | cut -d ',' -f 1)
nTips=$(sed "$((i+1))q;d" "tree_data.csv" | cut -d ',' -f 2)
origin=$(sed "$((i+1))q;d" "tree_data.csv" | cut -d ',' -f 3)
string="treeNr=$treeNr,nTips=$nTips,origin=$origin"
echo "$string"

java -Dglass.platform=Monocle -Dmonocle.platform=Headless \
--module-path $HOME/javafx-sdk-17.0.6-linux-monocle/lib --add-modules=javafx.base,javafx.fxml \
-jar $HOME/beast_runs/ADB_06-08-25.jar \
-version_file $HOME/beast_runs/version_files/adbp_version.xml \
-version_file $HOME/beast_runs/version_files/feast_version.xml \
-version_file $HOME/beast_runs/version_files/beast_version.xml \
-overwrite \
-seed 1 \
-loglevel debug \
-D "$string" \
-statefile inference/inference_n${nTips}_${treeNr}.state \
inference.xml >> inference/inference_n${nTips}_${treeNr}.out



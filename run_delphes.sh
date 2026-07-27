#!/bin/bash
#SBATCH --job-name=delphes_mVII
#SBATCH --partition=tartarus
#SBATCH --output=delphes_mVII.out
#SBATCH --error=delphes_mVII.err
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

# --- Configuration ---
DELPHES_EXE=~/Titan0/madgraph/MG5_aMC_v3_7_0/Delphes/DelphesHepMC2
# DELPHES_CARD=~/Titan0/delphes/cards/ATLAS_Run2/delphes_card_atlas_bare.tcl
DELPHES_CARD=~/Titan0/delphes/cards/ATLAS_Run2/delphes_card_atlas_dm_mIII.tcl
# DELPHES_CARD=~/Titan0/delphes/cards/ATLAS_Run2/delphes_card_atlas_conf_2020_002.tcl
# DELPHES_CARD=~/Titan0/delphes/cards/ATLAS_Run2/delphes_card_ATLAS_Run2.tcl
HEPMC_DIR=~/Titan0/pythia2/pythia8245/examples/Research/hepmc_lxplus
ROOT_DIR=~/Titan0/madgraph/MG5_aMC_v3_7_0/bin/ttbar_validation/Events/BSM_events

# List of tags to process
tags=("1000_6" "1000_8" "2000_4" "2000_6" "3000_2" "3000_4")

# --- Execution Loop ---
for tag in "${tags[@]}"
do
    hepmc_file="${HEPMC_DIR}/${tag}_events.hepmc"
    root_file="${ROOT_DIR}/${tag}_events.root"

    echo "-------------------------------------------"
    echo " - Processing Tag: $tag"
    echo "   Input:  $hepmc_file"
    echo "   Output: $root_file"
    
    # Check if input file exists before running
    if [ -f "$hepmc_file" ]; then
        $DELPHES_EXE "$DELPHES_CARD" "$root_file" "$hepmc_file"
    else
        echo "   [ERROR] Input file $hepmc_file not found. Skipping..."
    fi
done

echo "-------------------------------------------"
echo "Done!"
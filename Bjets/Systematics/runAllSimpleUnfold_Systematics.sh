#!/bin/bash

cd ../


commands=(
#JESJER
  "root -l -b -q 'MakeCorrections.C(-1, 91599, 0, 1, 0, 0, 0, 0)'"
#JetID
  "root -l -b -q 'MakeCorrections.C(-1, 91599, 0, 0, 1, 0, 0, 0)'"
#RecSelEff  
  "root -l -b -q 'MakeCorrections.C(-1, 91599, 0, 0, 0, 1, 0, 0)'"
#FitModel  
  "root -l -b -q 'MakeCorrections.C(-1, 91599, 0, 0, 0, 0, 1, 0)'"
#UnfoldPrior  
  "root -l -b -q 'MakeCorrections.C(-1, 91599, 0, 0, 0, 0, 0, 1)'"
)
	
# First create the MC weights
root -l -b -q MakeRMReweightFactors.C
	
# Create a new tmux session
tmux new-session -d -s mysessionMakeCorrections_Sys

# Loop through the command strings and create split windows in tmux
for command in "${commands[@]}"; do
    tmux split-window -h "$command"
    tmux select-layout tiled
done


cd Systematics/

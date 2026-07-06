#!/bin/bash

#srcroot8
#srcRoox8

commands=(
  "root -l -b -q 'MCMakeVarTree.C(-1, 61590)'"
  "root -l -b -q 'MCMakeVarTree.C(-1, 61591)'"
  "root -l -b -q 'MCMakeVarTree.C(-1, 71590)'"
  "root -l -b -q 'MCMakeVarTree.C(-1, 71591)'"
  "root -l -b -q 'MCMakeVarTree.C(-1, 81590)'"
  "root -l -b -q 'MCMakeVarTree.C(-1, 81591)'"
	)
# Create a new tmux session
tmux new-session -d -s mysessionMCMakeVarTree

# Loop through the command strings and create split windows in tmux
for command in "${commands[@]}"; do
    tmux split-window -h "$command"
    tmux select-layout tiled
done

# Attach to the tmux session
# tmux attach-session -t mysessionMCSimple



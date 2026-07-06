#!/bin/bash

commands=(
  "root -l -b -q 'MakeVarTree.C(-1, 61590, 0)'"
  "root -l -b -q 'MakeVarTree.C(-1, 61591, 0)'"
  "root -l -b -q 'MakeVarTree.C(-1, 71590, 0)'"
  "root -l -b -q 'MakeVarTree.C(-1, 71591, 0)'"
  "root -l -b -q 'MakeVarTree.C(-1, 81590, 0)'"
  "root -l -b -q 'MakeVarTree.C(-1, 81591, 0)'"
	)
# Create a new tmux session
tmux new-session -d -s mysessionMakeVarTree

# Loop through the command strings and create split windows in tmux
for command in "${commands[@]}"; do
    tmux split-window -h "$command"
    tmux select-layout tiled
done

# Attach to the tmux session
# tmux attach-session -t mysessionMCSimple



#!/bin/bash
# This runs on NUISANCE flat-trees and gets 3D cross-sections

for i in $(ls *.root_FlatTree.root); do
  root -l -b -q 'fsianal.cpp("'$i'")' &
done

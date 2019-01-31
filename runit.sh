#!/bin/bash
# This runs on NUISANCE flat-trees and gets 3D cross-sections

for i in $(ls *.root_FlatTree.root); do
  for j in CCQE 2p2h 1pi MpiDIS; do
    root -l -b -q 'fsianal.cpp("'$i'", "'$j'")' &
    #root -l -b -q 'fsianal_skew.cpp("'$i'")' &
  done
done

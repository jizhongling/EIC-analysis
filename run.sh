#!/bin/bash

proc="$1"
particle="$2"
N=10

#root -l -b -q "AnaClusters.C(${proc}, \"${particle}\")"
#root -l -b -q "AnaClustersForSiPM.C(${proc}, \"${particle}\")"
#root -l -b -q "AnaSiPM.C(${proc})"
for (( i = proc*N ; i < (proc+1)*N ; i++ )) ; do
  #root -l -b -q "AnaClusters.C(${i}, \"${particle}\")"
  #root -l -b -q "AnaHFJet.C(${i})"
  for (( j = 1 ; j <= 5 ; j++ )) ; do
    root -l -b -q "AnaJetID.C(${i}, ${j})"
  done
  #root -l -b -q "AnaRate.C(${i})"
done

#!/bin/bash

proc="$1"
particle="$2"
N=20

#root -l -b -q "AnaClusters.C(${proc}, \"${particle}\")"
#root -l -b -q "AnaClustersForSiPM.C(${proc}, \"${particle}\")"
#root -l -b -q "AnaSiPM.C(${proc})"
for (( i = proc*N ; i < (proc+1)*N ; i++ )) ; do
  root -l -b -q "AnaHFJet.C(${i})"
done

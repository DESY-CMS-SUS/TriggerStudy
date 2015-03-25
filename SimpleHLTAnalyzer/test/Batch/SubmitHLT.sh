#!/bin/zsh

Njobs=500
Nevents=500
OutDir=Output/TT_nojets
#OutDir=test

#qsub -t 1-$Njobs -o logs  HLTtupleJob.sh $Nevents $OutDir
qsub -t 2001-6000 -o logs  HLTtupleJob.sh $Nevents $OutDir

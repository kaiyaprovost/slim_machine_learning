#!/bin/bash

programPath="slim_model_selection_master_with_export.py"

## calculate for simulated data, train neural network fully, output
cd "~/simulated_statistics"
ls *.STATS | xargs -t python3 "${programPath}"

## calculate for empirical data after training neural network
cd "~/empirical_statistics"
nnfile="DEM-IBD-AGEtrained_neural_network.joblib"
empfile="cardcard16.fasta.FULLSTATS.COMBO.STATS"
trained="X_train.temp"
python3 "${programPath}" $nnfile $empfile $trained

#!/usr/bin/python -u
# coding: utf-8

# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/NONCORRELATED/"
# nnfile="trained_neural_network1602265327.joblib"
# empfile="DEM-IBD-AGE_COMBINED_PAN-YES-120K.COMBO.STATS" ## no spaces
# python3 "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/2.summary statistics to neural network/slim_model_selection_master_EMPIRICAL.py" $nnfile $empfile


import numpy as np
from joblib import dump, load
import sys
import os

import slim_model_selection_data_processing as sms_pre
import slim_model_selection_build_NN as sms_bnn
import slim_model_selection_summarize_NN as sms_snn

def main():

	try:
		trained_nn_file=sys.argv[1]
	except:
		print("Please give trained NN .joblib file")
		exit()
		trained_nn_file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/NONCORRELATED/trained_neural_network1602265327.joblib"
	
	try:
		empirical_data_file = sys.argv[2]
	except:
		print("Please give empircal data to fit")
		exit()
		empirical_data_file = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/NONCORRELATED/DEM-IBD-AGE_COMBINED_PAN-YES-120K.COMBO.STATS"
	
	clf = load(trained_nn_file)
	
	empirical_data = np.loadtxt(empirical_data_file)
	
	empirical_predicted = clf.predict(empirical_data)
	
	outfileString = trained_nn_file+"_"+os.path.basename(empirical_data_file)+"_PREDICTED.txt"
	np.savetxt(outfileString, empirical_predicted, delimiter="\t",fmt='%s')
	
if __name__ == "__main__":
	main()
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

from sklearn.preprocessing import StandardScaler ## 


def featureScalingEmpirical(X_train,X_empirical):
	# Don't cheat - fit only on training data
	scaler = StandardScaler(copy=False, with_mean=True, with_std=True).fit(X_train) 
	X_empirical_scaled = scaler.transform(X_empirical) 
	return(X_empirical_scaled)

def main():

	try:
		trained_nn_file=sys.argv[1]
	except:
		print("Please give trained NN .joblib file")
		exit()
		trained_nn_file="/home/kprovost/nas5/slim_osg/train_nn/noncorrelated/DEM-IBD-AGEtrained_neural_network_**NN_1602700103.joblib"
	try:
		empirical_data_file = sys.argv[2]
	except:
		print("Please give empircal data to fit")
		exit()
		empirical_data_file = "/home/kprovost/nas5/slim_osg/train_nn/MERGED_snakes_14oct2020_noncollinear_11_complete.txt"
	try:
		X_train_file = sys.argv[3]
	except:
		print("Please give trained data to scale on")
		exit()
		X_train_file = "/home/kprovost/nas5/slim_osg/train_nn/noncorrelated/X_train.**NN_1602700103.temp"
		
	clf = load(trained_nn_file)
	
	#clf2 = clf.estimators_[0].best_estimator_.coefs_

	empirical_data = np.loadtxt(empirical_data_file)
	X_train = np.loadtxt(X_train_file,delimiter=",")
	
	X_empirical_scaled = featureScalingEmpirical(X_train,empirical_data)
	
	empirical_predicted = clf.predict(X_empirical_scaled)
	
	outfileString = trained_nn_file+"_"+os.path.basename(X_train_file)+"_"+os.path.basename(empirical_data_file)+"_PREDICTED.txt"
	np.savetxt(outfileString, empirical_predicted, delimiter="\t",fmt='%s')
	
if __name__ == "__main__":
	main()
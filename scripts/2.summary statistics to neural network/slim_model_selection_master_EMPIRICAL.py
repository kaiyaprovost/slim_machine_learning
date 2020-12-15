#!/usr/bin/python -u
# coding: utf-8

# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/NONCORRELATED/"
# nnfile="DEM-IBD-AGEtrained_neural_network_**NN_1602701237.joblib"
# empfile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/cardcard16.fasta.FULLSTATS.COMBO.STATS" ## no spaces
# trained="X_train.**NN_1602701237.temp"
# python3 "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/2.summary statistics to neural network/slim_model_selection_master_EMPIRICAL.py" $nnfile $empfile $trained

import numpy as np
from joblib import dump, load
import sys
import os

import slim_model_selection_data_processing as sms_pre
import slim_model_selection_build_NN as sms_bnn
import slim_model_selection_summarize_NN as sms_snn

from sklearn.preprocessing import StandardScaler ## 
from sklearn.model_selection import learning_curve
from sklearn.model_selection import ShuffleSplit

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
		#exit()
		trained_nn_file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib"
	try:
		empirical_data_file = sys.argv[2]
	except:
		print("Please give empirical data to fit")
		#exit()
		#empirical_data_file = "/home/kprovost/nas5/train_nn/DEM-IBD-AGE_empiricalcomparison_SCUTULATUS.COMBO.STATS"
		empirical_data_file = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/cardcard16.fasta.FULLSTATS.COMBO.STATS"
		## atrox, cardinalis, catenifer, getula, scutulatus
	try:
		X_train_file = sys.argv[3]
	except:
		print("Please give trained data to scale on")
		#exit()
		X_train_file = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/X_train.**NN_expand_1604531680.temp"
	clf = load(trained_nn_file)
	#clf2 = clf.estimators_[0].best_estimator_.coefs_
	empirical_data = np.loadtxt(empirical_data_file)
	X_train = np.loadtxt(X_train_file,delimiter=",")
	X_scaled_file = X_train_file.replace("train.","train_scaled.")
	X_scaled = np.loadtxt(X_scaled_file,delimiter=",")
	Y_train_file = X_train_file.replace("/X","/Y")
	Y_train = np.loadtxt(Y_train_file,delimiter=",",dtype="str")
	cv = ShuffleSplit(n_splits=100, test_size=0.2, random_state=0)
	#train_sizes, train_scores, test_scores, fit_times, _ =learning_curve(clf, X_scaled, Y_train, cv=cv,return_times=True)
	X_empirical_scaled = featureScalingEmpirical(X_train,empirical_data)
	outfileString_base = trained_nn_file+"_"+os.path.basename(X_train_file)+"_"+os.path.basename(empirical_data_file)
	empirical_predicted = clf.predict(X_empirical_scaled)
	empirical_predicted_prob = clf.predict_proba(X_empirical_scaled)
	demog = empirical_predicted_prob[0]
	np.savetxt(outfileString_base+"_DEMOG_PROB.txt", demog, delimiter="\t",fmt='%s',header=str(clf.classes_[0]))
	ibd = empirical_predicted_prob[1]
	np.savetxt(outfileString_base+"_IBD_PROB.txt", ibd, delimiter="\t",fmt='%s',header=str(clf.classes_[1]))
	age = empirical_predicted_prob[2]
	np.savetxt(outfileString_base+"_AGE_PROB.txt", age, delimiter="\t",fmt='%s',header=str(clf.classes_[2]))	
	np.savetxt(outfileString_base+"_PREDICTED.txt", empirical_predicted, delimiter="\t",fmt='%s')

main()

#if __name__ == "__main__":
#	main()
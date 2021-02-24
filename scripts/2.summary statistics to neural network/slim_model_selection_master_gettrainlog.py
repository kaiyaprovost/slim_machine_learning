#!/usr/bin/python -u
# coding: utf-8

# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW//NN_Results_Empirical/"
# trained="X_train.**NN_expand_1604531680.temp"
# for i in {2..3}; do
# echo $i
# python3 -u "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/2.summary statistics to neural network/slim_model_selection_master_gettrainlog.py" $trained > logtile_${trained}_${i}.out.txt 2>&1
# done

# cd "/home/kprovost/nas5/train_nn/FINAL_VERSION/"
# trained="X_train.**NN_expand_1604531680.temp"
# for i in {1..2}; do
# echo $i
# python3 -u "/home/kprovost/nas5/train_nn/slim_model_selection_master_gettrainlog.py" $trained > /home/kprovost/nas5/train_nn/EMPIRICAL_LOG_RETRAIN/logfile_${trained}_${i}.out.txt 2>&1
# done 

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

from sklearn.model_selection import GridSearchCV ## 
from sklearn.model_selection import ParameterGrid ## 
from sklearn.multioutput import MultiOutputClassifier ## 
from sklearn.neural_network import MLPClassifier ## 
from sklearn.metrics import classification_report ## 


def main():
	try:
		X_train_file = sys.argv[1]
	except:
		print("Please give trained data to scale on")
		#exit()
		X_train_file = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/X_train.**NN_expand_1604531680.temp"
	X_train = np.loadtxt(X_train_file,delimiter=",")
	X_scaled_file = X_train_file.replace("train.","train_scaled.")
	X_scaled = np.loadtxt(X_scaled_file,delimiter=",")
	Y_train_file = X_train_file.replace("/X","/Y")
	Y_train = np.loadtxt(Y_train_file,delimiter=",",dtype="str")

	estimator=MLPClassifier(activation='relu',alpha=0.0001,batch_size='auto',beta_1=0.9,beta_2=0.999,early_stopping=False,
		epsilon=1e-08,hidden_layer_sizes=(100,),learning_rate='constant',learning_rate_init=0.001,max_fun=15000,max_iter=1000,momentum=0.9,
		n_iter_no_change=10,verbose=True)

	estimator=MultiOutputClassifier(estimator)
	estimator.fit(X_scaled,Y_train)

main()

#if __name__ == "__main__":
#	main()
#!/usr/bin/python -u
# coding: utf-8

import numpy as np
import pandas as pd
from joblib import dump, load
import time

import slim_model_selection_data_processing as sms_pre
import slim_model_selection_build_NN as sms_bnn
import slim_model_selection_summarize_NN as sms_snn

# source activate py36
# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/NONCORRELATED"
# programPath="/Users/kprovost/Documents/Github/slim_machine_learning/scripts/2.summary statistics to neural network/slim_model_selection_master_with_export.py"
# echo "COMBO ALL" 2>&1 > log_pop-demography-ibd-age-colfix.temp 
# ls DEM-IBD-AGE_*1000*.COMBO.STATS | xargs -t python3 "${programPath}" >> log_pop-demography-ibd-age-colfix_smallshort.temp 2>&1

def main():

	## consider making the outputs of each of these modules a class instead?
	
	current_time=int(time.time())
	
	print("\n------------------------------------\nSTART DATA PROCESSING MODULE\n------------------------------------")

	filenames = sms_pre.readInFilenames()
	
	[X,Y,numcategories,classOrderLs,outputString] = sms_pre.inputDataFiles(filenames)
	
	suffix1="*"
	suffix2="*NNexp"
	suffix=suffix1+suffix2
	
	[X_train,X_test,X_validate,X_cv_list,Y_train,Y_test,Y_validate,Y_cv_list] = sms_pre.manualTrainTestSplit(X,Y,N=10)
	
	sms_pre.outputXYdatasets(dataset=X_train,writename="X_train."+str(suffix)+"_"+str(current_time)+".temp")
	sms_pre.outputXYdatasets(dataset=X_validate,writename="X_validate."+str(suffix)+"_"+str(current_time)+".temp")
	sms_pre.outputXYdatasets(dataset=X_test,writename="X_test."+str(suffix)+"_"+str(current_time)+".temp")
	sms_pre.outputXYdatasets(dataset=Y_train,writename="Y_train."+str(suffix)+"_"+str(current_time)+".temp")
	sms_pre.outputXYdatasets(dataset=Y_validate,writename="Y_validate."+str(suffix)+"_"+str(current_time)+".temp")
	sms_pre.outputXYdatasets(dataset=Y_test,writename="Y_test."+str(suffix)+"_"+str(current_time)+".temp")
	
	print("\nShapes of features:")
	print("Full dataset X, Y")
	print(X.shape,Y.shape)
	print("\nTrain X, Y")
	print(X_train.shape, Y_train.shape)
	print("CV X, Y")
	print(X_validate.shape, Y_validate.shape)
	print("Test X, Y")
	print(X_test.shape, Y_test.shape)
	
	[X_train_scaled,X_test_scaled,X_validate_scaled] = sms_pre.featureScaling(X_train,X_test,X_validate)
	
	sms_pre.outputXYdatasets(dataset=X_train_scaled,writename="X_train_scaled."+str(suffix)+"_"+str(current_time)+".temp")
	sms_pre.outputXYdatasets(dataset=X_validate_scaled,writename="X_validate_scaled."+str(suffix)+"_"+str(current_time)+".temp")
	sms_pre.outputXYdatasets(dataset=X_test_scaled,writename="X_test_scaled."+str(suffix)+"_"+str(current_time)+".temp")
	
	print("\n---------------------------------------\nEND DATA PROCESSING MODULE\n---------------------------------------")
	
	print("\n---------------------------------------\nSTART BUILDING NN MODULE\n---------------------------------------")

	X_train = X_train_scaled
	X_test = X_test_scaled
	X_validate = X_validate_scaled
	
	## BELOW FOR NEURAL NET

	print("\nOptimizing Neural Network")
	param_grid = [
		{
		'alpha': [0,0.01,0.1,1,10], ## for 1000k, 1 and 10 are trash, and rest roughly the same. For 6k and 21k, 10 is trash, rest roughly same. 
		#'alpha': [0,1],
		#'estimator__alpha': [1], ## adding this is helpful when you have a clf inside a clf
		'activation': ['logistic'],
		#'activation': ['identity','logistic','relu','tanh'],
		#'hidden_layer_sizes': [(100),(100,100),(100,100,100)],
		'hidden_layer_sizes': [(5,5,5),(25,25,25),(100,100,100),(5,5),(25,25),(100,100),(5),(25),(100)],
		#'hidden_layer_sizes': [(25),(25,25),(25,25,25)],
		#'hidden_layer_sizes': [(5),(5,5),(5,5,5)],
		#'solver': ['sgd','lbfgs','adam'],
		'solver': ['lbfgs'],
		#'learning_rate': ['adaptive','constant','invscaling'],
		'learning_rate': ['adaptive'],
		#'learning_rate_init': [1e-1,1e-4,1e-7],
		'learning_rate_init': [1e-4],
		'shuffle': [True],
		'verbose': [False],
		'early_stopping': [True],
		'max_iter': [10000], ## this was 1000 before and decidedly not high enough 
		'warm_start': [True]
		}
	]

	#clf = sms_bnn.makeSimpleNN(X_train,Y_train)
	clf = sms_bnn.makeGridCVNN(X_train=X_train,Y_train=Y_train,param_grid=param_grid,X_test=X_validate,Y_test=Y_validate,n_jobs=64)
	
	train_predicted = sms_bnn.summarizePredictScores(X_train,Y_train,clf,label="Training")
	validate_predicted = sms_bnn.summarizePredictScores(X_validate,Y_validate,clf,label="Validation")
	test_predicted = sms_bnn.summarizePredictScores(X_test,Y_test,clf,label="Test")
	
	## export clf object for later usage
	
	dump(clf,outputString+"trained_neural_network_"+str(suffix)+"_"+str(current_time)+".joblib")
	
	print("\n---------------------------------------\nEND BUILDING NN MODULE\n---------------------------------------")

	print("\n---------------------------------------\nSTART SUMMARIZE NN MODULE\n---------------------------------------")

	
	#[test_predicted_merged,Y_test_merged] = sms_snn.mergeMultioutputToSingleoutput2(test_predicted,Y_test)
	#[train_predicted_merged,Y_train_merged] = sms_snn.mergeMultioutputToSingleoutput2(train_predicted,Y_train)
	#[validate_predicted_merged,Y_validate_merged] = sms_snn.mergeMultioutputToSingleoutput2(validate_predicted,Y_validate)

	test_predicted_merged = sms_snn.mergeMultioutputToSingleoutput(test_predicted)
	train_predicted_merged = sms_snn.mergeMultioutputToSingleoutput(train_predicted)
	validate_predicted_merged = sms_snn.mergeMultioutputToSingleoutput(validate_predicted)
	Y_test_merged = sms_snn.mergeMultioutputToSingleoutput(Y_test)
	Y_train_merged = sms_snn.mergeMultioutputToSingleoutput(Y_train)
	Y_validate_merged = sms_snn.mergeMultioutputToSingleoutput(Y_validate)
	
	sms_snn.summarizeMergedScores(Y_train_merged,train_predicted_merged,numcategories,Y_train,train_predicted,label="Training")
	sms_snn.summarizeMergedScores(Y_validate_merged,validate_predicted_merged,numcategories,Y_validate,validate_predicted,label="Validation")
	sms_snn.summarizeMergedScores(Y_test_merged,test_predicted_merged,numcategories,Y_test,test_predicted,label="Test")
	
	## SOURCE OF THIS: http://sdsawtelle.github.io/blog/output/week4-andrew-ng-machine-learning-with-python.html
	
	print("\nOutput the following counts for train:")
	print(classOrderLs)
	counts_train = sms_snn.generateCountsForConfusionMatrix(classOrderLs,Y_train_merged,train_predicted_merged)
	print(counts_train)
	df_train = pd.DataFrame(counts_train)

	print("\nOutput the following counts for validate:")
	print(classOrderLs)
	counts_validate = sms_snn.generateCountsForConfusionMatrix(classOrderLs,Y_validate_merged,validate_predicted_merged)
	print(counts_validate)
	df_validate = pd.DataFrame(counts_validate)

	print("\nOutput the following counts for test:")
	print(classOrderLs)
	counts = sms_snn.generateCountsForConfusionMatrix(classOrderLs,Y_test_merged,test_predicted_merged)
	print(counts)
	df = pd.DataFrame(counts)

	[true_positives,false_positives,true_negatives,false_negatives] = sms_snn.trueFalsePosNegValues(classOrderLs,counts)

	print("\n",classOrderLs)
	print("TP:",true_positives)
	print("FP:",false_positives)
	print("TN:",true_negatives)
	print("FN:",false_negatives)
	
	sms_snn.outputCountsToFile(classOrderLs,df,str(suffix)+"_"+str(current_time)+"_"+outputString)
	
	## output the weights 
	
	print("\n---------------------------------------\nEND SUMMARIZE NN MODULE\n---------------------------------------")


	
if __name__ == "__main__":
	main()
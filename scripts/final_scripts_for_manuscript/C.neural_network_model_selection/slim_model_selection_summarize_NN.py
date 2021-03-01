#!/usr/bin/python -u
# coding: utf-8

import numpy as np ## 
import pandas as pd ## 

from sklearn.metrics import accuracy_score ## 
from sklearn.metrics import precision_recall_fscore_support ## 
from sklearn.model_selection import GridSearchCV ## 
from sklearn.model_selection import ParameterGrid ## 
from sklearn.multioutput import MultiOutputClassifier ## 
from sklearn.neural_network import MLPClassifier ## 

import slim_model_selection_data_processing as sms_pre
import slim_model_selection_build_NN as sms_bnn


def mergeMultioutputToSingleoutput2(data_predicted,Y_data):
	data_predicted_merged = []
	Y_data_merged = []
	for i in range((len(data_predicted))):

		joined = "-".join(data_predicted[i,])
		data_predicted_merged += [joined]
	
		joined2 = "-".join(Y_data[i,])
		Y_data_merged += [joined2]
	
	return([data_predicted_merged,Y_data_merged])
	
def mergeMultioutputToSingleoutput(data):
	data_merged = []
	for i in range((len(data))):

		joined = "-".join(data[i,])
		data_merged += [joined]
	
	return(data_merged)

def summarizeMergedScores(Y_data_merged,data_predicted_merged,numcategories,Y_data,data_predicted,label=""):

	print("\n***"+label+" dataset***")
	[precision,recall,fscore,support] = precision_recall_fscore_support(Y_data_merged, data_predicted_merged, average='weighted')
	print("Precision\t",round(precision,4),"\nRecall\t",round(recall,4),"\nFscore\t",round(fscore,4))

	if numcategories > 1:
		print("\n"+label+" accuracy per axis looked at")
		for c in range(numcategories):
			tot_bins = "-".join(np.sort(np.unique(Y_data[:,c])))
			#print(tot_bins)
			#print(Y_data[0,c],end="\t")
			Y_col = Y_data[:,c]
			data_col = data_predicted[:,c]
			print(tot_bins,"\t",round(accuracy_score(Y_col,data_col),4))
			
	return()

def generateCountsForConfusionMatrix(classOrderLs,Y_test_merged,test_predicted_merged):
	counts = np.zeros((len(classOrderLs),len(classOrderLs)))
	
	d = dict([(y,x) for x,y in enumerate(sorted(set(classOrderLs)))])

	for i in range(len(Y_test_merged)):
		Y_test_key = (Y_test_merged[i])
		test_predicted_key = (test_predicted_merged[i])
		Y_test_value = d.get(Y_test_key)
		test_predicted_value = d.get(test_predicted_key)
	
		counts[Y_test_value][test_predicted_value] += 1
		
	return(counts)

## IMPLEMENT TRUE AND FALSE POSITIVES

def trueFalsePosNegValues(classOrderLs,counts):

	true_positives = []
	false_positives = []
	true_negatives = []
	false_negatives = []

	TP_temp = 0
	FP_temp = 0
	TN_temp = 0
	FN_temp = 0

	for lab in range(len(classOrderLs)):
		for row in range(int(counts.shape[0])):
			for col in range(int(counts.shape[1])):
				#print(lab,row,col)
				value = counts[int(row)][int(col)]
				if row == lab and col == lab: 
					TP_temp += value
				elif row == lab and col != lab:
					FN_temp += value
				elif row != lab and col == lab:
					FP_temp += value
				elif row != lab and col != lab:
					TN_temp += value
	
		sum_cat = TP_temp+FP_temp+TN_temp+FN_temp
	
		true_positives += [round(TP_temp/sum_cat,2)]
		false_positives += [round(FP_temp/sum_cat,2)]
		true_negatives += [round(TN_temp/sum_cat,2)]
		false_negatives += [round(FN_temp/sum_cat,2)]
	
		TP_temp = 0
		FP_temp = 0
		TN_temp = 0
		FN_temp = 0 
		
	return([true_positives,false_positives,true_negatives,false_negatives])

def outputCountsToFile(classOrderLs,df,outputString):

	index_dict = {}
	header_list = classOrderLs
	for j in range(len(classOrderLs)):
		label = classOrderLs[j]
		index_dict[j] = label+"-true"
		header_list[j] = label+"-pred"

	df.rename(index=index_dict,inplace=True)

	csv_name = "Confusion_Matrix_"+outputString

	print("\nWriting to "+csv_name+".csv")
	df.to_csv(csv_name+".csv",header=header_list)

def main():
	
	print("\n------------------------------------\nSTART DATA PROCESSING MODULE\n------------------------------------\n")

	filenames = sms_pre.readInFilenames()
	
	[X,Y,numcategories,classOrderLs,outputString] = sms_pre.inputDataFiles(filenames)
	
	[X_train,X_test,X_validate,X_cv_list,Y_train,Y_test,Y_validate,Y_cv_list] = sms_pre.manualTrainTestSplit(X,Y,N=10)
	
	sms_pre.outputXYdatasets(dataset=X_train,writename="X_train.temp")
	sms_pre.outputXYdatasets(dataset=X_validate,writename="X_validate.temp")
	sms_pre.outputXYdatasets(dataset=X_test,writename="X_test.temp")
	sms_pre.outputXYdatasets(dataset=Y_train,writename="Y_train.temp")
	sms_pre.outputXYdatasets(dataset=Y_validate,writename="Y_validate.temp")
	sms_pre.outputXYdatasets(dataset=Y_test,writename="Y_test.temp")
	
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
	
	sms_pre.outputXYdatasets(dataset=X_train_scaled,writename="X_train_scaled.temp")
	sms_pre.outputXYdatasets(dataset=X_validate_scaled,writename="X_validate_scaled.temp")
	sms_pre.outputXYdatasets(dataset=X_test_scaled,writename="X_test_scaled.temp")
	
	print("\n---------------------------------------\nEND DATA PROCESSING MODULE\n---------------------------------------\n")
	
	print("\n---------------------------------------\nSTART BUILDING NN MODULE\n---------------------------------------\n")

	X_train = X_train_scaled
	X_test = X_test_scaled
	X_validate = X_validate_scaled
	
	## BELOW FOR NEURAL NET

	print("\nOptimizing Neural Network")
	param_grid = [
		{
		'alpha': [0,0.01,0.1,1,10],
		#'alpha': [0.1,0.01],
		'activation': ['logistic'],
		#'activation': ['identity','logistic','relu','tanh'],
		'hidden_layer_sizes': [(100,100,100)],
		#'hidden_layer_sizes': [(5,5,5),(25,25,25),(100,100,100)],
		#'solver': ['sgd','lbfgs','adam'],
		'solver': ['lbfgs'],
		#'learning_rate': ['adaptive','constant','invscaling'],
		'learning_rate': ['adaptive'],
		#'learning_rate_init': [1e-1,1e-4,1e-7],
		'learning_rate_init': [1e-4],
		'shuffle': [True],
		'verbose': [True],
		'early_stopping': [False],
		'max_iter': [1000],
		'warm_start': [True]
		}
	]

	#clf = sms_bnn.makeSimpleNN(X_train,Y_train)
	clf = sms_bnn.makeGridCVNN(X_train,Y_train,param_grid)
	
	train_predicted = sms_bnn.summarizePredictScores(X_train,Y_train,clf,label="Training")
	validate_predicted = sms_bnn.summarizePredictScores(X_validate,Y_validate,clf,label="Validation")
	test_predicted = sms_bnn.summarizePredictScores(X_test,Y_test,clf,label="Test")
	
	print("\n---------------------------------------\nEND BUILDING NN MODULE\n---------------------------------------\n")

	print("\n---------------------------------------\nSTART SUMMARIZE NN ALONE\n---------------------------------------\n")

	
	test_predicted_merged = mergeMultioutputToSingleoutput(test_predicted)
	train_predicted_merged = mergeMultioutputToSingleoutput(train_predicted)
	validate_predicted_merged = mergeMultioutputToSingleoutput(validate_predicted)
	Y_test_merged = mergeMultioutputToSingleoutput(Y_test)
	Y_train_merged = mergeMultioutputToSingleoutput(Y_train)
	Y_validate_merged = mergeMultioutputToSingleoutput(Y_validate)
	
	summarizeMergedScores(Y_train_merged,train_predicted_merged,numcategories,Y_train,train_predicted,label="Training")
	summarizeMergedScores(Y_validate_merged,validate_predicted_merged,numcategories,Y_validate,validate_predicted,label="Validation")
	summarizeMergedScores(Y_test_merged,test_predicted_merged,numcategories,Y_test,test_predicted,label="Test")
	
	## SOURCE OF THIS: http://sdsawtelle.github.io/blog/output/week4-andrew-ng-machine-learning-with-python.html
	
	print("\nOutput the following counts for train:")
	print(classOrderLs)
	counts_train = generateCountsForConfusionMatrix(classOrderLs,Y_train_merged,train_predicted_merged)
	print(counts_train)
	df_train = pd.DataFrame(counts_train)

	print("\nOutput the following counts for validate:")
	print(classOrderLs)
	counts_validate = generateCountsForConfusionMatrix(classOrderLs,Y_validate_merged,validate_predicted_merged)
	print(counts_validate)
	df_validate = pd.DataFrame(counts_validate)

	print("\nOutput the following counts for test:")
	print(classOrderLs)
	counts = generateCountsForConfusionMatrix(classOrderLs,Y_test_merged,test_predicted_merged)
	print(counts)
	df = pd.DataFrame(counts)

	[true_positives,false_positives,true_negatives,false_negatives] = trueFalsePosNegValues(classOrderLs,counts)

	print("\n",classOrderLs)
	print("TP:",true_positives)
	print("FP:",false_positives)
	print("TN:",true_negatives)
	print("FN:",false_negatives)
	
	outputCountsToFile(classOrderLs,df,outputString)
	
	print("\n---------------------------------------\nEND SUMMARIZE NN ALONE\n---------------------------------------\n")

	
if __name__ == "__main__":
	main()
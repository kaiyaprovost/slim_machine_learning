#!/usr/bin/python -u
# coding: utf-8

import numpy as np ## 
import pandas as pd ## 

from sklearn.model_selection import GridSearchCV ## 
from sklearn.model_selection import ParameterGrid ## 
from sklearn.multioutput import MultiOutputClassifier ## 
from sklearn.neural_network import MLPClassifier ## 

import slim_model_selection_data_processing as sms_pre

def makeSimpleNN(X_train,Y_train,item=None):
	print("\tSimple NN")
	clf = MLPClassifier() 
	
	if item != None:
		clf.set_params(**item)
	clf = MultiOutputClassifier(clf,n_jobs=-1) ## if you turn this off it splits the data 50/50 
	clf = clf.fit(X_train, Y_train) 
	#print(clf.get_params())
	
	return(clf)

def makeGridCVNN(X_train,Y_train,param_grid):
	print("\tGrid Search NN")
	print("Fitting the following parameters:")
	print(param_grid,end="\n")
	clf = MLPClassifier() ## bring to developers? or check with RF? 
	## see if you can use your data in R or matlab 
	
	clf = GridSearchCV(clf, param_grid,refit=True) ## tweak this 
	clf = MultiOutputClassifier(clf,n_jobs=-1) ## if you turn this off it splits the data 50/50 
	clf = clf.fit(X_train, Y_train) 
	print("\nBest parameters:")
	print(clf.get_params())
	
	return(clf)
	
def manualSingleValidationParameters(param_grid,X_train,X_validate,Y_train,Y_validate):
	print("\tNN with manual single validation")
	print("Fitting the following parameters:")
	print(param_grid,end="\n")

	param_grid_list = list(ParameterGrid(param_grid))

	print("\nGrid searching, single validation")

	## generate obj called classifier_score
	classifier_scores = []
	## generate obj called best_classifier_score
	best_classifier_score = 0
	## generate obj called best_parameters
	best_parameters = []
	
	classifier_trains = []
	best_classifier_train = 0
	best_parameters_train = []
	
	classifier_harmonic = []
	best_classifier_harmonic = 0
	best_parameters_harmonic = []

	#for item in param_grid_list:
	print("Total param combos to check:",len(param_grid_list))
	for index in range(len(param_grid_list)):
		print("Param",str(index+1),end=" ")
		item = param_grid_list[index]
	
		X_test_cv = X_train
		X_train_cv = X_validate
	
		Y_test_cv = Y_train
		Y_train_cv = Y_validate
	
		#print(X_train_cv.shape)
		#print(X_test_cv.shape)
		#print(Y_train_cv.shape)
		#print(Y_test_cv.shape)

		## set the item to be the parameters for the MLP Classifier and MultiOutputClassifier
		#print(item)
		clf = MLPClassifier()
		clf.set_params(**item)

		## train the classifier normally on X_train_cv
		clf = MultiOutputClassifier(clf,n_jobs=-1)
		#print(clf.get_params(),end="\n###\n")
		clf = clf.fit(X_train_cv, Y_train_cv) 

		## evaluate the classifier on X_test_cv
		train_pred_cv = clf.predict(X_train_cv)
		test_pred_cv = clf.predict(X_test_cv)
	
		## this really needs to be f-score not accuracy its finding local minima
		## however the fscore function being used here does not account for multilabel class
		score = clf.score(X_test_cv, Y_test_cv)
		train = clf.score(X_train_cv,Y_train_cv)
		harmonic = 2*((score*train)/(score+train))
	
		## store the score in temp_classifier_scores
		
		#print(train,score)

		## store the average score of the classifier in classifier_score
	
		classifier_scores += [round(score,4)]
		classifier_trains += [round(train,4)]
		classifier_harmonic += [round(harmonic,4)]
	
		## evaluate the current and best score to see which one to keep
		if best_classifier_score < score:
			best_classifier_score = score
			best_parameters = [item]
			best_model = clf
			
		if best_classifier_train < train:
			best_classifier_train = train
			best_parameters_train = [item]
			best_model_train = clf
			
		if best_classifier_harmonic < harmonic:
			best_classifier_harmonic = harmonic
			best_parameters_harmonic = [item]
			best_model_harmonic = clf
	
	print("\n---\nAll average train scores:")
	print(classifier_trains)

	print("\nBest Train Scores belong to:")
	print(best_classifier_train)
	print(best_parameters_train)
	
	print("\n---\nAll average test scores:")
	print(classifier_scores)

	print("\nBest Test Scores belong to:")
	print(best_classifier_score)
	print(best_parameters)

	print("\n---\nAll harmonic scores:")
	print(classifier_harmonic)
	
	print("\nBest Harmonic Scores belong to:")
	print(best_classifier_harmonic)
	print(best_parameters_harmonic)

	print("\nSetting clf to be best test model")
	clf = best_model
	
	return(clf)

def manualCrossValidationParameters(param_grid,X_cv_list,Y_cv_list):
	print("\tNN with manual cross validation")
	print("Fitting the following parameters:")
	print(param_grid,end="\n")

	param_grid_list = list(ParameterGrid(param_grid))

	print("\nGrid searching, cross-validation")

	## generate obj called classifier_score
	classifier_scores = []
	## generate obj called best_classifier_score
	best_classifier_score = 0
	## generate obj called best_parameters
	best_parameters = []
	
	classifier_trains = []
	best_classifier_train = 0
	best_parameters_train = []
	
	classifier_harmonic = []
	best_classifier_harmonic = 0
	best_parameters_harmonic = []

	print("Total param combos to check:",len(param_grid_list))
	for index in range(len(param_grid_list)):
		print("\nParam",str(index+1),end=" ")
		item = param_grid_list[index]

		temp_classifier_scores = []
		temp_classifier_trains = []
		temp_classifier_harmonic = []

		print("")

		for i in range(len(X_cv_list)):    
			print("\tRep",i, end = " ")
			## may need to do this k-fold 
		
			X_test_cv = X_cv_list[i]
			X_train_cv = np.concatenate(X_cv_list[:i] + X_cv_list[i+1:],axis=0)
		
			Y_test_cv = Y_cv_list[i]
			Y_train_cv = np.concatenate(Y_cv_list[:i] + Y_cv_list[i+1:],axis=0)
		
			#print(X_train_cv.shape)
			#print(X_test_cv.shape)
			#print(Y_train_cv.shape)
			#print(Y_test_cv.shape)

			## set the item to be the parameters for the MLP Classifier and MultiOutputClassifier
			#print(item)
			clf = MLPClassifier()
			clf.set_params(**item)
	
			## train the classifier normally on X_train_cv
			clf = MultiOutputClassifier(clf,n_jobs=-1)
			#print(clf.get_params(),end="\n###\n")
			clf = clf.fit(X_train_cv, Y_train_cv) 
	
			## evaluate the classifier on X_test_cv
			train_pred_cv = clf.predict(X_train_cv)
			test_pred_cv = clf.predict(X_test_cv)
		
			## this really needs to be f-score not accuracy its finding local minima
			## however the fscore function being used here does not account for multilabel class
			score = clf.score(X_test_cv, Y_test_cv)
			train = clf.score(X_train_cv,Y_train_cv)
			harmonic = 2*((score*train)/(score+train))
		
			## store the score in temp_classifier_scores
			temp_classifier_scores += [score]
			temp_classifier_trains += [train]
			temp_classifier_harmonic += [harmonic]

		## calculate the average score 
	
		#print("\n\tTEMP TRAIN")
		#print("\t",temp_classifier_trains)
		#print("\tTEMP TEST")
		#print("\t",temp_classifier_scores)
		average_score = np.mean(temp_classifier_scores)
		average_train = np.mean(temp_classifier_trains)
		average_harmonic = np.mean(temp_classifier_harmonic)
	
		## store the average score of the classifier in classifier_score
	
		classifier_scores += [round(average_score,4)]
		classifier_trains += [round(average_train,4)]
		classifier_harmonic += [round(average_harmonic,4)]
	
		## evaluate the current and best score to see which one to keep
		if best_classifier_score < average_score:
			best_classifier_score = average_score
			best_parameters = [item]
			best_model = clf
		if best_classifier_train < average_train:
			best_classifier_train = average_train
			best_parameters_train = [item]
			best_model_train = clf			
		if best_classifier_harmonic < average_harmonic:
			best_classifier_harmonic = average_harmonic
			best_parameters_harmonic = [item]
			best_model_harmonic = clf	
			
	print("\nAll average train scores:")
	print(classifier_trains)

	print("\nBest Train Scores belong to:")
	print(best_classifier_train)
	print(best_parameters_train)
	
	print("\nAll average test scores:")
	print(classifier_scores)

	print("\nBest Test Scores belong to:")
	print(best_classifier_score)
	print(best_parameters)
	
	print("\nAll average harmonic scores:")
	print(classifier_harmonic)

	print("\nBest harmonic Scores belong to:")
	print(best_classifier_harmonic)
	print(best_parameters_harmonic)

	print("\nSetting clf to be best test model")
	clf = best_model
	
	return(clf)

## REDO THIS AS ONE DATASET NOT TWO

def summarizePredictScores(X_data,Y_data,clf,label):
	data_predicted = clf.predict(X_data) ## "Call fit before using this method"
	print("\n"+label+" set score: %f" % clf.score(X_data, Y_data))
	return(data_predicted)

def summarizePredictScores2(X_train,X_test,X_validate,clf):
	
	## predict the values of X_train and X_test
	train_predicted = clf.predict(X_train) ## "Call fit before using this method"
	test_predicted = clf.predict(X_test)
	validate_predicted = clf.predict(X_validate)

	print("\nTraining set score: %f" % clf.score(X_train, Y_train))
	print("CV set score: %f" % clf.score(X_validate, Y_validate))
	print("Test set score: %f" % clf.score(X_test, Y_test))
	
	return([train_predicted,test_predicted,validate_predicted])

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
	
	print("\n---------------------------------------\nSTART BUILDING NN ALONE\n---------------------------------------\n")

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
		'verbose': [False],
		'early_stopping': [False],
		'max_iter': [1000],
		'warm_start': [True]
		}
	]

	#clf = makeSimpleNN(X_train,Y_train)
	clf = makeGridCVNN(X_train,Y_train,param_grid)
	
	train_predicted = summarizePredictScores(X_train,Y_train,clf,label="Training")
	validate_predicted = summarizePredictScores(X_validate,Y_validate,clf,label="Validation")
	test_predicted = summarizePredictScores(X_test,Y_test,clf,label="Test")
	
	print("\n---------------------------------------\nEND BUILDING NN ALONE\n---------------------------------------\n")

	
if __name__ == "__main__":
	main()
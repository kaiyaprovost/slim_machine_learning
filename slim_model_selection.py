#!/usr/bin/python -u
# coding: utf-8

import numpy as np
import pandas as pd
import sys
import glob
import os
import matplotlib.pyplot as plt
import snips as snp

#from sklearn import svm
#from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_recall_fscore_support
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import ParameterGrid
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import normalize
from sklearn.preprocessing import StandardScaler  



## TODO: figure out why CV vd TEST not working
## TODO: functionize repetitive code 

# python slim_model_selection.py "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/test_results/p1_ms_exponential_sumstats_small.stats" "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/test_results/p1_ms_decline_sumstats_small.stats"
# python slim_model_selection.py ""

## read in arguments
try:
	filenames = sys.argv[1:]
	print("\nRead the following filenames:")
	print(filenames)
except:
	print("\nNo filenames given, quitting")
	sys.exit()
if filenames == []:
	print("\nInvalid filenames given, quitting")
	sys.exit()

## suck up the data
print("\nReading data")

def inputDataFiles(filenames):

	tempname = filenames[0].split("_")[-1].split(".")[0]
	type = filenames[0].split("_")[-1].split(".")[1] ## popgenome or sumstats or combo
	outputString = filenames[0].split("_")[0]+"_"+type

	numcategories = len(tempname.split("-"))
	#print(numcategories)


	if type == "POPGENOME":
		X = np.array([]).reshape(0,36)
	elif type == "SUMSTATS":
		X = np.array([]).reshape(0,5)
	elif type == "COMBO":
		X = np.array([]).reshape(0,41)
	else:
		X = np.array([]).reshape(0,0)
	#X = np.array([])


	y = []
	classOrderLs = []
	outputString = ""
	for i in range(len(filenames)):
		## read the file 
		file = filenames[i]
	
		## count the lines of the file
	
		if type == "SUMSTATS":
			Xi = np.loadtxt(file,usecols=(1,3,5,7,9))
		else:
			Xi = np.loadtxt(file)
		#
		## concatenate the file to the main dataframe
		try:
			X = np.concatenate((X,Xi))
		except:
			X.reshape(X.shape[0],Xi.shape[1])
			X = np.concatenate((X,Xi))
	
		## add the name to classOrderLs
		name = file.split("_")[-1].split(".")[0]
		splitname = name.split("-")
	
		classOrderLs = classOrderLs + [name]
		outputString = file.split("_")[0]
	
		y = y + ([splitname]*len(Xi))
	Y = np.array(y)
	
	return([X,Y,numcategories,classOrderLs,outputString])

[X,Y,numcategories,classOrderLs,outputString] = inputDataFiles(filenames)

## split data into training and test 
print("\nCreating test, validate, and training data")

def manualTrainTestSplit(X,Y,N=10):

	## create 10 sets, one test, one validation, eight train/cv
	
	## NOTE: REMOVED STRATIFY
	
	X_train, X_test, Y_train, Y_test = train_test_split(X,Y,test_size=1/N,shuffle=True)
	X_train, X_validate, Y_train, Y_validate = train_test_split(X_train,Y_train,test_size=1/(N-1),shuffle=True)

	X_cv_list = []
	Y_cv_list = []

	X_train_cv = X_train
	Y_train_cv = Y_train

	print("Splitting training data into",str(N-2),"Train/cv sets")
	for i in list(range((N-2),1,-1)):
		if i == 2:
			X_cvA, X_cvB, Y_cvA, Y_cvB = train_test_split(X_train_cv,Y_train_cv,test_size=(1/i),shuffle=True)
			X_cv_list.append(X_cvA)
			Y_cv_list.append(Y_cvA)
			X_cv_list.append(X_cvB)
			Y_cv_list.append(Y_cvB)
		
		else:
			X_train_cv, X_cvA, Y_train_cv, Y_cvA = train_test_split(X_train_cv,Y_train_cv,test_size=(1/i),shuffle=True,stratify=Y_train_cv)
			X_cv_list.append(X_cvA)
			Y_cv_list.append(Y_cvA)	

	## this will result in 80% train, 10% CV, 10% test and even splits of the Y parameters 
	
	toReturn = [X_train,X_test,X_validate,X_cv_list,Y_train,Y_test,Y_validate,Y_cv_list]
	return(toReturn)

[X_train,X_test,X_validate,X_cv_list,Y_train,Y_test,Y_validate,Y_cv_list] = manualTrainTestSplit(X,Y,N=10)

## TRY RE-OUTPUT AND RE-INPUT THE DATA
pd.DataFrame(X_train).to_csv("X_train.temp",header=False,index=False)
pd.DataFrame(X_validate).to_csv("X_validate.temp",header=False,index=False)
pd.DataFrame(X_test).to_csv("X_test.temp",header=False,index=False)

pd.DataFrame(Y_train).to_csv("Y_train.temp",header=False,index=False)
pd.DataFrame(Y_validate).to_csv("Y_validate.temp",header=False,index=False)
pd.DataFrame(Y_test).to_csv("Y_test.temp",header=False,index=False)

X_train = np.loadtxt("X_train.temp",delimiter=",")
X_validate = np.loadtxt("X_test.temp",delimiter=",")
X_test = np.loadtxt("X_validate.temp",delimiter=",")

Y_train = np.loadtxt("Y_train.temp",delimiter=",",dtype=str).reshape(Y_train.shape[0],numcategories)
Y_validate = np.loadtxt("Y_test.temp",delimiter=",",dtype=str).reshape(Y_validate.shape[0],numcategories)
Y_test = np.loadtxt("Y_validate.temp",delimiter=",",dtype=str).reshape(Y_test.shape[0],numcategories)

## note: you don't need to separate the CV out! the GridSearchCV algorithm does this automatically 

print("\nShapes of features:")
print("Full dataset X, Y")
print(X.shape,Y.shape)
print("\nTrain X, Y")
print(X_train.shape, Y_train.shape)

#X_train_cv = X_train
#Y_train_cv = Y_train

print("CV X, Y")
print(X_validate.shape, Y_validate.shape)
print("Test X, Y")
print(X_test.shape, Y_test.shape)

## implement feature scaling
scaler = StandardScaler()  
# Don't cheat - fit only on training data
scaler.fit(X_train)  
X_train = scaler.transform(X_train)  
#X_train_cv = X_train

#X_test_cv = scaler.transform(X_test_cv)  
# apply same transformation to test data
#X_test_cv = scaler.transform(X_test) 
X_test = scaler.transform(X_test) 

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

item = list(ParameterGrid(param_grid))[0]


def makeSimpleNN(X_train,Y_train,item=None):
	clf = MLPClassifier() ## bring to developers? or check with RF? 
	## see if you can use your data in R or matlab 
	
	#if item != None:
		#clf.set_params(**item)
	clf = MultiOutputClassifier(clf,n_jobs=-1) ## if you turn this off it splits the data 50/50 
	clf = clf.fit(X_train, Y_train) 
	print(clf.get_params())
	
	return(clf)
	
def manualSingleValidationParameters(param_grid,X_train,X_validate,Y_train,Y_validate):
	print("Fitting the following parameters:")
	print(param_grid,end="\n\n")

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

	print("Fitting the following parameters:")
	print(param_grid,end="\n\n")

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
		print("Param",str(index+1),end=" ")
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
			
	print("\n\nAll average train scores:")
	print(classifier_trains)

	print("\nBest Train Scores belong to:")
	print(best_classifier_train)
	print(best_parameters_train)
	
	print("\n\nAll average test scores:")
	print(classifier_scores)

	print("\nBest Test Scores belong to:")
	print(best_classifier_score)
	print(best_parameters)
	
	print("\n\nAll average harmonic scores:")
	print(classifier_harmonic)

	print("\nBest harmonic Scores belong to:")
	print(best_classifier_harmonic)
	print(best_parameters_harmonic)

	print("\nSetting clf to be best test model")
	clf = best_model
	
	return(clf)

#clf = MLPClassifier()

clf = makeSimpleNN(X_train,Y_train,item=None)

#clf = manualSingleValidationParameters(param_grid,X_train,X_validate,Y_train,Y_validate)

#print("\n###xxx###\n")

#clf2 = MLPClassifier()
#clf2 = manualCrossValidationParameters(param_grid,X_cv_list,Y_cv_list)

## predict the values of X_train and X_test
train_predicted = clf.predict(X_train) ## "Call fit before using this method"
test_predicted = clf.predict(X_test)
validate_predicted = clf.predict(X_validate)


print("\nTraining set score: %f" % clf.score(X_train, Y_train))
print("CV set score: %f" % clf.score(X_validate, Y_validate))
print("Test set score: %f" % clf.score(X_test, Y_test))

## SOURCE OF THIS: http://sdsawtelle.github.io/blog/output/week4-andrew-ng-machine-learning-with-python.html
#plt.figure()
#plt(clf.estimator.loss_curve_)
#plt.show()

## REMOVED PLOTTING FUNCTIONS

## REDO THIS AS ONE DATASET NOT TWO
def mergeMultioutputToSingleoutput(data_predicted,Y_data):
	data_predicted_merged = []
	Y_data_merged = []
	for i in range((len(data_predicted))):

		joined = "-".join(data_predicted[i,])
		data_predicted_merged += [joined]
	
		joined2 = "-".join(Y_data[i,])
		Y_data_merged += [joined2]
	
	return([data_predicted_merged,Y_data_merged])
	
def mergeMultioutputToSingleoutput2(data):
	data_merged = []
	for i in range((len(data))):

		joined = "-".join(data[i,])
		data_merged += [joined]
	
	return(data_merged)

[test_predicted_merged,Y_test_merged] = mergeMultioutputToSingleoutput(test_predicted,Y_test)
[train_predicted_merged,Y_train_merged] = mergeMultioutputToSingleoutput(train_predicted,Y_train)
[validate_predicted_merged,Y_validate_merged] = mergeMultioutputToSingleoutput(validate_predicted,Y_validate)

test_predicted_merged = mergeMultioutputToSingleoutput2(test_predicted)
train_predicted_merged = mergeMultioutputToSingleoutput2(train_predicted)
validate_predicted_merged = mergeMultioutputToSingleoutput2(validate_predicted)
Y_test_merged = mergeMultioutputToSingleoutput2(Y_test)
Y_train_merged = mergeMultioutputToSingleoutput2(Y_train)
Y_validate_merged = mergeMultioutputToSingleoutput2(Y_validate)


## precision and recall 
print("\nTraining Precision-Recall-Fscore")
print(precision_recall_fscore_support(Y_train_merged, train_predicted_merged, average='weighted'))
print("\nCV Precision-Recall-Fscore")
print(precision_recall_fscore_support(Y_validate_merged, validate_predicted_merged, average='weighted'))
print("\nTest Precision-Recall-Fscore")
print(precision_recall_fscore_support(Y_test_merged, test_predicted_merged, average='weighted'))

if numcategories > 1:
	print("\nTRAIN ACCURACY PER AXIS LOOKED AT")
	for c in range(numcategories):
		print(Y_train[0,c],end="\t")
		Y_col = Y_train[:,c]
		train_col = train_predicted[:,c]
		print(accuracy_score(Y_col,train_col))
	print("\nCV ACCURACY PER AXIS LOOKED AT")
	for c in range(numcategories):
		print(Y_validate[0,c],end="\t")
		Y_col = Y_validate[:,c]
		cv_col = validate_predicted[:,c]
		print(accuracy_score(Y_col,cv_col))
	print("\nTEST ACCURACY PER AXIS LOOKED AT")
	for c in range(numcategories):
		print(Y_test[0,c],end="\t")
		Y_col = Y_test[:,c]
		test_col = test_predicted[:,c]
		print(round(accuracy_score(Y_col,test_col),2))


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

[true_positives,false_positives,true_negatives,false_negatives] = trueFalsePosNegValues(classOrderLs,counts)

print("\n",classOrderLs)
print("TP:",true_positives)
print("FP:",false_positives)
print("TN:",true_negatives)
print("FN:",false_negatives)

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

outputCountsToFile(classOrderLs,df,outputString)
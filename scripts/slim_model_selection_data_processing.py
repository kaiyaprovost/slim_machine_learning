#!/usr/bin/python -u
# coding: utf-8

import numpy as np ##
import pandas as pd ## 
import sys ## 
from sklearn.model_selection import train_test_split ##
from sklearn.preprocessing import StandardScaler ## 

def readInFilenames():
	print("\nReading data")
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
	return(filenames)

def inputDataFiles(filenames):

	## suck up the data

	tempname = filenames[0].split("_")[-1].split(".")[0]
	print(tempname)
	type = filenames[0].split("_")[-1].split(".")[1] ## popgenome or sumstats or combo
	print(type)
	outputString = filenames[0].split("_")[0]+"_"+type
	print(outputString)

	numcategories = len(tempname.split("-"))

	if type == "POPGENOME":
		X = np.array([]).reshape(0,46) ## changed from 36 to 46
	elif type == "SUMSTATS":
		X = np.array([]).reshape(0,5)
	elif type == "COMBO":
		X = np.array([]).reshape(0,11) ## 33 or 11
	else:
		X = np.array([]).reshape(0,0)

	y = []
	classOrderLs = []
	outputString = ""
	for i in range(len(filenames)):
		## read the file 
		file = filenames[i]
		print(file)
	
		## count the lines of the file
	
		if type == "SUMSTATS":
			Xi = np.loadtxt(file,usecols=(1,3,5,7,9))
		else:
			Xi = np.loadtxt(file)
		
		## print the shapes
		print(X.shape)
		print(Xi.shape)

		
		## concatenate the file to the main dataframe
		try:
			X = np.concatenate((X,Xi))
		except:
			X.reshape(X.shape[0],Xi.shape[1]) ## causing an error 
			X = np.concatenate((X,Xi))
	
		## add the name to classOrderLs
		name = file.split("_")[-1].split(".")[0]
		splitname = name.split("-")
	
		classOrderLs = classOrderLs + [name]
		outputString = file.split("_")[0]
	
		y = y + ([splitname]*len(Xi))
	Y = np.array(y)
	
	return([X,Y,numcategories,classOrderLs,outputString])

def manualTrainTestSplit(X,Y,N=10):
	## split data into training and test 
	print("\nCreating test, validate, and training data")

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

def outputXYdatasets(dataset,writename):
	pd.DataFrame(dataset).to_csv(writename,header=False,index=False)

def outputXYdatasetsMult(X_train,X_validate,X_test,Y_train,Y_validate,Y_test):

	## TRY RE-OUTPUT AND RE-INPUT THE DATA
	pd.DataFrame(X_train).to_csv("X_train.temp",header=False,index=False)
	pd.DataFrame(X_validate).to_csv("X_validate.temp",header=False,index=False)
	pd.DataFrame(X_test).to_csv("X_test.temp",header=False,index=False)

	pd.DataFrame(Y_train).to_csv("Y_train.temp",header=False,index=False)
	pd.DataFrame(Y_validate).to_csv("Y_validate.temp",header=False,index=False)
	pd.DataFrame(Y_test).to_csv("Y_test.temp",header=False,index=False)

def inputXYdatasets(numcategories,readname,xy="X"):
	if xy == "X":
		dataset = np.loadtxt(readname,delimiter=",")
	elif xy == "Y":
		dataset = np.loadtxt(readname,delimiter=",",dtype=str).reshape(dataset.shape[0],numcategories)
	else:
		raise Exception("xy must be either X or Y")
	return(dataset)

def inputXYdatasetsMult(numcategories):

	X_train = np.loadtxt("X_train.temp",delimiter=",")
	X_validate = np.loadtxt("X_test.temp",delimiter=",")
	X_test = np.loadtxt("X_validate.temp",delimiter=",")

	Y_train = np.loadtxt("Y_train.temp",delimiter=",",dtype=str).reshape(Y_train.shape[0],numcategories)
	Y_validate = np.loadtxt("Y_test.temp",delimiter=",",dtype=str).reshape(Y_validate.shape[0],numcategories)
	Y_test = np.loadtxt("Y_validate.temp",delimiter=",",dtype=str).reshape(Y_test.shape[0],numcategories)
	
	return([X_train,X_test,X_validate,Y_train,Y_test,Y_validate])

def featureScaling(X_train,X_test,X_validate):
	# Don't cheat - fit only on training data
	scaler = StandardScaler(copy=False, with_mean=True, with_std=True).fit(X_train)
	X_train_scaled = scaler.transform(X_train)  
	X_test_scaled = scaler.transform(X_test) 
	X_validate_scaled = scaler.transform(X_validate)
	
	return([X_train_scaled,X_test_scaled,X_validate_scaled])

def main():

	print("\n------------------------------------\nRUNNING DATA PROCESSING MODULE ALONE\n------------------------------------\n")

	filenames = readInFilenames()
	
	[X,Y,numcategories,classOrderLs,outputString] = inputDataFiles(filenames)
	
	[X_train,X_test,X_validate,X_cv_list,Y_train,Y_test,Y_validate,Y_cv_list] = manualTrainTestSplit(X,Y,N=10)
	
	outputXYdatasets(dataset=X_train,writename="X_train.temp")
	outputXYdatasets(dataset=X_validate,writename="X_validate.temp")
	outputXYdatasets(dataset=X_test,writename="X_test.temp")
	outputXYdatasets(dataset=Y_train,writename="Y_train.temp")
	outputXYdatasets(dataset=Y_validate,writename="Y_validate.temp")
	outputXYdatasets(dataset=Y_test,writename="Y_test.temp")
	
	print("\nShapes of features:")
	print("Full dataset X, Y")
	print(X.shape,Y.shape)
	print("\nTrain X, Y")
	print(X_train.shape, Y_train.shape)
	print("CV X, Y")
	print(X_validate.shape, Y_validate.shape)
	print("Test X, Y")
	print(X_test.shape, Y_test.shape)
	
	[X_train_scaled,X_test_scaled,X_validate_scaled] = featureScaling(X_train,X_test,X_validate)
	
	outputXYdatasets(dataset=X_train_scaled,writename="X_train_scaled.temp")
	outputXYdatasets(dataset=X_validate_scaled,writename="X_validate_scaled.temp")
	outputXYdatasets(dataset=X_test_scaled,writename="X_test_scaled.temp")
	
	print("\n---------------------------------------\nFINISHING DATA PROCESSING MODULE, ALONE\n---------------------------------------\n")

	
if __name__ == "__main__":
	main()



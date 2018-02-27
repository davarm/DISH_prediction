import pandas as pd 
from sklearn import svm
import numpy as np
from sklearn.svm import SVC
from sklearn.externals import joblib
import pickle
from pandas import Series


# Read Database
df = pd.read_csv('DISH_database.csv',  sep = ',',skipinitialspace = False)

# Function to split hot arrays into individual columns
def split_columns(dataframe,hot_array):
	 length=dataframe[hot_array][0]
	 length=len(length.split(','))
	 s                     = dataframe[hot_array].str.split('').apply(Series, 1).stack	()
	 s.index               = s.index.droplevel(-1)
	 s.name                = hot_array
	 del dataframe[hot_array]
	 dataframe =dataframe.join(s.apply(lambda x: Series(x.split(','))))
	 for i in range(length):
		 x=str(hot_array)+str(i)
		 dataframe=dataframe.rename(columns = {i:x})
	 return dataframe

# Split the secondary structure array
df = split_columns(df,'cys1_ss_array')


################################################
#----------------------------------
# For X1 Prediction, Stage 1
#----------------------------------
################################################

# Required inputs for x1_stage_1 [60, 'other']
configuration = [60,'other']

x1_stage1_df  = df[['x1',
					'cys1_phi',
					'cys1_psi',
					'cys1_ss_array0',
					'cys1_ss_array1',
					'cys1_ss_array2',
					'cys1_Ha',
					'cys1_Ca',
					'cys1_Cb',
					'cys2_Ca',
					'cys2_Cb']]


# Convert X1 into a binary
def x1_stage1_index(x1):
	if x1 == 60.0:
		x1_index = configuration.index(x1)
	if x1 != 60.0:
		x1_index = configuration.index('other')
	return (x1_index)

# Y is the target
y = x1_stage1_df['x1'].apply(x1_stage1_index)

# X represents the inputs
x = x1_stage1_df[[	'cys1_phi',
					'cys1_psi',
					'cys1_ss_array0',
					'cys1_ss_array1',
					'cys1_ss_array2',
					'cys1_Ha',
					'cys1_Ca',
					'cys1_Cb',
					'cys2_Ca',
					'cys2_Cb']]


clf = svm.SVC(kernel = 'rbf',probability=True,gamma=0.004,C=5)
clf.fit(x,y)
joblib.dump(clf, 'x1_stage1.pkl')

################################################
#----------------------------------
# For X1 Prediction, Stage 2
#----------------------------------
################################################
y = []
x = []
configuration = []

# Predict as [-60, 180]
configuration=[-60.0,180.0]

x1_stage2_df  = df[['x1',
					'cys1_phi',
					'cys1_psi',
					'cys2_psi',
					'cys1_ss_array0',
					'cys1_ss_array1',
					'cys1_ss_array2',
					'cys1_Ha',
					'cys1_Ca',
					'cys1_Cb',
					'cys1_N',
					'cys1_Hn',
					'cys1_b_vdw',
					'cys2_Ca',
					'cys2_Cb',
					'cys2_N']]
					
x1_stage2_df  = x1_stage2_df.loc[x1_stage2_df['x1'] != 60.0]
x1_stage2_df  = x1_stage2_df.reset_index(drop = True)

def x1_stage2_index(x1):
	x1_index = configuration.index(x1)
	return (x1_index)

y = x1_stage2_df['x1'].apply(x1_stage2_index)
x = x1_stage2_df[[
				'cys1_phi',
				'cys1_psi',
				'cys2_psi',
				'cys1_ss_array0',
				'cys1_ss_array1',
				'cys1_ss_array2',
				'cys1_Ha',
				'cys1_Ca',
				'cys1_Cb',
				'cys1_N',
				'cys1_Hn',
				'cys1_b_vdw',
				'cys2_Ca',
				'cys2_Cb',
				'cys2_N']]

clf2=svm.SVC(kernel='rbf',probability=True,gamma=9e-6,C=3000)
clf2.fit(x,y)
joblib.dump(clf2, 'x1_stage2.pkl')

################################################
#----------------------------------
# For X2 Prediction, Stage 1
#----------------------------------
################################################

configuration = []
x             = []
y             = []

configuration = ['other',-60.0]

x2_stage1_df  = df[['x2',
					'x1',
					'cys1_ss_array0',
					'cys1_ss_array1',
					'cys1_ss_array2',	
					'cys1_N',
					'cys1_Ca',
					'cys1_Cb',
					'cys1_Hn',
					'cys2_Ca',
					'cys2_Cb']]

def x2_stage1_index(x2):
	if x2 == -60.0:
		x2_index = configuration.index(x2)
	if x2 != -60.0:
		x2_index = configuration.index('other')
	return (x2_index)

y = x2_stage1_df['x2'].apply(x2_stage1_index)
x = x2_stage1_df[['x1',
					'cys1_ss_array0',
					'cys1_ss_array1',
					'cys1_ss_array2',	
					'cys1_N',
					'cys1_Ca',
					'cys1_Cb',
					'cys1_Hn',
					'cys2_Ca',
					'cys2_Cb']]

clf3=svm.SVC(kernel='rbf',probability=True,gamma=0.006,C=40)
clf3.fit(x,y)
joblib.dump(clf3, 'x2_stage1.pkl')


################################################
#----------------------------------
# For X2 Prediction, Stage 2
#----------------------------------
################################################

configuration = []
x             = []
y             = []
x2_stage2_df = df[[	'x2',
					'x1',
					'cys1_ss_array0',
					'cys1_ss_array1',
					'cys1_ss_array2',					
					'cys1_N',
					'cys1_Ca',
					'cys1_Cb',
					'cys1_Ha',
					'cys2_Ca',
					'cys2_Cb']]

x2_stage2_df  = x2_stage2_df.loc[x2_stage2_df['x2'] != -60.0]
x2_stage2_df  = x2_stage2_df.reset_index(drop = True)


configuration=[60.0,180.0]
def x2_stage2_index(x2):
	x2_index = configuration.index(x2)
	return (x2_index)

y = x2_stage2_df['x2'].apply(x2_stage2_index)
x = x2_stage2_df[[	'x1',
					'cys1_ss_array0',
					'cys1_ss_array1',
					'cys1_ss_array2',					
					'cys1_N',
					'cys1_Ca',
					'cys1_Cb',
					'cys1_Ha',
					'cys2_Ca',
					'cys2_Cb']]

clf4=svm.SVC(kernel='rbf',probability=True,gamma=0.002,C=90)
clf4.fit(x,y)
joblib.dump(clf4, 'x2_stage2.pkl')#
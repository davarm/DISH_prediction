import pandas as pd 
from sklearn import svm
import numpy as np
from sklearn.svm import SVC
from sklearn.externals import joblib
import pickle

df = pd.read_csv('DISH_database.csv',  sep = ',',skipinitialspace = False)

################################################
#----------------------------------
# For X1 Prediction, Stage 1
#----------------------------------
################################################
print df[['cys1_N','cys1_Ca','cys1_Cb','cys2_Cb','cys2_Ca','cys1_before_psi','cys1_phi','cys1_psi','x1','PDB','cys1_residue','cys2_residue']]
x1_stage1_df  = df[['cys1_N','cys1_Ca','cys1_Cb','cys2_Cb','cys2_Ca','cys1_before_psi','cys1_phi','cys1_psi','x1']]

configuration = ['other',60.0]
def x1_stage1_index(x1):
	if x1 == 60.0:
		x1_index = configuration.index(x1)
	if x1 != 60.0:
		x1_index = configuration.index('other')
	return (x1_index)

y = x1_stage1_df['x1'].apply(x1_stage1_index)
x = x1_stage1_df[['cys1_N','cys1_Ca','cys1_Cb','cys2_Cb','cys2_Ca','cys1_before_psi','cys1_phi','cys1_psi']]
clf = svm.SVC(kernel = 'rbf',probability=True,gamma=0.001,C=50)
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
x1_stage2_df  = df[['cys1_before_vdw','cys1_N','cys1_Ca','cys1_Cb','cys1_HN','cys2_Cb','cys2_Ca','cys1_before_psi','cys1_phi','cys1_psi','cys2_psi','x1']]
x1_stage2_df  = x1_stage2_df.loc[x1_stage2_df['x1'] != 60.0]
x1_stage2_df  = x1_stage2_df.reset_index(drop = True)

configuration=[-60.0,180.0]
def x1_stage2_index(x1):
	x1_index = configuration.index(x1)
	return (x1_index)

y = x1_stage2_df['x1'].apply(x1_stage2_index)
x = x1_stage2_df[['cys1_before_vdw','cys1_N','cys1_Ca','cys1_Cb','cys1_HN','cys2_Cb','cys2_Ca','cys1_before_psi','cys1_phi','cys1_psi','cys2_psi']]

clf2=svm.SVC(kernel='rbf',probability=True,gamma=7e-6,C=1000)
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

x2_stage1_df  = df[['cys1_N','cys1_Ca','cys1_Cb','cys2_Cb','cys2_Ca','x1','x2']]

def x2_stage1_index(x2):
	if x2 == -60.0:
		x2_index = configuration.index(x2)
	if x2 != -60.0:
		x2_index = configuration.index('other')
	return (x2_index)

y = x2_stage1_df['x2'].apply(x2_stage1_index)
x = x2_stage1_df[['cys1_N','cys1_Ca','cys1_Cb','cys2_Cb','cys2_Ca','x1']]

clf3=svm.SVC(kernel='rbf',probability=True,gamma=0.008,C=30)
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
x2_stage2_df = df[['cys1_N','cys1_Ca','cys1_Cb','cys2_Cb','cys2_Ca','x1','x2']]
x2_stage2_df  = x2_stage2_df.loc[x2_stage2_df['x2'] != -60.0]
x2_stage2_df  = x2_stage2_df.reset_index(drop = True)


configuration=[60.0,180.0]
def x2_stage2_index(x2):
	x2_index = configuration.index(x2)
	return (x2_index)

y = x2_stage2_df['x2'].apply(x2_stage2_index)
x = x2_stage2_df[['cys1_N','cys1_Ca','cys1_Cb','cys2_Cb','cys2_Ca','x1']]

clf4=svm.SVC(kernel='rbf',probability=True,gamma=2e-5,C=1e6)
clf4.fit(x,y)
joblib.dump(clf4, 'x2_stage2.pkl')#
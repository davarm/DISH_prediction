#------------------------------------------------------------------------
# Prediction of X1 and X2 angles using DISH (Armstrong, et al. 2017)
# Generate the correct inpurt format from 'generate_inputs.py'
#------------------------------------------------------------------------

import pandas as pd
import numpy  as np 
import sys
from sklearn.externals import joblib
import os
os.getcwd()

path    = ('./SVMs/')

# Load the SVMs, saved as .pkl files in 'SVMs' directory
clf_x1_stage1     = joblib.load(path+'x1_stage1.pkl')
clf_x1_stage2     = joblib.load(path+'x1_stage2.pkl')
clf_x2_stage1     = joblib.load(path+'x2_stage1.pkl')
clf_x2_stage2     = joblib.load(path+'x2_stage2.pkl')

# Read the desired inputs into a dataframe
df                = pd.read_csv('DISH.csv',sep = ',', skipinitialspace = False)



def x1_prediction(inputs):

	#----------------------------------
	# Stage 1 is the prediciton of X1 as +60 or 'Other'
	# Desired inputs in the order of:
	# 		[N, CA, CB, Hemi-CB, Hemi-CA, before_psi, phi, psi]
	# Prediction: 0 = Other, 1 = +60
	#----------------------------------

	x1_stage1_array  = ['Other', 60]

	x1_stage1_inputs = inputs[['N'		   ,
					   		   'CA'		   ,
					   		   'CB'		   ,
					   		   'hemi_CB'   ,
					   		   'hemi_CA'   ,
					   		   'before_psi',
					   		   'phi'       ,
					   		   'psi']]

	x1_stage1_inputs            = np.array(x1_stage1_inputs).reshape((1, -1))
	x1_stage1_prediction_scores = (clf_x1_stage1.predict_proba(x1_stage1_inputs)[0]).tolist()
	x1_stage1_prediction_index  = x1_stage1_prediction_scores.index(max(x1_stage1_prediction_scores))
	x1_stage1_prediction        = x1_stage1_array[x1_stage1_prediction_index]

	if x1_stage1_prediction == 60:
		x1_prediction = 60
		x1_prob       = max(x1_stage1_prediction_scores)
	#----------------------------------
	# If the prediction is 'Other', move to 'STAGE 2' SVM
	# Predciton of -60 or 180
	#----------------------------------

	if x1_stage1_prediction == 'Other':
		x1_stage2_array  = [-60, 180]

		x1_stage2_inputs = inputs[['before_vdw_radi',
								   'N'		        ,
								   'CA'		        ,
								   'CB'		        ,
								   'HN'             ,
								   'hemi_CB'        ,
								   'hemi_CA'        ,
								   'before_psi'     ,
								   'phi'	        ,
								   'psi'            ,
								   'hemi_psi']]

		x1_stage2_inputs            = np.array(x1_stage2_inputs).reshape((1, -1))
		x1_stage2_prediction_scores = (clf_x1_stage2.predict_proba(x1_stage2_inputs)[0]).tolist()
		x1_stage2_prediction_index  = x1_stage2_prediction_scores.index(max(x1_stage2_prediction_scores))
		x1_prediction               = x1_stage2_array[x1_stage2_prediction_index]
		x1_prob                     = max(x1_stage2_prediction_scores)

	return(x1_prediction,x1_prob)



#----------------------------------
# X2 PREDICTION: Uses X1 as an Input
#----------------------------------
def x2_prediction(inputs):

	# Stage 1 is the prediciton of x2 as -+60 or 'Other'
	# Desired inputs in the order of:
	# 		[N, CA, CB, Hemi-CB, Hemi-CA, before_psi, phi, psi]
	# Prediction: 0 = Other, 1 = -60
	x2_stage1_array  = ['Other', -60]

	x2_stage1_inputs = inputs[['N'      ,
							   'CA'     ,
							   'CB'     ,
							   'hemi_CB',
							   'hemi_CA',
							   'x1']]


	x2_stage1_inputs            = np.array(x2_stage1_inputs).reshape((1, -1))
	x2_stage1_prediction_scores = (clf_x2_stage1.predict_proba(x2_stage1_inputs)[0]).tolist()
	x2_stage1_prediction_index  = x2_stage1_prediction_scores.index(max(x2_stage1_prediction_scores))
	x2_stage1_prediction        = x2_stage1_array[x2_stage1_prediction_index]

	if x2_stage1_prediction == -60:
		x2_prediction = -60
		x2_prob       = max(x2_stage1_prediction_scores)

	#----------------------------------
	# If the prediction is 'Other', move to 'STAGE 2' SVM
	# Predciton of -60 or 180
	#----------------------------------

	if x2_stage1_prediction == 'Other':
		x2_stage2_array  = [60, 180]

		x2_stage2_inputs = inputs[['N'      ,
								   'CA'     ,
								   'CB'     ,
								   'hemi_CB',
								   'hemi_CA',
								   'x1']]

		x2_stage2_inputs            = np.array(x2_stage2_inputs).reshape((1, -1))
		x2_stage2_prediction_scores = (clf_x2_stage2.predict_proba(x2_stage2_inputs)[0]).tolist()
		x2_stage2_prediction_index  = x2_stage2_prediction_scores.index(max(x2_stage2_prediction_scores))
		x2_prediction               = x2_stage2_array[x2_stage2_prediction_index]
		x2_prob                     = max(x2_stage2_prediction_scores)

	return(x2_prediction, x2_prob)



print '{:5s} {:4s} {:8s} {:4s} {:2s}'.format('Res','X1',"Prob.",'X2',"Prob.")
for index,row in df.iterrows():
	x1_pred   = x1_prediction(row)
	row['x1'] = x1_pred[0]
	x2_pred   = x2_prediction(row)
	print '{:3d} {:4d} {:8.4f} {:4d} {:8.4f}'.format(row['residue_number'],x1_pred[0],x1_pred[1],x2_pred[0],x2_pred[1])








































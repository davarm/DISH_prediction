#------------------------------------------------------------------------
# Prediction of X1 and X2 angles using DISH (Armstrong, et al. 2017)
# Generate the correct inpurt format from 'generate_inputs.py'
#------------------------------------------------------------------------

import pandas as pd
from pandas import Series
import numpy  as np 
import sys
from sklearn.externals import joblib
import os
os.getcwd()

path    = ('./compile_SVMs/')
get = open("./peptides/"+sys.argv[1]+"/DISH_prediction.txt",'w')

# Load the SVMs, saved as .pkl files in 'SVMs' directory
clf_x1_stage1     = joblib.load(path+'x1_stage1.pkl')
clf_x1_stage2     = joblib.load(path+'x1_stage2.pkl')
clf_x2_stage1     = joblib.load(path+'x2_stage1.pkl')
clf_x2_stage2     = joblib.load(path+'x2_stage2.pkl')

# Read the desired inputs into a dataframe
df                = pd.read_csv("./peptides/"+sys.argv[1]+'/DISH_inputs.csv',sep = ',', skipinitialspace = False)
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


df = split_columns(df,'ss_array')
# print df 


def x1_prediction(inputs):

	#----------------------------------
	# Stage 1 is the prediciton of X1 as +60 or 'Other'
	# Desired inputs in the order of:
	# 		[N, CA, CB, Hemi-CB, Hemi-CA, before_psi, phi, psi]
	# Prediction: 0 = Other, 1 = +60
	#----------------------------------

	x1_stage1_array  = [60, 'Other']

	x1_stage1_inputs = inputs [[	
									'phi',
									'psi',
									'ss_array0',
									'ss_array1',
									'ss_array2',				
									'HA',
									'CA',
									'CB',
									'hemi_CA',
									'hemi_CB']]

	x1_stage1_inputs                 = np.array(x1_stage1_inputs).reshape((1, -1))
	x1_stage1_prediction_index       = (clf_x1_stage1.predict(x1_stage1_inputs)[0])
	x1_stage1_prediction_probability = ((clf_x1_stage1.predict_proba(x1_stage1_inputs)[0]).tolist())
	x1_stage1_prediction             = x1_stage1_array[x1_stage1_prediction_index]

	if x1_stage1_prediction == 60:
		x1_prediction = 60
		x1_prob       = x1_stage1_prediction_probability[x1_stage1_prediction_index]
	#----------------------------------
	# If the prediction is 'Other', move to 'STAGE 2' SVM
	# Predciton of -60 or 180
	#----------------------------------

	if x1_stage1_prediction == 'Other':
		x1_stage2_array  = [-60, 180]


		x1_stage2_inputs = inputs[[	'phi',
									'psi',
									'hemi_psi',
									'ss_array0',
									'ss_array1',
									'ss_array2',					
									'HA',
									'CA',
									'CB',
									'N',
									'HN',
									'before_vdw_radi',
									'hemi_CA',
									'hemi_CB',
									'hemi_N']]

		# print x1_stage2_inputs
		x1_stage2_inputs                 = np.array(x1_stage2_inputs).reshape((1, -1))
		
		x1_stage2_inputs                 = np.array(x1_stage2_inputs).reshape((1, -1))
		x1_stage2_prediction_index       = (clf_x1_stage2.predict(x1_stage2_inputs)[0])
		x1_stage2_prediction_probability = ((clf_x1_stage2.predict_proba(x1_stage2_inputs)[0]).tolist())
		x1_stage2_prediction             = x1_stage2_array[x1_stage2_prediction_index]
		
		x1_prob                          = x1_stage2_prediction_probability[x1_stage2_prediction_index]
		x1_prediction                    = x1_stage2_prediction  
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

	x2_stage1_inputs =  inputs[['x1',
								'ss_array0',
								'ss_array1',
								'ss_array2',					
								'N',
								'CA',
								'CB',
								'HN',
								'hemi_CA',
								'hemi_CB']]


	x2_stage1_inputs                 = np.array(x2_stage1_inputs).reshape((1, -1))
	x2_stage1_prediction_index       = (clf_x2_stage1.predict(x2_stage1_inputs)[0])
	x2_stage1_prediction_probability = ((clf_x2_stage1.predict_proba(x2_stage1_inputs)[0]).tolist())
	x2_stage1_prediction             = x2_stage1_array[x2_stage1_prediction_index]

	if x2_stage1_prediction == -60:
		x2_prediction = -60
		x2_prob       = x2_stage1_prediction_probability[x2_stage1_prediction_index]

	#----------------------------------
	# If the prediction is 'Other', move to 'STAGE 2' SVM
	# Predciton of -60 or 180
	#----------------------------------

	if x2_stage1_prediction == 'Other':
		x2_stage2_array  = [60, 180]

		x2_stage2_inputs = inputs [['x1',
					'ss_array0',
					'ss_array1',
					'ss_array2',					
					'N',
					'CA',
					'CB',
					'HA',
					'hemi_CA',
					'hemi_CB']]

		x2_stage2_inputs                 = np.array(x2_stage2_inputs).reshape((1, -1))
		x2_stage2_inputs                 = np.array(x2_stage2_inputs).reshape((1, -1))
		x2_stage2_prediction_index       = (clf_x2_stage2.predict(x2_stage2_inputs)[0])
		x2_stage2_prediction_probability = ((clf_x2_stage2.predict_proba(x2_stage2_inputs)[0]).tolist())
		x2_stage2_prediction             = x2_stage2_array[x2_stage2_prediction_index]
		x2_prediction                    = x2_stage2_prediction
		x2_prob                          = x2_stage2_prediction_probability[x2_stage2_prediction_index]

	return(x2_prediction, x2_prob)



print '{:5s} {:4s} {:8s} {:4s} {:2s}'.format('Res','X1',"Prob.",'X2',"Prob.")
get.write('{:5s} {:4s} {:8s} {:4s} {:2s}'.format('Res','X1',"Prob.",'X2',"Prob."))
get.write('\n')
for index,row in df.iterrows():
	# print row 
	x1_pred   = x1_prediction(row)
	row['x1'] = x1_pred[0]
	x2_pred   = x2_prediction(row)
	results = ('{:3.0f} {:4d} {:8.4f} {:4d} {:8.4f}'.format(row['residue_number'],x1_pred[0],x1_pred[1],x2_pred[0],x2_pred[1]))
	print results
	get.write(str(results))
	get.write('\n')







































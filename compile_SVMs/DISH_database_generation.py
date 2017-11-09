import pandas as pd 
import numpy as np
import itertools
import os	
module_dir = os.path.dirname(os.path.abspath(__file__))																				
#-----------
#To be used for creating hot arrays-------------------------------------------------------------
#Read the disulfide database into a pandas dataframe
#------------------------------------------------------------------------




df = pd.read_csv(os.path.join(module_dir,'angles_dssp.csv'),  sep = ',',skipinitialspace = False)
# Changing all unassigned chemical shifts to 0
df = df.replace([9999.000],['0'])

#------------------------------------
#Generate a list of PDB numbers
#------------------------------------
pdb_list       =  df['PDB'].tolist()
pdb_list       =  list(set(pdb_list))


#----------------------------------
# VDW RADI dict
#----------------------------------
vdw_radi_dict={}
vdw_radi_dict['A']=67
vdw_radi_dict['R']=148
vdw_radi_dict['N']=96
vdw_radi_dict['D']=91
vdw_radi_dict['C']=86
vdw_radi_dict['c']=86
vdw_radi_dict['Q']=114
vdw_radi_dict['E']=109
vdw_radi_dict['G']=48
vdw_radi_dict['H']=118
vdw_radi_dict['I']=124
vdw_radi_dict['L']=124
vdw_radi_dict['K']=135
vdw_radi_dict['M']=124
vdw_radi_dict['F']=135
vdw_radi_dict['P']=90
vdw_radi_dict['S']=73
vdw_radi_dict['T']=93
vdw_radi_dict['W']=163
vdw_radi_dict['Y']=141
vdw_radi_dict['V']=105

def vdw_radi(residue):
	vdw_radi_value = vdw_radi_dict[residue]
	return (vdw_radi_value)

#------------------------------------------------------------------------
# Change the residue before to a vdw radi value
#------------------------------------------------------------------------	

df['Cys1 before residue']  = df['Cys1 before residue'].apply(vdw_radi)
df['Cys2 before residue']  = df['Cys2 before residue'].apply(vdw_radi)
#------------------------------------------------------------------------
#Applying rounding functions to the cysteine dihedral angles
#------------------------------------------------------------------------
def x1_rounded(x1):
			x1=float(x1)
			if (x1 <=  90) & (x1 >= 30):
				x1= 60
			if (x1 >= -90)  & (x1 <= -30):
				x1 = -60   
			if (x1 <=  180) & (x1 >= 150):
				x1= 180
			if (x1 >= -180) & (x1 <= -150):
				x1=180
			return(x1)
	
def x2_rounded(x2):
			x2=float(x2)
			if (x2 <=  120) & (x2 >= 30):
				x2= 60
			if (x2 >= -120) & (x2 <= -30):
				x2 = -60   
			if (x2 <=  180) & (x2 >= 150):
				x2= 180
			if (x2 >= -180) & (x2 <= -150):
				x2=180
			return(x2)
	
def x3_rounded(x3):
			x3=float(x3)
			if (x3  <=  120)  & (x3 >= 60):
				x3   = 90
			if (x3  >=  -120) & (x3 <= -60):
				x3   = -90
			return(x3)



#Apply funcitons
df['x1' ]  = df['x1'].apply(x1_rounded)
df['x1b'] = df['x1b'].apply(x1_rounded)
df['x2' ]  = df['x2'].apply(x2_rounded)
df['x2b'] = df['x2b'].apply(x2_rounded)
df['x3' ]  = df['x3'].apply(x3_rounded)


	
#------------------------------------------------------------------------
# REMOVE ALL Cysteine residues where dihedral angles are outliers
#-----------------------------------------------------1-------------------
configuration=[60,-60,180,90,-90]
df=df[df['x3' ].isin(configuration)]
df=df[df['x1' ].isin(configuration)]
df=df[df['x2' ].isin(configuration)]
df=df[df['x1b'].isin(configuration)]
df=df[df['x2b'].isin(configuration)]


#------------------------------------------------------------------------
# ADDING IN HOT ARRAYS REQUIRED AS INPUTS AND FOR TESTING
#------------------------------------------------------------------------

#X2 Angles
def x2_array(x2):
	configuration                         = [-60,180,60]
	x2_hot_array                          = [0. for _ in range(3)]
	x2_hot_array[configuration.index(x2)] = 1
	x2_hot_array                          = [str(x) for x in x2_hot_array]
	x2_hot_array                          = ",".join(x2_hot_array)
	return (x2_hot_array)

#df['x2_array']  = df['x2'].apply(x2_array)
#df['x2b_array']  = df['x2b'].apply(x2_array)

#X3 Angles
def x3_array(x3):
	if x3 == -90:
		x3_hot_array = '1'
	if x3 == 90:
		x3_hot_array = '0'
	return (x3_hot_array)

df['x3_array']  = df['x3'].apply(x3_array)
	
#X1 Angles
def x1_array(x1):
	configuration                         = [-60,180,60]
	x1_hot_array                          = [0. for _ in range(3)]
	x1_hot_array[configuration.index(x1)] = 1
	x1_hot_array                          = [str(x) for x in x1_hot_array]
	x1_hot_array                          = ",".join(x1_hot_array)
	return (x1_hot_array)

df['x1_array']  = df['x1'].apply(x1_array)
df['x1b_array']  = df['x1b'].apply(x1_array)

df['x2_array']  = df['x2'].apply(x2_array)
df['x2b_array']  = df['x2b'].apply(x2_array)


	
	
#------------------------------------------------------------------------
# Create df_split, where divide the cysteines into invidiual residues
#Use dictionary to make sure columns are the same
#------------------------------------------------------------------------
#FOR THE FIRST CYS

new_columns=['PDB',
			'cys1_residue',
			'cys2_residue',
			'cys1_phi',
			'cys1_psi',
			'cys1_before_psi',
			'cys1_before_vdw',
			'cys2_psi'
			'cys1_N',
			'cys1_HN',
			'cys1_Ca',
			'cys1_Cb',
			'cys2_Ca',
			'cys2_Cb',
			'x1',
			'x1_array'
			'x2',
			'x2_array']
	
	
	
df_  = pd.DataFrame(index = None,columns = new_columns)
df_  =  pd.DataFrame({ 
				'PDB' 				: df['PDB'],
				'cys1_residue'      : df['residue 1'],
				'cys2_residue'      : df['residue 2'],
				'cys1_phi'  		: df['phi'],
				'cys1_psi'  		: df['psi'],
				'cys1_before_psi'   : df['Cys1 before psi'],
				'cys1_before_vdw'   : df['Cys1 before residue'],
				'cys2_psi'          : df['psi_x'],
				'cys1_N'            : df['Cys1 N'],
				'cys1_HN'            : df['Cys1 HN'],
				'cys1_Ca'           : df['Cys1 CA'],
				'cys1_Cb'           : df['Cys1 CB'],
				'cys2_Ca'           : df['Cys2 CA'],
				'cys2_Cb'           : df['Cys2 CB'],
				'x1'  				: df['x1'],
				'x1_array'  		: df['x1_array'],
				'x2'  				: df['x2'],
				'x2_array'  		: df['x2_array'],
				})	

 		                                  
	
# df_ = df_.astype(str)

df__  = pd.DataFrame(index = None,columns = new_columns)
df__  =  pd.DataFrame({                                                
				'PDB' 				: df['PDB'],
				'cys1_residue'      : df['residue 2'],
				'cys2_residue'      : df['residue 1'],
				'cys1_phi'  		: df['phi_x'],
				'cys1_psi'  		: df['psi_x'],
				'cys1_before_psi'   : df['Cys2 before psi'],
				'cys1_before_vdw'   : df['Cys2 before residue'],
				'cys2_psi'          : df['psi'],
				'cys1_N'            : df['Cys2 N'],
				'cys1_HN'            : df['Cys2 HN'],
				'cys1_Ca'           : df['Cys2 CA'],
				'cys1_Cb'           : df['Cys2 CB'],
				'cys2_Ca'           : df['Cys1 CA'],
				'cys2_Cb'           : df['Cys1 CB'],
				'x1'  				: df['x1b'],
				'x1_array'  		: df['x1b_array'],
				'x2'  				: df['x2b'],
				'x2_array'  		: df['x2b_array']
				})	


# df__ = df__.astype(str)

df_split = df_.append(df__, ignore_index = True)
df_split = df_split[['PDB',
					'cys1_residue',
					'cys2_residue',
					'cys1_phi',
					'cys1_psi',
					'cys1_before_psi',
					'cys1_before_vdw',
					'cys2_psi',
					'cys1_N',
					'cys1_HN',
					'cys1_Ca',
					'cys1_Cb',
					'cys2_Ca',
					'cys2_Cb',
					'x1',
					'x1_array',
					'x2',
					'x2_array']]
					
					
df_split = df_split.astype(str)
df_split.to_csv('DISH_database.csv',index=False)












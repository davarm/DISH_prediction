import os
from glob import glob
import itertools
import pandas as pd 
import numpy as np
from collections import OrderedDict
import os,sys

# Script to generate the required inputs for the DISH program
# Chemical shift, sequence and backbone angles are extracted from TALOS-N files (README)
# Must have at minimum the 'pred.tab', 'predAdjCS.tab' and 'predSS.tab' files
# These should NOT be changed. The sequence is extracted from the pred.tab file
# Must also have an additional 'connectivity.txt' file that defines correct cysteine pairings 
# The program does not predict termini cysteines and therefore they are excluded
# If chemical shifts are unassigned they are designated as '0'
# 'hemi' means an input from the corresponding hemi-cysteine of a disulfide bond (as derived from connectivity.txt)
# Backbone angles are only taken from TALOS-N if they are 'STRONG or GENEROUS' predictions
# If they are 'None' or 'Warn', inputs are substituted with the average phi and psi angles derived from the training database.
# To run specifiy folder name as first argv 
#
#				: ./generating_DISH_inputs.py 2n8e
#
#

peptide = sys.argv[1]
path = "./peptides/"+peptide
#------------------------------------------------------------------------
# Create ordered dictionaries to store all relevant information
#------------------------------------------------------------------------
nuclei_list   = ['HA','CB','CA','N','HN']
cysteine_dict =  OrderedDict()

for nuclei in nuclei_list:
	cysteine_dict['Cys1_'        + nuclei] = [] 
	cysteine_dict['Cys2_'        + nuclei] = []
	
	
cysteine_dict['PDB'                ] = []
cysteine_dict['Cys1_phi'           ] = []
cysteine_dict['Cys1_psi'           ] = []
cysteine_dict['Cys2_phi'           ] = []
cysteine_dict['Cys2_psi'           ] = []
cysteine_dict['residue_1'          ] = []
cysteine_dict['residue_2'          ] = []
cysteine_dict['Cys1_after_residue' ] = []
cysteine_dict['Cys2_after_residue' ] = []
cysteine_dict['Cys1_before_residue'] = []
cysteine_dict['Cys2_before_residue'] = []
cysteine_dict['Cys1_SS'] = []
cysteine_dict['Cys2_SS'] = []


# Initiate dataframe to store information
df = pd.DataFrame()


#---------------------------------------------------
# Defining connetivity in a 'connectivity.txt' file
# Store in connectivity_list
#---------------------------------------------------
get               = open(path+'/connectivity.txt','r')
chi1_dict         = {}
cys_list          = []
connectivity_list = []
for line in get:
	line = line.strip('\n')
	if ':' in line:
			line = line.split(':')
			cys_list.append(line[0])
			cys_list.append(line[1])
			connectivity_list.append(line[0:2])
get.close()        


#--------------------------------------------------------------------------------
# Storing chemical shifts by residue number and nuclei in experimenatl_shift dict
# Chemical shifts stored in the TALOS-N format '.tab' file
# Store in the dicitonary as a combination of residue number and nuclei type
#Calculating secondary shift values is done later on to simplify process as must 
#accommodate for Proline residues changing RC values
#---------------------------------------------------------------------------------

adjusted_chemical_shift_dict = {}
sequence                     = []
sequence_dict                = {}
get                          = open(path+'/predAdjCS.tab','r') 
for line in get:
		line  = " ".join(line.split())
		if 'REMARK' in line:
			continue 
		if 'VARS' in line:
			continue
		if 'FORMAT' in line:
			continue 
		if len(line) == 0:
			continue

			
		# In TALOS-N sequence can be spread over two lines
		# Will convert and append into a singular line
		if 'DATA SEQUENCE' in line and len(sequence) != 0:
				sequence2 = line.split(' ')
				sequence2 = sequence2[2:]
				sequence2 = "".join(sequence2)
				sequence  = sequence + sequence2

		
		
		# EXTRACTING THE SEQUENCE FROM THE 'DATA SEQUENCE LINE'
		if 'DATA SEQUENCE' in line and len(sequence) == 0:
				sequence = line.split(' ')
				sequence = sequence[2:]
				sequence = "".join(sequence)
		# TAke chemical shifts and strore into dictionary 
		lines = line.split(' ')
		try:

			adjusted_chemical_shift_dict[lines[0]+','+lines[2]]=lines[3]
		except IndexError:
			continue
	

#----------------------------
# Storing sequence information into a dictionrary
# Residue number is the key, residue the value
# Add terminal to each
#----------------------------

for k,residue in enumerate(sequence):
		sequence_dict[str(k+1)            ] = residue
		sequence_dict['0'                 ] = 'TERMINAL'
		sequence_dict[str(len(sequence)+1)] = 'TERMINAL'



#---------------------------------------------------------------------	
# Gaining the phi and psi angles from TALOS-N results 'pred.tab'
# Store in a dictionary based on residue number (storing all residues)
# If TALOS-N does not make a prediction store as 'NaN'
#---------------------------------------------------------------------	

phi_dict      = {}
psi_dict      = {}
phi_dict['0'] = 'NaN'
psi_dict['0'] = 'NaN'

phi_dict[str(len(sequence)+1)] = 'N/A'
psi_dict[str(len(sequence)+1)] = 'N/A'
get           = open(path+'/pred.tab','r')

for line in get:
		line  = " ".join(line.split())
		lines = line.split(' ')
		if 'REMARK' in line:
			continue 
		if 'VARS' in line:
			continue
		if 'FORMAT' in line:
			continue
		if 'DATA' in line:
			continue  
		if len(line) == 0:
			continue


		### YOU CAN  CHOOSE TO COMMENT OUT HERE WHAT TALOS_N BACKBONE PREDICTIONS THAT YOU IGNORE
		if lines[10] != 'Strong' and lines[10] != 'Generous':
 			phi_dict[lines[0]] = 'N/A'
			psi_dict[lines[0]] = 'N/A'
			if lines[1] == 'c':
				print 'WARNING: Backbone dihedral angles for', lines[0], 'are not predicted and must be input manually'
		
		if lines[1] == 'c' and lines[10] == 'Generous':
			print 'WARNING: Backbone dihedral angles for', lines[0], 'is a Generous prediction. Should be checked against structure'		
		if lines[10] == 'Strong' or lines[10] == 'Generous': 
			phi_dict[lines[0]] = lines[2]
			psi_dict[lines[0]] = lines[3]
get.close()   

#---------------------------------------------------------------------	
#### SECONDARY STRUCTURE
#---------------------------------------------------------------------	

ss_dict      = {}
ss_dict['0'] = 'NaN'
ss_dict[str(len(sequence)+1)] = 'N/A'
get           = open(path+'/predSS.tab','r')

for line in get:
		line  = " ".join(line.split())
		lines = line.split(' ')
		if 'REMARK' in line:
			continue 
		if 'VARS' in line:
			continue
		if 'FORMAT' in line:
			continue
		if 'DATA' in line:
			continue  
		if len(line) == 0:
			continue

		ss_dict[lines[0]] = lines[8]
#---------------------------------------------------------------------------------
#All of the shift and angle information has been stored in dictionaries
#Can now go through each connectivity in the previously stored 'Connectivity list'
#Start to assemble information
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------

#								REQUIRED FORMAT FOR DATABASE

# Nuclei order = [N,HA,C,CA,CB,HN]

# Peptide name: Chain: Cys1: Chain: 
# Cys2: Cys1 Nuclei: Cys2 Nuclei:
# X1,X2,X3,X2',X1'
# Cys1 Phi,Psi,X1
# Cys1 secondary structure  
# Cys2 Phi,Psi,X1

#  NEIGHBOURING INFORMATION For Cys1
# Cys1_before vdw
#---------------------------------------------------------------------------------

average_dict       = {}
average_dict['HA'] = 0.117
average_dict['HN'] = 0.067
average_dict['CA'] = 0.453
average_dict['CB'] = 1.422
average_dict['N' ] = 1.307
#For each connectivity
for connectivity in connectivity_list:
			connectivity = [int(x) for x in connectivity]
			#---------------------------------------------
			#By default, the larger residue number is Cys1
			#---------------------------------------------
			connectivity = sorted(connectivity, key=int, reverse=True) 
			cys1         = str(connectivity[0])
			cys2         = str(connectivity[1])
			cys1_before  = str(int(cys1)-1    )
			cys1_after   = str(int(cys1)+1    )
			cys2_before  = str(int(cys2)-1    )
			cys2_after   = str(int(cys2)+1    )

			#------------------------------------------------------------------------
			# CHEMICAL SHIFTS FOR CYSTEINES 
			#------------------------------------------------------------------------
			nuclei_list = ['HA','CB','CA','N','HN']


			#------------------------------------------------------------------------
			# IF CHEMICAL SHIFT ISN"T ASSIGNED RETURN AS ZERO
			#------------------------------------------------------------------------

			def chemical_shift(residue,nuclei):
				try:
					value = adjusted_chemical_shift_dict[residue+','+nuclei]
				except KeyError:
					print 'Warning: Unassigned chemical shift Replaced with average', residue, nuclei
					value =  average_dict[nuclei]
				return (value)

			#------------------------------------------------------------------------
			# STORE IN DICT
			#------------------------------------------------------------------------

			for nuclei in nuclei_list:

				cysteine_dict['Cys1_'       +nuclei].append(chemical_shift(cys1,nuclei       ))
				cysteine_dict['Cys2_'		+nuclei].append(chemical_shift(cys2,nuclei       ))
			#------------------------------------------------------------------------
			# Adding in the Backbones
			#------------------------------------------------------------------------
			cysteine_dict['Cys1_phi'       ].append(phi_dict[cys1		])
			cysteine_dict['Cys1_psi'       ].append(psi_dict[cys1		])
			cysteine_dict['Cys2_phi'       ].append(phi_dict[cys2		])
			cysteine_dict['Cys2_psi'       ].append(psi_dict[cys2		])
			cysteine_dict['Cys1_SS'        ].append(ss_dict[cys1		])
			cysteine_dict['Cys2_SS'        ].append(ss_dict[cys2		])

			#------------------------------------------------------------------------
			# Adding in residues before
			#------------------------------------------------------------------------
			
			cysteine_dict['Cys1_before_residue'].append(sequence_dict[cys1_before])
			cysteine_dict['Cys1_after_residue' ].append(sequence_dict[cys1_after ])
			cysteine_dict['Cys2_before_residue'].append(sequence_dict[cys2_before])
			cysteine_dict['Cys2_after_residue' ].append(sequence_dict[cys2_after ])
			
			#------------------------------------------------------------------------
			# Add in PDB and residues
			#------------------------------------------------------------------------
			cysteine_dict['PDB'].append(peptide)
			cysteine_dict['residue_1'].append(cys1)
			cysteine_dict['residue_2'].append(cys2)

#------------------------------------------------------------------------
# ADD all stored in the cysteine dict to DF
#------------------------------------------------------------------------

for key in cysteine_dict:
	xx = np.array(cysteine_dict[key])
	df[key]=	(xx)


#----------------------------------
# Add in VDW_ RADI OF NEIGHBOURS 
#----------------------------------
#------------------------------------------------------------------------
# VDW Radi dictionray, to be used as 
#------------------------------------------------------------------------
vdw_radi_dict = {}
vdw_radi_dict['A'       ] = 67
vdw_radi_dict['R'       ] = 148
vdw_radi_dict['N'       ] = 96
vdw_radi_dict['D'       ] = 91
vdw_radi_dict['C'       ] = 86
vdw_radi_dict['c'       ] = 86
vdw_radi_dict['Q'       ] = 114
vdw_radi_dict['E'       ] = 109
vdw_radi_dict['G'       ] = 48
vdw_radi_dict['H'       ] = 118
vdw_radi_dict['I'       ] = 124
vdw_radi_dict['L'       ] = 124
vdw_radi_dict['K'       ] = 135
vdw_radi_dict['M'       ] = 124
vdw_radi_dict['F'       ] = 135
vdw_radi_dict['P'       ] = 90
vdw_radi_dict['S'       ] = 73
vdw_radi_dict['T'       ] = 93
vdw_radi_dict['W'       ] = 163
vdw_radi_dict['Y'       ] = 141
vdw_radi_dict['V'       ] = 105
vdw_radi_dict['TERMINAL'] = 0


def vdw_radi(residue):
	vdw_radi_value = vdw_radi_dict[residue]
	return (vdw_radi_value)

#------------------------------------------------------------------------
# Change the residue before to a vdw radi value
#------------------------------------------------------------------------	

df['Cys1_before_residue_VDW']  = df['Cys1_before_residue'].apply(vdw_radi)
df['Cys2_before_residue_VDW']  = df['Cys2_before_residue'].apply(vdw_radi)



#------------------------------------------------------------------------
# Convert SS into array
#------------------------------------------------------------------------	
def ss_array(ss):
		
		secondary_structure_list = [
			'H',
			'E',
			'L'] 
		#print ss
		ss_hot_array 									 = [0. for _ in range(3)]  
		ss_hot_array[secondary_structure_list.index(ss)] = 1
		ss_hot_array                                     = [str(x) for x in ss_hot_array]
		ss_hot_array                                     = ",".join(ss_hot_array)
		return (ss_hot_array)

df['Cys1_SS_array']  = df['Cys1_SS'].apply(ss_array)
df['Cys2_SS_array']  = df['Cys2_SS'].apply(ss_array)

#----------------------------------
# SPLIT INTO INDIVIDUAL HEMI-CYSTEINES, WITH INFORMATION OF CORRESPONDING HEMI-CYSTEINE REQUIRED FOR DISH INPUTS
#----------------------------------

new_columns = ["PDB"				
"residue_number"	
"hemi_number"	
"N"					
"HN"				
"HA"				
"CA"				
"CB"				
"psi"				
"phi"
"ss_array"				
"before_vdw_radi"
"hemi_CA"			
"hemi_CB"			
"hemi_psi"			
"ss"
"ss_array"]	

df_  = pd.DataFrame(index = None,columns = new_columns)
df_  =  pd.DataFrame({                                                
			"PDB"				:df['PDB'                    ],
			"residue_number"	:df["residue_1"              ],
			"hemi_number"	    :df["residue_2"              ],
			"N"					:df["Cys1_N"                 ],
			"HN"				:df["Cys1_HN"                ],
			"HA"				:df["Cys1_HA"                ],
			"CA"				:df["Cys1_CA"                ],
			"CB"				:df["Cys1_CB"                ],
			"psi"				:df["Cys1_psi"               ],
			"phi"				:df["Cys1_phi"               ],
			"before_vdw_radi"	:df["Cys1_before_residue_VDW"],
			"ss"				:df["Cys1_SS" 				 ],
			"ss_array"			:df["Cys1_SS_array" 		 ],
			"hemi_CA"			:df["Cys2_CA"                ],
			"hemi_CB"			:df["Cys2_CB"                ],
			"hemi_N"			:df["Cys2_N"                ],
			"hemi_psi"			:df["Cys2_psi"               ]})
	
	
	
df__ = pd.DataFrame(index = None,columns = new_columns)
df__  =  pd.DataFrame({                                                
			"PDB"				:df['PDB'                    ],
			"residue_number"	:df["residue_2"              ],
			"hemi_number"		:df["residue_1"              ],
			"N"					:df["Cys2_N"                 ],
			"HN"				:df["Cys2_HN"                ],
			"HA"				:df["Cys2_HA"                ],
			"CA"				:df["Cys2_CA"                ],
			"CB"				:df["Cys2_CB"                ],
			"psi"				:df["Cys2_psi"               ],
			"phi"				:df["Cys2_phi"               ],
			"before_vdw_radi"	:df["Cys2_before_residue_VDW"],
			"ss"				:df["Cys2_SS" 				 ],
			"ss_array"			:df["Cys2_SS_array" 		 ],
			"hemi_CA"			:df["Cys1_CA"                ],
			"hemi_CB"			:df["Cys1_CB"                ],
			"hemi_N"			:df["Cys1_N"                ],
			"hemi_psi"			:df["Cys1_psi"               ]})
	

df_split = df_.append(df__, ignore_index = True)


df_split = df_split[["PDB",				
"residue_number",	
"hemi_number",	
"N",					
"HN",				
"HA",				
"CA",				
"CB",				
"psi",				
"phi",				
"before_vdw_radi",
"ss",
"ss_array",
"hemi_CA",			
"hemi_CB",			
"hemi_N",			
"hemi_psi"]]




df_split = df_split.loc[df_split['before_vdw_radi'] != 0]


df_split['residue_number'] = df_split['residue_number'].astype(int)
df_split = df_split.sort_values("residue_number")
print df_split
df_split                        = df_split.replace(r'\s+', np.nan, regex=True)
df_split.to_csv(path+'/DISH_inputs.csv',index = False)
			# print df 

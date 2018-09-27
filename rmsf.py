#!/home/rbdavid/bin/python
# ----------------------------------------
# USAGE:

# ./rmsd.ref.py config_file

# ----------------------------------------
# PREAMBLE:

import sys
import importlib
import numpy as np
import MDAnalysis

config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS: 

necessary_parameters = ['reference_pdb','analysis_pdb','traj_loc','start','end','data_output_filename','selection_file','distance_functions_file']
all_parameters = ['reference_pdb','analysis_pdb','traj_loc','start','end','data_output_filename','selection_file','distance_functions_file','wrapping_boolean','wrapping_selection','alignment_selection','analysis_selection','protein_selection','summary_boolean','summary_filename','selection_output_filename']
def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['wrapping_boolean'] = False
	parameters['wrapping_selection'] = None
	parameters['alignment_selection'] = 'protein'
	parameters['analysis_selection'] = 'protein'
        parameters['protein_selection'] = 'all'
	parameters['summary_boolean'] = False 
	parameters['summary_filename'] = None
	parameters['selection_output_filename'] = 'selections.txt'

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

def summary(summary_filename):
	with open(summary_filename,'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			f.write('%s ' %(sys.argv[i]))
		f.write('\n\n')
		f.write('Parameters used:\n')
                for key, value in parameters.iteritems():
                        if key == '__builtins__':
                                continue
                        else:
			        f.write("%s = '%s'\n" %(key,value))
		f.write('\n')

def main():

        # ----------------------------------------
        # LOAD IN THE REFERENCE UNIVERSE
        ref = MDAnalysis.Universe(parameters['reference_pdb'])
        ref_align = ref.select_atoms(parameters['alignment_selection'])
        ref_all = ref.select_atoms('all')
        ref_analysis = ref.select_atoms(parameters['analysis_selection'])
        ref_all.translate(-ref_align.center_of_mass())
        pos0 = ref_align.positions
        
        # ----------------------------------------
        # LOAD IN THE ANALYSIS UNIVERSE
        u = MDAnalysis.Universe(parameters['analysis_pdb'])
        u_align = u.select_atoms(parameters['alignment_selection'])
        u_all = u.select_atoms('all')
        u_analysis = u.select_atoms(parameters['analysis_selection'])
        nResidues_analysis = u_analysis.n_residues

        if u_align.n_atoms != ref_align.n_atoms:
        	print 'Alignment atom selections do not have the same number of atoms.'
        	sys.exit()

        # ----------------------------------------
        # MAKE ATOM SELECTIONS FOR SUBSEQUENT RMSF ANALYSIS
        selection_list = []
        nAtoms = []
        reference_positions = []
        count = 0
        nResidues_analysis_range = range(nResidues_analysis)
        with open(parameters['selection_output_filename'],'w') as f:
                for i in nResidues_analysis_range:
                        temp_resname = u_analysis.residues[i].resname
                        temp_resid = u_analysis.residues[i].resid
        		num_selections_made = make_selections(u,ref,temp_resname,temp_resid,f,selection_list,nAtoms,reference_positions,count,parameters['protein_selection'])
                        count += num_selections_made
        
        # ----------------------------------------
        # RMSF ANALYSIS OF ATOM SELECTIONS
        trajectory_list = [parameters['traj_loc']%(i) for i in range(parameters['start'],parameters['end']+1)]

        if parameters['wrapping_boolean']:
                RMSF_calc(u,u_all,u_align,pos0,selection_list,nAtoms,reference_positions,trajectory_list,parameters['data_output_filename'],wrapping_boolean = True,wrapping_selection = parameters['wrapping_selection'])
        else:
                RMSF_calc(u,u_all,u_align,pos0,selection_list,nAtoms,reference_positions,trajectory_list,parameters['data_output_filename'])

        if parameters['summary_boolean']:
                summary(parameters['summary_filename'])

# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOADING IN NECESSARY FUNCTIONS FROM MODULE FILES
make_selections = importlib.import_module(parameters['selection_file'].split('.')[0],package=None).make_selections
RMSF_calc = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).RMSF_calc
if parameters['wrapping_boolean']:
        wrapping = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).wrapping

# ----------------------------------------
# MAIN
if __name__ == '__main__':
	main()


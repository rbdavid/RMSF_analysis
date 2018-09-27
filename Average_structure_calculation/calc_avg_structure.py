#!/home/rbdavid/bin/python
# ----------------------------------------
# USAGE:

# ./calc_avg_structure.py calc_avg_structure.config

# ----------------------------------------
# PREAMBLE:

import sys
import importlib
import MDAnalysis

config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS:

necessary_parameters = ['distance_functions_file','pdb','traj_loc','start','end','alignment_selection','analysis_selection','average_structure_file_name']
all_parameters = ['distance_functions_file','pdb','traj_loc','start','end','alignment_selection','analysis_selection','average_structure_file_name','wrapping_boolean','summary_boolean','summary_file_name']

def config_parser(config_file):	
        """ Function to take config file and create/fill the parameter dictionary (created before function call). 
        
        Usage: 
            parameters = {}     # initialize the dictionary to be filled with keys and values
            config_parser(config_file)

        Arguments:
            config_file: string object that corresponds to the local or global position of the config file to be used for this analysis.

        """
        
        # NECESSARY PARAMETERS ARE INITIALIZED IN DICTIONARY WITH EMPTY STRINGS:
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''
	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['wrapping_boolean'] = False
        parameters['summary_boolean'] = False
        parameters['summary_file_name'] = None
	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
        # TESTING IF ANY PARAMETER HAS BEEN LEFT EMPTY:
        for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

def summary(summary_filename):
        """ Function to create a text file that holds important information about the analysis that was just performed. Outputs the version of MDAnalysis, how to rerun the analysis, and the parameters used in the analysis.

        Usage:
            summary(summary_filename)

        Arguments:
            summary_filename: string object of the file name to be written that holds the summary information.

        """
                
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
                        #if type(value) == int or type(value) == float:
			#        f.write("%s = %s\n" %(key,value))
                        else:
			        f.write("%s = '%s'\n" %(key,value))
		f.write('\n')

def main():

        trajectory_list = [parameters['traj_loc']%(i) for i in range(parameters['start'],parameters['end']+1)]

        if parameters['wrapping_boolean']:
                average_structure_calc(parameters['pdb'],parameters['alignment_selection'],parameters['analysis_selection'],trajectory_list,parameters['average_structure_file_name'],wrapping_boolean = True)
        else:
                average_structure_calc(parameters['pdb'],parameters['alignment_selection'],parameters['analysis_selection'],trajectory_list,parameters['average_structure_file_name'])

        if parameters['summary_boolean']:
                summary(parameters['summary_file_name'])

# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOADING IN NECESSARY FUNCTIONS FROM MODULE FILES
average_structure_calc = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).average_structure_calc
if parameters['wrapping_boolean']:
        wrapping = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).wrapping

# ----------------------------------------
# MAIN
if __name__ == '__main__':
	main()


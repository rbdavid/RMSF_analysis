
# USAGE:
# from distance_functions import *

# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np
from numpy.linalg import *

sums = np.sum
square = np.square
zeros = np.zeros
dot_prod = np.dot

# SUBROUTINES/FUNCTIONS:

def RMSD(x,y,n):
	""" Calculates the Root Mean Squared Distance between two arrays of the same size

	Usage: rmsd = RMSD(x,y,n)

	Arguments:
	x, y: numpy arrays with the same shape (n X 3)
	n: number of particles being summed over; ex: number of atoms in the atom selection being analyzed;
		if n = 1, this function calculates the distance between x and y arrays

	"""
	
	return (sums(square(x-y))/n)**0.5	# the MSD should never be negative, so using **0.5 rather than np.sqrt is safe

def MSD(x,y,n):
	""" Calculates the Mean Squared Distance between two arrays of the same size

	Usage: msd = MSD(x,y,n)

	Arguments:
	x, y: numpy arrays with the same shape
	n: number of particles being summed over; ex: number of atoms in the atom selection being analyzed;
		if n = 1, this function calculates the distance squared between x and y arrays

	"""

	return sums(square(x-y))/n

def wrapping(x,dim):
	""" Calculates the translation matrix needed to wrap a particle back into the original periodic box
	
	Usage: t = wrapping(x,dim)

	Arguments:
	x: a numpy array of size (3) that corresponds to the xyz coordinates of an ATOM/COM/COG of a residue
	dim: a numpy array of size (3) that holds the xyz dimensions of the periodic box at that timestep

	"""
	
	t = zeros(3)
	dim2 = dim/2.
	for i in range(3):
		if (x[i]<-dim2[i]) or (x[i]>dim2[i]):
			t[i] = -dim[i]*round(x[i]/dim[i])
	return t

def euclid_dist(x,y):
	""" Calculates the Euclidian Distance between two arrays of the same size
	Usage: dist,dist2 = euclid_dist(x,y)
		
	Arguments:
	x, y: numpy arrays with the same size
	"""
	
	dist2 = sums(square(x-y))
	dist = dist2**0.5	# the MSD should never be negative, so using **0.5 rather than np.sqrt is safe
	return dist, dist2

def average_structure_calc(pdb,alignment_selection,analysis_selection,traj_list,avg_structure_file_name, convergence_threshold = 1E-5, maximum_num_iterations = 100, wrapping_boolean = False):
        """
        """

        # ----------------------------------------
        # LOAD IN AND CREATE ATOM SELECTIONS IN THE ANALYSIS UNIVERSE OBJECT
        u = MDAnalysis.Universe(pdb)
        u_all = u.select_atoms('all')
        u_align = u.select_atoms(alignment_selection)  # MDAnalysis atom selection string formatting required.
        u_analysis = u.select_atoms(analysis_selection)  # MDAnalysis atom selection string formatting required.
        nAtoms_align = u_align.n_atoms
        nAtoms_analysis = u_analysis.n_atoms

        # ----------------------------------------
        # ANALYZE TRAJECTORIES TO COLLECT THE NECESSARY POSITION DATA
        all_pos_align = []
        all_pos_analysis = []
        nSteps = 0
        nResidues_range = range(u_analysis.n_residues)
        for i in traj_list:
                print 'Loading trajectory', i
                u.load_new(i)
                nSteps += len(u.trajectory)
                for ts in u.trajectory:
                        dimensions = u.dimensions[:3]
                        u_all.translate(-u_align.center_of_mass())
                        
                        if wrapping_boolean:
                                for j in nResidues_range:
                                        COM = u_analysis.residues[j].center_of_mass()
                                        t = wrapping(COM,dimensions)
                                        u_analysis.residues[j].atoms.translate(t)
                        
                        all_pos_align.append(u_align.positions)
                        all_pos_analysis.append(u_analysis.positions)
        
        print 'Analyzed', nSteps, 'frames. Does this match up with expectation?'

        all_pos_align = np.array(all_pos_align)
        all_pos_analysis = np.array(all_pos_analysis)

        avg_pos_align = np.sum(all_pos_align,axis=0)/nSteps
        avg_pos_analysis = np.sum(all_pos_analysis,axis=0)/nSteps

        # ----------------------------------------
        # ITERATIVE ALIGNMENT TO AVERAGE ALIGNMENT POSITIONS
        iteration = 0 
        residual = convergence_threshold + 9999.
        nSteps_range = range(nSteps)
        print 'Beginning the iterative process of aligning to the average alignment positions, calculating new positions, and recalculating the average positions'
        while residual > convergence_threshold and iteration < maximum_num_iterations:
                temp_avg_pos_align = np.zeros((nAtoms_align,3),dtype=np.float32)
                temp_avg_pos_analysis = np.zeros((nAtoms_analysis,3),dtype=np.float32)

                for i in nSteps_range:
                        R, d = rotation_matrix(all_pos_align[i,:,:],avg_pos_align)      # calculate the rotation matrix (and distance) between frame i's alignment postions to the average alignment positions
                        all_pos_align[i,:,:] = dot_prod(all_pos_align[i,:,:],R.T)       # take the dot product between frame i's alignment positions and the calculated rotation matrix; overwrite frame i's positions with the rotated postions
                        all_pos_analysis[i,:,:] = dot_prod(all_pos_analysis[i,:,:],R.T) # take the dot product between frame i's analysis positions and the calculated rotation matrix; overwrite frame i's positions with the rotated postions
                        temp_avg_pos_align += all_pos_align[i,:,:]          # running sum of alignment positions to calculate a new average
                        temp_avg_pos_analysis += all_pos_analysis[i,:,:]    # running sum of analysis positions to calculate a new average
               
                temp_avg_pos_align /= nSteps
                temp_avg_pos_analysis /= nSteps
                residual = RMSD(avg_pos_align,temp_avg_pos_align,nAtoms_align)
                analysis_RMSD = RMSD(avg_pos_analysis,temp_avg_pos_analysis,nAtoms_analysis)
                iteration += 1
                avg_pos_align = temp_avg_pos_align
                avg_pos_analysis = temp_avg_pos_analysis
                print 'Iteration ', iteration, ': RMSD btw alignment landmarks: ', residual,', RMSD btw analysis atoms: ', analysis_RMSD
        
        print 'Finished calculating the average structure using the iterative averaging approach. Outputting the average structure to a pdb now.'
        # ----------------------------------------
        # LOAD IN AND CREATE ATOM SELECTIONS IN THE RESULTS UNIVERSE OBJECT
        avg = MDAnalysis.Universe(pdb)
        avg_analysis = avg.select_atoms(analysis_selection)
        avg_analysis.positions = avg_pos_analysis
        avg_analysis.write(avg_structure_file_name)

        return avg_pos_analysis, all_pos_analysis

def dist_matrix_calc(pdb,atom_selections,traj_loc,start,end,system_descriptor,ignore_n_nearest_neighbors=0,step=1):
        """
        """

        # ----------------------------------------
        # CREATE OUTPUT FILE NAMING VARIABLES
        node_output_filename = system_descriptor + '.nodes.txt'
        selection_output_filename = system_descriptor + '.selections.txt'
        data_output_filename = system_descriptor + '.col_var.dat'

        # ----------------------------------------
        # LOAD IN AND CREATE ATOM SELECTIONS IN THE ANALYSIS UNIVERSE OBJECT
        u = MDAnalysis.Universe(pdb)
        col_var_selections = u.select_atoms(atom_selections)  # MDAnalysis atom selection string formatting required.
        nNodes = col_var_selections.n_atoms     # assumes each col var is a distance between a pair of atoms, i and j. 
        nNodes_range = range(nNodes)
	i_max = nNodes-1-ignore_n_nearest_neighbors   # max value of atom looping index i; to be used multiple times so why calculate it multiple times
        boolean_matrix = np.full((nNodes,nNodes),False)  # 2D matrix of False values; elements of this matrix will be set to True as we loop over our i,j atom pairs. This will be used later for the plotting of 1D collective variable vectors onto the respective 2D atom pair matrix.

        count = 0
        with open(selection_output_filename,'w') as W, open(node_output_filename,'w') as Y:
                Y.write('# Node description: atom name, index, residresname\n')
                for i in nNodes_range:
                        Y.write('Nodes: %s %d %s%d\n'%(col_var_selections[i].name,col_var_selections[i].index+1,col_var_selections[i].resname,col_var_selections[i].resid))
                
                W.write('# Collective variable description: atom name, index, residresname to atom name, index, residresname\n')
		for i in nNodes_range[:-1-ignore_n_nearest_neighbors]:
                        for j in nNodes_range[i+1+ignore_n_nearest_neighbors:]:
                            boolean_matrix[i,j] = True
                            W.write('Collective Variable: %s %d %s%d   to    %s %d %s%d\n'%(col_var_selections[i].name,col_var_selections[i].index+1,col_var_selections[i].resname,col_var_selections[i].resid,col_var_selections[j].name,col_var_selections[j].index+1,col_var_selections[j].resname,col_var_selections[j].resid)) 
                            count += 1
        
        print 'Number of nodes:', nNodes,', while skipping ', ignore_n_nearest_neighbors, ' nearest neighbors, creates ', count, 'number of collective variables to be analyzed.'

        # ----------------------------------------
        # TRAJECTORY ANALYSIS
        start = int(start)
        end = int(end)
        print 'Beginning trajectory analysis'
        with open(data_output_filename,'w') as W:
                while start <= end:
	                print 'Loading trajectory ', start, ', located at ', traj_loc%(start)
                        u.load_new(traj_loc%(start))
                        for ts in u.trajectory[::step]:
                                temp_positions = col_var_selections.positions
                                for i in nNodes_range[:-1-ignore_n_nearest_neighbors]:
                                        for j in nNodes_range[i+1+ignore_n_nearest_neighbors:]:
                                                dist,dist2 = euclid_dist(temp_positions[i],temp_positions[j])
                                                W.write('%f '%(dist))
                                W.write('\n')
                        start += 1

        print 'Finished analyzing trajectory and outputting raw data. Onto calculating the covariance and correlation matrices for the collective variables analyzed.'
        return data_output_filename, boolean_matrix

def calc_collective_variable_boolean_matrix(nColVars,ignore_n_nearest_neighbors=0):
        """
        """
        nNodes = int(np.round(np.max(np.roots([1,-1,1-nColVars*2]))))+ignore_n_nearest_neighbors
        nNodes_range = range(nNodes)
        boolean_matrix = np.full((nNodes,nNodes),False)  # 2D matrix of False values; elements of this matrix will be set to True as we loop over our i,j atom pairs. This will be used later for the plotting of 1D collective variable vectors onto the respective 2D atom pair matrix.
	for i in nNodes_range[:-1-ignore_n_nearest_neighbors]:
                for j in nNodes_range[i+1+ignore_n_nearest_neighbors:]:
                        boolean_matrix[i,j] = True

        print boolean_matrix
        print 'Number of nodes (aka atoms) in the original atom selection is', nNodes, 'with', ignore_n_nearest_neighbors, 'nearest neighbors in this atom selection being ignored in the collective variable space. This creates', nColVars, 'collective variables. Do these numbers check out with the past analysis?'

        return boolean_matrix



import MDAnalysis
import sys

# ----------------------------------------
# STANDARD RESIDUES:
protein = ['ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP','SER','THR','ASN','GLN','ARG','LYS','ASP','GLU','GLY','PRO','CYS','HIE','HID','HIP','HIS']
nucleic = ['A5','A3','A','G5','G3','G','C5','C3','C','T5','T3','T','U5','U3','U']
triphosphate = ['atp','adp','PHX']
other = ['MG']

# ----------------------------------------
# HOMEMADE ATOM SELECTION STRINGS FOR THE STANDARD RESIDUES:
sugar = "name C5' C4' O4' C1' C3' C2' O2' " + " C5* C4* O4* C1* C3* O3* C2* O2* "		# NO HYDROGENS; DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...; the atoms with * are found in triphosphates;
sugar_5= sugar + " O5'"		# NO HYDROGENS
sugar_3= sugar + " O3' "	# NO HYDROGENS
base = 'name N9 C8 N7 C5 C6 N6 N1 C2 N3 C4 O6 N4 C2 O2 O4'	# NO HYDROGENS; selection string that will select all appropriate atoms for any of the nucleic residues...
a_phos = 'name O5* O2A O1A PA O3A'
b_phos = 'name PB O1B O2B O3B'
g_phos = 'name PG O1G O2G O3G'
inorg_phos = 'name P O1 O2 O3 O4'	# NO HYDROGENS

# ----------------------------------------
# FUNCTION USED TO MAKE ANY OF THE HOMEMADE ATOM SELECTIONS FOR THE STANDARD RESIDUES

def make_selections(analysis_universe,ref_universe,resname,resid,output_file,selection_list,nAtoms,ref_pos,count,protein_selection):
	"""A function that takes in a residue name and creates a non-standard MDAnalysis atom selection
	
	Usage: make_selection(........)
	
	Arguments:
		analysis_universe: MDAnalysis Universe object to be used as the analysis universe.
		reference_universe: MDAnalysis Universe object to be used as the reference universe.
		resname: string of the residue name;
		resid: int of the residue ID number;
		output_file: file object that is to be written to;
        Returns:
                n_selections_made
        """
        

        # ----------------------------------------
        # CREATING THE PROTEIN SELECTIONS 
        if resname in protein:
		# CREATING THE SELECTION FOR THE PROTEIN RESIDUES
                sel_string = 'resname %s and resid %d and %s' %(resname,resid,protein_selection)
                u_temp = analysis_universe.select_atoms(sel_string)
                selection_list.append(u_temp)
                nAtoms.append(u_temp.n_atoms)
                ref_temp = ref_universe.select_atoms(sel_string)
                ref_pos.append(ref_temp.positions)
		if u_temp.n_atoms != ref_temp.n_atoms:
			print 'Number of atoms do not match for selection %d, %s, %s' %(count,resname,sel_string)
			sys.exit()
		output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
		count +=1
                return 1

	# ----------------------------------------
	# CREATING THE NUCLEIC SELECTIONS
	elif resname in nucleic:
		# CREATING THE SELECTION FOR THE BASE OF NUCLEIC RESIDUES
		sel_string = 'resname %s and resid %d and %s' %(resname,resid,base)
		u_temp = analysis_universe.select_atoms(sel_string)
		selection_list.append(u_temp)
		nAtoms.append(u_temp.n_atoms)
		ref_temp = ref_universe.select_atoms(sel_string)
		ref_pos.append(ref_temp.positions)
		if u_temp.n_atoms != ref_temp.n_atoms:
			print 'Number of atoms do not match for selection %d, %s, %s' %(count,resname,sel_string)
			sys.exit()
		output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
		count +=1

		# CREATING THE SLECTION FOR THE SUGAR OF NUCLEIC RESIDUES
		if resname in ['A5','G5','C5','T5','C5']:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar_5)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %d, %s, %s' %(count,resname,sel_string)
				sys.exit()
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1
			return 2

		elif resname in ['A3','U3','C3','G3']:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar_3)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %d, %s, %s' %(count,resname,sel_string)
                                sys.exit()
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1

		else:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %d, %s, %s' %(count,resname,sel_string)
				sys.exit()
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1

		# CREATING THE SLECTION FOR THE PHOSPHATE OF NUCLEIC RESIDUES
		sel_string = "(resname %s and resid %s and name P OP1 OP2 O5') or (resid %s and name O3')" %(resname,resid,analysis_universe.residues[resid-1].resid) 
		u_temp = analysis_universe.select_atoms(sel_string)
		selection_list.append(u_temp)
		nAtoms.append(u_temp.n_atoms)
		ref_temp = ref_universe.select_atoms(sel_string)
		ref_pos.append(ref_temp.positions)
		output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
		count += 1
		return 3

	# ----------------------------------------
	# CREATING THE TRIPHOSPHATE ATOM SELECTIONS
	elif resname in triphosphate:
		if resname in ['atp','adp']:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,base)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %d, %s, %s' %(count,resname,sel_string)
				sys.exit()
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %d, %s, %s' %(count,resname,sel_string)
				sys.exit()
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1

		if resname == 'atp':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,a_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,b_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,g_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1
			return 5

		elif resname == 'adp':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,a_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,b_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1
			return 4

		elif resname == 'PHX':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,inorg_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
			count +=1
			return 1

	# ----------------------------------------
	# CREATING ANY REMAINING ATOM SELECTIONS...
	elif resname in other:
		sel_string = 'resname %s and resid %d' %(resname,resid)
		u_temp = analysis_universe.select_atoms(sel_string)
		selection_list.append(u_temp)
		nAtoms.append(u_temp.n_atoms)
		ref_temp = ref_universe.select_atoms(sel_string)
		ref_pos.append(ref_temp.positions)
		output_file.write('%d   %s   %s   %d\n' %(count,resname,sel_string,u_temp.n_atoms))
		count +=1
		return 1

# ----------------------------------------



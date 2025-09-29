import numpy as np
from scipy.spatial.distance import pdist, squareform
from collections import Counter
from Bio.PDB import PDBParser, Select, PDBIO
import freesasa
import os
import pandas as pd
from pathlib import Path

from utils import load_param, calculate_pairwise_distances_intramolecular, calculate_pairwise_distances_intermolecular

class VDWAnalyzer:
    def __init__(self, structure, ligand_chains, receptor_chains, short_threshold=5, long_threshold=10):
        self.structure = structure
        self.ligand_chains = ligand_chains
        self.receptor_chains = receptor_chains
        self.calculate_distances = calculate_pairwise_distances_intermolecular
        self.short_threshold = short_threshold
        self.long_threshold = long_threshold
        
    def __call__(self):
        
        self.coordinates_1 = self.extract_coordinates(self.receptor_chains)
        self.coordinates_2 = self.extract_coordinates(self.ligand_chains)
        self.distances = self.calculate_distances(self.coordinates_1, self.coordinates_2, return_dict=True)
        self.pairs = np.asarray(list(self.distances.keys()))
        self.radii_1 = self.pairs[:,0]
        self.radii_2 = self.pairs[:,1]
        self.calculate_repulsive_london()
        self.calculate_vdw_total()
        self.calculate_lennard_jones()
        self.analyze_electrostatic_interactions()
        
    def extract_coordinates(self, target_chains):
        coordinates = []
        for model in self.structure:
            for chain in model:
                if chain.id not in target_chains:
                    continue
                for residue in chain:
                    for atom in residue:
                        coordinates.append(atom.get_coord())
        return np.array(coordinates)
       
    def calculate_repulsive_london(self):
        """Calculate repulsive interactions and London dispersion forces."""
        distances = np.asarray(list(self.distances.values()))
        allowed = self.radii_1 + self.radii_2
        
        # Avoid division by zero
        distances = np.maximum(distances, 1e-6)
        
        repulsive_mask = distances <= allowed
        self.repulsive = np.sum(repulsive_mask)
        self.london = np.sum(1 / distances)
        
        self.repulsive_adj = self.repulsive / (len(self.radii_1)+len(self.radii_2))
        self.london_adj = self.london / (len(self.radii_1)+len(self.radii_2))
        
    def calculate_vdw_total(self):
        """Calculate the total van der Waals interaction for multiple atom pairs."""
        distances = np.asarray(list(self.distances.values()))
        vdw_radii = self.radii_1 + self.radii_2
        ratio = vdw_radii / distances
        vdw_interactions = np.clip(ratio**8 - 2 * ratio**4, a_min=None, a_max=100)
        self.vdw_total = np.sum(vdw_interactions)
        self.vdw_adj = self.vdw_total/(len(self.radii_1)+len(self.radii_2))
        
    def calculate_lennard_jones(self, epsilon=1.0, sigma=1.0):
        """Calculate Lennard-Jones potential and force."""
        distances = np.asarray(list(self.distances.values()))
        sr6 = (sigma/distances)**6
        self.lj_potential = np.sum(4 * epsilon * (sr6**2 - sr6))
        self.lj_force = np.sum(24 * epsilon * (2 * sr6**2 - sr6) / distances)
        self.lj_potential_adj = self.lj_potential/(len(self.radii_1)+len(self.radii_2))
        self.lj_force_adj = self.lj_force/(len(self.radii_1)+len(self.radii_2))
    
    def _is_charged(self, residue):
        """Identify charged residues including histidine (always treated as potentially positive)."""
        res_name = residue.get_resname().strip().upper()
        
        # Standard charged residues
        pos_charges = {'ARG', 'LYS', 'HIS'}  # Added HIS to positive charges
        neg_charges = {'ASP', 'GLU'}
        
        if res_name in pos_charges:
            return 'positive'
        if res_name in neg_charges:
            return 'negative'
        return None 
    
    def get_charged_residues(self, target_chains):
        """Extract charged residues with their atoms and charge type."""
        charged = []
        for model in self.structure:
            for chain in model:
                if chain.id not in target_chains:
                    continue
                for residue in chain:
                    charge = self._is_charged(residue)
                    if charge:
                        atoms = [atom for atom in residue]
                        charged.append((residue, charge, atoms))
        return charged
    
    def analyze_electrostatic_interactions(self):
        """Calculate electrostatic interactions between charged residues."""
        receptor_charged = self.get_charged_residues(self.receptor_chains)
        ligand_charged = self.get_charged_residues(self.ligand_chains)
        
        counts = {
            'attractive_short': 0,
            'attractive_long': 0,
            'repulsive_short': 0,
            'repulsive_long': 0,
            'his_involved_short': 0,  # Track histidine participation
            'his_involved_long': 0  # Track histidine participation
        }

        for rec_res, rec_charge, rec_atoms in receptor_charged:
            for lig_res, lig_charge, lig_atoms in ligand_charged:
                # Check if histidine is involved
                his_present = (rec_res.get_resname().strip().upper() == 'HIS' or 
                              lig_res.get_resname().strip().upper() == 'HIS')
                
                # Calculate all atom-atom distances between residue pair
                rec_coords = np.array([atom.get_coord() for atom in rec_atoms])
                lig_coords = np.array([atom.get_coord() for atom in lig_atoms])
                distances = self.calculate_distances(rec_coords, lig_coords, return_dict=False)
                
                # Determine interaction type
                if rec_charge != lig_charge:
                    interaction = 'attractive'
                else:
                    interaction = 'repulsive'
                
                # Categorize distances
                short_mask = (distances <= self.short_threshold)
                long_mask = ((distances > self.short_threshold) & 
                            (distances <= self.long_threshold))
                
                counts[f'{interaction}_short'] += np.sum(short_mask)
                counts[f'{interaction}_long'] += np.sum(long_mask)
                if his_present and np.any(short_mask):
                    counts['his_involved_short'] += 1
                if his_present and np.any(long_mask):
                    counts['his_involved_long'] += 1
                    
        
        self.electrostatic_counts = counts    
    
    def get_results(self):
        """Return a dictionary of all calculated results."""
        results =  {
            'repulsive interactions': self.repulsive,
            'repulsive adj': self.repulsive_adj,
            'london dispersion': self.london,
            'london adj': self.london_adj,
            'vdw total': self.vdw_total,
            'vdw adj': self.vdw_adj,
            'lj potential': self.lj_potential,
            'lj padj': self.lj_potential_adj,
            'lj force': self.lj_force,
            'lj fadj': self.lj_force_adj
        }
        results.update(self.electrostatic_counts)
        return results

class HCAnalyzer:
    
    def __init__(self, structure, ligand_chains, receptor_chains):
        
        self.structure = structure
        self.ligand_chains = ligand_chains
        self.receptor_chains = receptor_chains

        self.load_param = load_param
        self.calculate_distances = calculate_pairwise_distances_intermolecular
        
        self.atom_radii = self.load_param('pdb_radii.param')
        
        self.short_dist_threshold = 0.5
        self.long_dist_threshold = 2.0
        
        self.hc_res_atom_map = {
            'ALA': ['CA', 'C', 'CB'],
            'ARG': ['CA', 'C', 'CB', 'CG', 'CD'],
            'ASN': ['CA', 'C', 'CB', 'CG'],
            'ASP': ['CA', 'C', 'CB'],
            'CYS': ['CA', 'C', 'CB', 'SG'],
            'GLN': ['CA', 'C', 'CB', 'CG', 'CD'],
            'GLU': ['CA', 'C', 'CB', 'CG'],
            'GLY': ['CA', 'C'],
            'HIS': ['CA', 'C', 'CB', 'CG', 'CD2', 'CE1'],
            'ILE': ['CA', 'C', 'CB', 'CG1', 'CG2', 'CD1'],
            'LEU': ['CA', 'C', 'CB', 'CG', 'CD1', 'CD2'],
            'LYS': ['CA', 'C', 'CB', 'CG', 'CD', 'CE'],
            'MET': ['CA', 'C', 'CB', 'CG', 'SD', 'CE'],
            'PHE': ['CA', 'C', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'PRO': ['CA', 'C', 'CB', 'CG', 'CD'],
            'SER': ['CA', 'C', 'CB'],
            'THR': ['CA', 'C', 'CB', 'CG2'],
            'TRP': ['CA', 'C', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'TYR': ['CA', 'C', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2'],
            'VAL': ['CA', 'C', 'CB', 'CG1', 'CG2']
        }

    def __call__(self):
        
        self.receptor_data = self.extract_hc_atoms(self.receptor_chains)
        self.ligand_data = self.extract_hc_atoms(self.ligand_chains)
        self.hydrophobic_distances = self.calculate_distances(self.receptor_data['hc_coords'], self.ligand_data['hc_coords'], return_dict=True)
        self.compute_hc(self.receptor_data, self.ligand_data)
        
    def extract_hc_atoms(self, target_chains):
        # Initialize lists to store donor and acceptor information
        hc_atoms = []
        hc_radii = []
        hc_coords = []
        hc_res = []
        
        data_dict = {}
        data_dict['protein_length'] = 0

        # Iterate through all atoms in the structure
        for model in self.structure:
            for chain in model:
                if chain.id not in target_chains:
                    continue
                for residue in chain:
                    res_name = residue.get_resname()
                    res_id = residue.get_id()[1]  # Get residue number
                    for atom in residue:
                        atom_name = atom.get_name()
                        atom_id = atom.get_serial_number()
                        coord = atom.get_coord()

                        if atom_name not in self.atom_radii:
                            continue  # Skip this atom if its radius is not known

                        data_dict['protein_length'] += 1
                        # Check if atom is a potential acceptor
                        if res_name in self.hc_res_atom_map and atom_name in self.hc_res_atom_map[res_name]:
                            hc_atoms.append(atom_id)
                            hc_radii.append(self.atom_radii[atom_name])
                            hc_coords.append(coord)
                            hc_res.append(res_id)

        # Convert lists to numpy arrays
        data_dict['hc_atoms'] = np.array(hc_atoms, dtype=int)
        data_dict['hc_radii'] = np.array(hc_radii, dtype=float)
        data_dict['hc_coords'] = np.array(hc_coords, dtype=float)
        data_dict['hc_res'] = np.array(hc_res, dtype=int)
        
        return data_dict
        
    def compute_hc(self, data_dict_1, data_dict_2):
        
        hydrophobic_pairs = np.asarray(list(self.hydrophobic_distances.keys()))
        r1 = data_dict_1['hc_radii'][hydrophobic_pairs[:,0]]
        r2 = data_dict_2['hc_radii'][hydrophobic_pairs[:,1]]
        
        distances = np.asarray(list(self.hydrophobic_distances.values()))

        # Check that all input arrays have the same length
        if not (len(distances) == len(r1) == len(r2)):
            raise ValueError("All input arrays must have the same length", len(distances), len(r1), len(r2))

        # Calculate the sum of VDW radii
        vdw_radii = r1 + r2

        HC_allowed_1 = vdw_radii + self.short_dist_threshold
        HC_allowed_2 = vdw_radii + self.long_dist_threshold

        # Calculate HC2 scores
        HC2 = np.zeros_like(distances)
        mask1 = distances <= HC_allowed_1
        mask2 = (distances > HC_allowed_1) & (distances <= HC_allowed_2)            

        HC2[mask1] = 1
        HC2[mask2] = (1 / 1.5) * ((vdw_radii[mask2] + 2.0) ** 2 - distances[mask2] ** 2)

        self.HC_total = np.sum(HC2)
        self.HC_adj = self.HC_total/(data_dict_1['protein_length']+data_dict_2['protein_length'])

        #https://pubs.rsc.org/en/content/articlelanding/2014/cp/c4cp01131g
        #six ranges of distances (r1–r6) for each analyzed lattice. 
        tresholds = {'r1' : [0, 3.8], 'r2' : [3.8, 5.37], 
                     'r3' : [5.37, 6.5], 'r4' : [6.5, 7.6],
                     'r5' : [7.6, 8.5], 'r6' : [8.5, 9.5]}
        
        self.HC_counts = {}
        self.HC_scores = {}
        for k, v in tresholds.items():
            HC_tmp = np.zeros_like(distances)
            mask = (distances > v[0]) & (distances <= v[1])
            HC_tmp[mask] = 1
            self.HC_counts[k] = np.sum(HC_tmp)
            
            HC2 = np.zeros_like(distances)
            HC_allowed_1 = vdw_radii + v[0]
            HC_allowed_2 = vdw_radii + v[1]
            mask1 = distances <= HC_allowed_1
            mask2 = (distances > HC_allowed_1) & (distances <= HC_allowed_2)            
            HC2[mask1] = 1
            HC2[mask2] = (1 / 1.5) * ((vdw_radii[mask2] + 2.0) ** 2 - distances[mask2] ** 2)
            self.HC_scores[k] = np.sum(HC2)

        #hc DISTANCES
        #very short- (r1: 0 < x ≤ 3.8 Å, where x is the measured inter-Cαs distance), 
        #short- (r2: 3.8 < x ≤ 5.37 Å; r3: 5.37 < x ≤ 6.5 Å), 
        #and medium-range (r4: 6.5 < x ≤ 7.6 Å; r5: 7.6 < x ≤ 8.5 Å; r6: 8.5 < x ≤ 9.5 Å) interactions. 

    
    def get_results(self):
        """Return a dictionary of all calculated results."""
        results = {"HC total": self.HC_total, 
            "HC adj": self.HC_adj,
            "RN atoms": self.receptor_data['protein_length'],
            "LN atoms": self.ligand_data['protein_length']}
        for k, v in self.HC_counts.items():
            results['HC_count '+k] = v
        for k, v in self.HC_scores.items():
            results['HC_score '+k] = v
            
        return results
            

class HATAnalyzer():
    def __init__(self, structure, ligand_chains, receptor_chains):
        
        self.ligand_chains = ligand_chains
        self.receptor_chains = receptor_chains
        
        self.structure = structure
        self.load_param = load_param
        self.calculate_distances = calculate_pairwise_distances_intermolecular
        
        self.distance_threshold = 0.7
        
        self.element_types = ['N', 'O', 'S', 'P']
        #permitido = r1+r2;
		#permitido1 = permitido-0.7;  
        #DIS2>=permitido1 && DIS2<=permitido
        self.hydro_param = self.load_param('hydrophobicity.param')
        self.tension_param = self.load_param('tension.param')
        self.atom_radii = self.load_param('pdb_radii.param')
        
        self.allowed_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", 
                                 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        
    def __call__(self):
        #self.extract_atoms()        
        self.ligand_data = self.extract_atoms(self.ligand_chains)
        self.receptor_data = self.extract_atoms(self.receptor_chains)
        #self.distances = self.calculate_distances(self.coords, return_dict=True)
        self.distances = self.calculate_distances(self.receptor_data['coords'], self.ligand_data['coords'], return_dict=True)
        #atom_pairs = np.asarray(list(self.distances.keys()))
        #atom_radii_1 = self.radii[atom_pairs[:,0]]
        #atom_radii_2 = self.radii[atom_pairs[:,1]]
        self.receptor_data, self.ligand_data = self.filter_residues_by_contact_distance(self.receptor_data, self.ligand_data)
        
        self.contact_hydrophobicity_1 = self.sum_param_values(self.receptor_data['filtered_residue_list'], self.hydro_param)
        self.contact_hydrophobicity_2 = self.sum_param_values(self.ligand_data['filtered_residue_list'], self.hydro_param)
        
        self.contact_hydrophobicity = self.contact_hydrophobicity_1 + self.contact_hydrophobicity_2
        
        self.contact_tension_1 = self.sum_param_values(self.receptor_data['filtered_residue_list'], self.tension_param)
        self.contact_tension_2 = self.sum_param_values(self.ligand_data['filtered_residue_list'], self.tension_param)
        
        self.contact_tension = self.contact_tension_1 + self.contact_tension_2
        
        self.total_hydrophobicity_1 = self.sum_param_values(self.receptor_data['res_names'], self.hydro_param)
        self.total_hydrophobicity_2 = self.sum_param_values(self.ligand_data['res_names'], self.hydro_param)
        
        self.total_hydrophobicity = self.total_hydrophobicity_1 + self.total_hydrophobicity_2
        
        self.total_tension_1 = self.sum_param_values(self.receptor_data['res_names'], self.tension_param)
        self.total_tension_2 = self.sum_param_values(self.ligand_data['res_names'], self.tension_param)
        
        self.total_tension = self.total_tension_1 + self.total_tension_2
        
        self.receptor_data = self.calculate_ASA(self.receptor_data, self.receptor_chains)
        self.ligand_data = self.calculate_ASA(self.ligand_data, self.ligand_chains)
        
        self.total_ASA = self.receptor_data['total_ASA'] + self.ligand_data['total_ASA']
        self.polar_ASA = self.receptor_data['polar_ASA'] + self.ligand_data['polar_ASA']
        self.apolar_ASA = self.receptor_data['apolar_ASA'] + self.ligand_data['apolar_ASA']
        self.contact_ASA = self.receptor_data['contact_ASA'] + self.ligand_data['contact_ASA']
        
    def sum_param_values(self, res_name_list, param):
        
        param_sum = 0
        for res_name in res_name_list:
            if res_name in param:
                param_sum += param[res_name]
        
        return param_sum
                
    def extract_atoms(self, target_chains):
        
        # Initialize lists to store donor and acceptor information
        atoms = []
        radii = []
        coords = []
        res = []
        res_names = []
        protein_length = 0

        # Iterate through all atoms in the structure
        for model in self.structure:
            for chain in model:
                if chain.get_id() not in target_chains:
                    continue
                for residue in chain:
                    res_name = residue.get_resname()
                    if res_name not in self.allowed_residues:
                        continue
                    res_id = residue.get_id()[1]  # Get residue number
                    for atom in residue:
                        atom_name = atom.get_name()
                        element_name = atom.element
                        atom_id = atom.get_serial_number()
                        coord = atom.get_coord()

                        if atom_name not in self.atom_radii:
                            continue  # Skip this atom if its radius is not known

                        protein_length += 1
                        # Check if atom is a potential acceptor
                        if (res_name in self.hydro_param or res_name in self.tension_param) and element_name in self.element_types:
                            atoms.append(atom_id)
                            radii.append(self.atom_radii[atom_name])
                            coords.append(coord)
                            res.append(res_id)
                            res_names.append(res_name)

        # Convert lists to numpy arrays
        atoms = np.array(atoms, dtype=int)
        radii = np.array(radii, dtype=float)
        coords = np.array(coords, dtype=float)
        res = np.array(res, dtype=int)
        res_names = np.array(res_names, dtype=str)
        
        data_dict = {'atoms': atoms,
                     'radii': radii,
                     'coords': coords,
                     'res': res,
                     'res_names': res_names,
                     'protein_length': protein_length}
        return data_dict

    def filter_residues_by_contact_distance(self, data_dict_1, data_dict_2):
        # Convert inputs to numpy arrays if they aren't already
        atom_pairs = np.asarray(list(self.distances.keys()))
        atom_radii_1 = data_dict_1['radii'][atom_pairs[:,0]]
        atom_radii_2 = data_dict_2['radii'][atom_pairs[:,1]]
        atom_res_1 = data_dict_1['res'][atom_pairs[:,0]]
        atom_res_2 =  data_dict_2['res'][atom_pairs[:,1]]
        distances = np.asarray(list(self.distances.values()))
        
        # Calculate allowed ranges
        allowed_2 = atom_radii_1 + atom_radii_2
        allowed_1 = allowed_2 - self.distance_threshold
        
        # Create the mask
        mask = (distances >= allowed_1) & (distances <= allowed_2)
        valid_pairs = np.where(mask)[0]
        valid_indices = atom_pairs[valid_pairs]
        
        # Use the valid indices to filter the arrays
        data_dict_1['filtered_atoms'] = data_dict_1['atoms'][valid_indices[:,0]]
        data_dict_2['filtered_atoms'] = data_dict_2['atoms'][valid_indices[:,1]]
        
        data_dict_1['filtered_radii'] = data_dict_1['radii'][valid_indices[:,0]]
        data_dict_2['filtered_radii'] = data_dict_2['radii'][valid_indices[:,1]]
        
        data_dict_1['filtered_res'] = data_dict_1['res'][valid_indices[:,0]]
        data_dict_2['filtered_res'] = data_dict_2['res'][valid_indices[:,1]]
        
        unique_elements, indices = np.unique(data_dict_1['filtered_res'], return_index=True)
        data_dict_1['filtered_residue_list'] = data_dict_1['res_names'][indices]    
        
        unique_elements, indices = np.unique(data_dict_2['filtered_res'], return_index=True)
        data_dict_2['filtered_residue_list'] = data_dict_2['res_names'][indices]    
 
        return data_dict_1, data_dict_2
    

    def calculate_ASA(self, data_dict, target_chains):
        
        from copy import deepcopy
        from tempfile import NamedTemporaryFile
              
        def remove_hydrogens(structure, target_chains, allowed_residues):
            # Create a working copy to avoid modifying original structure
            io = PDBIO()
            io.set_structure(structure)

            class CustomSelect(Select):
                def accept_chain(self, chain):
                    return chain.id in target_chains

                def accept_residue(self, residue):
                    return residue.get_resname() in allowed_residues

                def accept_atom(self, atom):
                    return atom.element != 'H' and atom.get_serial_number() in data_dict['filtered_atoms']

            # Use temporary file to force full structure rebuild
            with NamedTemporaryFile(delete=True) as tmp:
                io.save(tmp.name, CustomSelect(), preserve_atom_numbering=True)
                parser = PDBParser(QUIET=True)
                return parser.get_structure("clean", tmp.name)       

        def keep_allowed_res(structure, allowed_residues):

            io = PDBIO()
            io.set_structure(structure)

            class CustomSelect(Select):
                def accept_residue(self, residue):
                    return residue.get_resname() in allowed_residues

            with NamedTemporaryFile(delete=True) as tmp:
                io.save(tmp.name, CustomSelect(), preserve_atom_numbering=True)
                parser = PDBParser(QUIET=True)
                return parser.get_structure("clean", tmp.name)       

        safe_structure = deepcopy(self.structure)
        clean_structure = keep_allowed_res(safe_structure, self.allowed_residues)

        asa_result, asa_classes = freesasa.calcBioPDB(clean_structure)
        total_ASA = asa_result.totalArea()
    
        data_dict['total_ASA'] = total_ASA
        if "Polar" in asa_classes.keys():
            data_dict['polar_ASA'] = asa_classes['Polar']
        else:
            data_dict['polar_ASA'] = 0
        if "Apolar" in asa_classes.keys():     
            data_dict['apolar_ASA'] = asa_classes['Apolar']
        else:
            data_dict['apolar_ASA'] = 0

        new_structure = remove_hydrogens(safe_structure, target_chains, self.allowed_residues)

        is_empty = True
        for model in new_structure:
            if is_empty == False:
                break
            for chain in model:
                if is_empty == False:
                    break
                for residue in chain:
                    if is_empty == False:
                        break
                    for atom in residue:
                        is_empty = False
                        break

        if is_empty == False:
            asa_result, asa_classes = freesasa.calcBioPDB(new_structure)
            total_ASA = asa_result.totalArea()
        else:
            total_ASA = 0
            asa_classes = {}

        if "Polar" in asa_classes.keys():
            data_dict['contact_polar_ASA'] = asa_classes['Polar']
        else:
            data_dict['contact_polar_ASA'] = 0
        if "Apolar" in asa_classes.keys():
            data_dict['contact_apolar_ASA'] = asa_classes['Apolar']
        else:
            data_dict['contact_apolar_ASA'] = 0

        data_dict['contact_ASA'] = total_ASA

        if len(data_dict['filtered_atoms']) > 0:
            data_dict['mean_contact_ASA'] = total_ASA/len(data_dict['filtered_atoms'])
        else:
            data_dict['mean_contact_ASA'] = total_ASA
        '''
        i = 0
        for model in new_structure:
            for chain in model:
                #print('chain', chain.id)
                for residue in chain:
                    for atom in residue:
                        atom_number = atom.get_serial_number()
                        atom_id = atom.get_full_id()
                        #print(i, residue.get_resname(), target_chains, atom_id)
                        asa = asa_result.atomArea(i)
                        i += 1 
                        if atom_number in data_dict['filtered_atoms']:
                            contact_ASA += asa
                            print('a',atom_number, 'asa', asa)
        
        #print("chain", target_chains, "contact asa", contact_ASA)
        data_dict['total_ASA'] = total_ASA
        data_dict['polar_ASA'] = polar_ASA
        data_dict['apolar_ASA'] = apolar_ASA
        data_dict['contact_ASA'] = contact_ASA
        '''
        return data_dict
                            
    def get_results(self):
        """Return a dictionary of all calculated results."""
        return {
            "Contact hydrophobicity": self.contact_hydrophobicity, 
            "Total hydrophobicity": self.total_hydrophobicity, 
            "Contact Tension": self.contact_tension,
            "Total Tension": self.total_tension,
            "Contact ASA": self.contact_ASA,
            "Total ASA": self.total_ASA,
            "Polar ASA": self.polar_ASA,
            "Apolar ASA": self.apolar_ASA,
            "RContact hydrophobicity": self.contact_hydrophobicity_1, 
            "RTotal hydrophobicity": self.total_hydrophobicity_1, 
            "RContact Tension": self.contact_tension_1,
            "RTotal Tension": self.total_tension_1,
            "RContact ASA": self.receptor_data['contact_ASA'],
            "RMeanContact ASA": self.receptor_data['mean_contact_ASA'],
            "RContactPolar ASA": self.receptor_data['contact_polar_ASA'],
            "RContactApolar ASA": self.receptor_data['contact_apolar_ASA'],
            "RTotal ASA": self.receptor_data['total_ASA'],
            "RPolar ASA": self.receptor_data['polar_ASA'],
            "RApolar ASA": self.receptor_data['apolar_ASA'],
            "LContact hydrophobicity": self.contact_hydrophobicity_2, 
            "LTotal hydrophobicity": self.total_hydrophobicity_2, 
            "LContact Tension": self.contact_tension_2,
            "LTotal Tension": self.total_tension_2,
            "LContact ASA": self.ligand_data['contact_ASA'],
            "LMeanContact ASA": self.ligand_data['mean_contact_ASA'],
            "LContactPolar ASA": self.ligand_data['contact_polar_ASA'],
            "LContactApolar ASA": self.ligand_data['contact_apolar_ASA'],
            "LTotal ASA": self.ligand_data['total_ASA'],
            "LPolar ASA": self.ligand_data['polar_ASA'],
            "LApolar ASA": self.ligand_data['apolar_ASA'],
        }                            

class HBAnalyzer:
    def __init__(self, structure, ligand_chains, receptor_chains, distance_threshold=0.7):

        self.structure = structure
        self.ligand_chains = ligand_chains
        self.receptor_chains = receptor_chains
        
        self.load_param = load_param
        self.calculate_distances = calculate_pairwise_distances_intermolecular
        
        self.atom_radii = self.load_param('pdb_radii.param')
        
        self.distance_threshold = distance_threshold
        self.angle_threshold_acceptor = 60
        self.angle_threshold_donor = 60
                
        self.hb_donors = {
            'ALA': ['N'],
            'ARG': ['N', 'NE', 'NH1', 'NH2'],
            'ASN': ['N', 'ND2'],
            'ASP': ['N'],
            'CYS': ['N', 'SG'],  # SG can be a weak donor when protonated
            'GLN': ['N', 'NE2'],
            'GLU': ['N'],
            'GLY': ['N'],
            'HIS': ['N', 'ND1', 'NE2'],  # ND1 or NE2 can be donors depending on protonation state
            'ILE': ['N'],
            'LEU': ['N'],
            'LYS': ['N', 'NZ'],
            'MET': ['N'],
            'PHE': ['N'],
            'PRO': [],  # Proline N is not typically a hydrogen bond donor
            'SER': ['N', 'OG'],
            'THR': ['N', 'OG1'],
            'TRP': ['N', 'NE1'],
            'TYR': ['N', 'OH'],
            'VAL': ['N'],
            'SEP': ['N', 'OG'],
            # Phosphothreonine (TPO)
            'TPO': ['N', 'OG1'],
            # Phosphotyrosine (PTR)
            'PTR': ['N', 'OH']

        }

        self.hb_acceptors = {
            'ALA': ['O'],
            'ARG': ['O'],
            'ASN': ['O', 'OD1'],
            'ASP': ['O', 'OD1', 'OD2'],
            'CYS': ['O'],
            'GLN': ['O', 'OE1'],
            'GLU': ['O', 'OE1', 'OE2'],
            'GLY': ['O'],
            'HIS': ['O', 'ND1', 'NE2'],
            'ILE': ['O'],
            'LEU': ['O'],
            'LYS': ['O'],
            'MET': ['O'],
            'PHE': ['O'],
            'PRO': ['O'],
            'SER': ['O', 'OG'],
            'THR': ['O', 'OG1'],
            'TRP': ['O'],
            'TYR': ['O', 'OH'],
            'VAL': ['O'],
            'SEP': ['O', 'OG', 'O1P', 'O2P', 'O3P'],
            # Phosphothreonine (TPO)
            'TPO': ['O', 'O1P', 'O2P', 'O3P'],
            # Phosphotyrosine (PTR)
            'PTR': ['O', 'OH', 'O1P', 'O2P', 'O3P']

        }
        
        self.max_hbonds = { #Maximum number of bonds that the donors and acceptors can take part in
            'N': 2,  # Amide nitrogen can typically donate two H-bonds
            'NE': 1,
            'NH1': 2,
            'NH2': 2,
            'ND2': 1,
            'SG': 1,
            'NE2': 1,
            'ND1': 1,
            'NZ': 3,
            'OG': 2,
            'OG1': 2,
            'NE1': 1,
            'OH': 1,
            'O': 2,   # Carbonyl oxygen can typically accept two H-bonds
            'OD1': 2,
            'OD2': 2,
            'OE1': 2,
            'OE2': 2
        }
 
        #ref: https://www.rcsb.org/ligand/ARG (replace ARG by the name - 3 letter symbol - of the residue you are interested)
        self.covalent_bonds = {
            'ALA': {'N': ['CA'], 'O': ['C']},
            'ARG': {'N': ['CA'], 'NE': ['CD', 'CZ'], 'NH1': ['CZ'], 'NH2': ['CZ'], 'O': ['C'], 'OXT': ['C']},
            'ASN': {'N': ['CA'], 'ND2': ['CG'], 'O': ['C'], 'OD1': ['CG'], 'OXT': ['C']},
            'ASP': {'N': ['CA'], 'O': ['C'], 'OD1': ['CG'], 'OD2': ['CG'], 'OXT': ['C']},
            'CYS': {'N': ['CA'], 'SG': ['CB'], 'O': ['C'], 'OXT': ['C']},
            'GLN': {'N': ['CA'], 'NE2': ['CD'], 'O': ['C'], 'OE1': ['CD'], 'OXT': ['C']},
            'GLU': {'N': ['CA'], 'O': ['C'], 'OE1': ['CD'], 'OE2': ['CD'], 'OXT': ['C']},
            'GLY': {'N': ['CA'], 'O': ['C'], 'OXT': ['C']},
            'HIS': {'N': ['CA'], 'ND1': ['CG', 'CE1'], 'NE2': ['CD2', 'CE1'], 'O': ['C'], 'OXT': ['C']},
            'ILE': {'N': ['CA'], 'O': ['C'], 'OXT': ['C']},
            'LEU': {'N': ['CA'], 'O': ['C'], 'OXT': ['C']},
            'LYS': {'N': ['CA'], 'NZ': ['CE'], 'O': ['C'], 'OXT': ['C']},
            'MET': {'N': ['CA'], 'O': ['C'], 'OXT': ['C']},
            'PHE': {'N': ['CA'], 'O': ['C'], 'OXT': ['C']},
            'PRO': {'O': ['C'], 'OXT': ['C']},  # Proline N is not typically a hydrogen bond donor (only is when protonated), Interactor was wrong here
            'SER': {'N': ['CA'], 'OG': ['CB'], 'O': ['C'], 'OXT': ['C']},
            'THR': {'N': ['CA'], 'OG1': ['CB'], 'O': ['C'], 'OXT': ['C']},
            'TRP': {'N': ['CA'], 'NE1': ['CD1', 'CE2'], 'O': ['C'], 'OXT': ['C']},
            'TYR': {'N': ['CA'], 'OH': ['CZ'], 'O': ['C'], 'OXT': ['C']},
            'VAL': {'N': ['CA'], 'O': ['C'], 'OXT': ['C']},
            # Phosphoserine (SEP)
            'SEP': {'N': ['CA'], 'O': ['C'],'OXT': ['C'],'OG': ['CB', 'P'], 'P': ['OG', 'O1P', 'O2P', 'O3P']},
            # Phosphothreonine (TPO)
            'TPO': {'N': ['CA'],'O': ['C'],'OXT': ['C'],'OG1': ['CB', 'P'],'P': ['OG1', 'O1P', 'O2P', 'O3P']},
            # Phosphotyrosine (PTR)
            'PTR': {'N': ['CA'],'O': ['C'],'OXT': ['C'],'OH': ['CZ', 'P'], 'P': ['OH', 'O1P', 'O2P', 'O3P']}
        }

        self.rotable_bonds = {
            'ALA': {'CA': 3, 'C': 0},
            'ARG': {'CA': 3, 'C': 0, 'CD': 2, 'CZ': 0}, 
            'ASN': {'CA': 3, 'C': 0, 'CG': 1},
            'ASP': {'CA': 3, 'C': 0, 'CG': 0},
            'CYS': {'CA': 3, 'C': 0, 'CB': 2},
            'GLN': {'CA': 2, 'C': 0, 'CD': 2},
            'GLU': {'CA': 2, 'C': 0, 'CD': 2},
            'GLY': {'CA': 4, 'C': 0},  # Glycine is more flexible
            'HIS': {'CA': 3, 'C': 0, 'CD2': 0, 'CE1': 0, 'CG': 0},
            'ILE': {'CA': 3, 'C': 0},
            'LEU': {'CA': 3, 'C': 0},
            'LYS': {'CA': 3, 'C': 0, 'CE': 2},
            'MET': {'CA': 3, 'C': 0},
            'PHE': {'CA': 3, 'C': 0},
            'PRO': {'C': 0},
            'SER': {'CA': 3, 'C': 0, 'CB': 2},
            'SEP': {'CA': 3, 'C': 0, 'CB': 2, 'P': 0},
            'THR': {'CA': 3, 'C': 0, 'CB': 3},
            'TPO': {'CA': 3, 'C': 0, 'CB': 3, 'P':0},
            'TRP': {'CA': 3, 'C': 0},
            'TYR': {'CA': 3, 'C': 0, 'CZ': 0},
            'PTR': {'CA': 3, 'C': 0, 'CZ': 0, 'P': 0},
            'VAL': {'CA': 3, 'C': 0}
        }
        #CONTINUE HERE
        self.empty_results = {
            "Total HB": 0,
            "MC-MC HB": 0,
            "MC-SD HB": 0,
            "SD-MC HB": 0,
            "SD-SD HB": 0,
            "Total RT": 0,
            "LMC-MC HB": 0,
            "LMC-SD HB": 0,
            "LSD-MC HB": 0,
            "LSD-SD HB": 0,
            "RMC-MC HB": 0,
            "RMC-SD HB": 0,
            "RSD-MC HB": 0,
            "RSD-SD HB": 0,
            "Ligand RT": 0,
            "Receptor RT": 0,
            "LDonor RT": 0,
            "RDonor RT": 0,
            "LAcceptor RT": 0,
            "RAcceptor RT": 0
        }
        #CONTINUE HERE
        self.empty_hb_list = {
                'Donor_Chain': 'NA',
                'Donor_Residue_Name': 'NA',
                'Donor_Residue_Num': 'NA',
                'Acceptor_Chain': 'NA',
                'Acceptor_Residue_Name': 'NA',
                'Acceptor_Residue_Num': 'NA',
                'Donor_Atom': 'NA',
                'Acceptor_Atom': 'NA',
                'Distance': 'NA',
                'Donor_Angle': 'NA',
                'Acceptor_Angle': 'NA',
                'Type': 'NA'
            }

    def __call__(self):
        
        self.ligand_data = self.extract_hb_atoms(self.ligand_chains)
        self.receptor_data = self.extract_hb_atoms(self.receptor_chains)

        distances_1 = self.calculate_distances(self.receptor_data['donor_coords'], self.ligand_data['acceptor_coords'], return_dict=True)
        distances_2 = self.calculate_distances(self.ligand_data['donor_coords'], self.receptor_data['acceptor_coords'], return_dict=True)
        print("0000000000000000", len(distances_1.values()), np.min(list(distances_1.values())))
        
        
        distances_1_f, self.receptor_data, self.ligand_data = self.filter_distances(distances_1, self.receptor_data, self.ligand_data)
        distances_2_f, self.ligand_data, self.receptor_data = self.filter_distances(distances_2, self.ligand_data, self.receptor_data)
        print("AAAAAAAAAAAAAAAA", len(distances_1_f.values()), np.min(list(distances_1_f.values())))

        self.total_rotable_bond_count_1, self.receptor_data, self.ligand_data = self.find_acceptor_donor_roots(self.receptor_data, self.ligand_data)
        self.total_rotable_bond_count_2, self.ligand_data, self.receptor_data = self.find_acceptor_donor_roots(self.ligand_data, self.receptor_data)
        
        self.total_rotable_bond_count = self.total_rotable_bond_count_1 + self.total_rotable_bond_count_2
        
        self.calculate_angles()

        if 'donor_root_angles' in self.receptor_data.keys():
            distances_1_fa, self.receptor_data, self.ligand_data = self.filter_angles(distances_1_f, self.receptor_data, self.ligand_data)
        if 'donor_root_angles' in  self.ligand_data.keys():
            distances_2_fa, self.ligand_data, self.receptor_data = self.filter_angles(distances_2_f, self.ligand_data, self.receptor_data)        
        
        self.total_hb_1, self.category_counts_1, self.hb_df1 = self.count_hydrogen_bonds(self.receptor_data, self.ligand_data, distances_1_fa)
        self.total_hb_2, self.category_counts_2, self.hb_df2  = self.count_hydrogen_bonds(self.ligand_data, self.receptor_data, distances_2_fa)
        
        print("BBBBBBBBBBBBBB",len(self.hb_df1), len(distances_1_fa.values()), np.min(list(distances_1_fa.values())))

        self.total_hb = self.total_hb_1 + self.total_hb_2
        
        self.category_counts = {k1:v1+self.category_counts_2[k1] for k1,v1 in self.category_counts_1.items()}
        

    def extract_hb_atoms(self, target_chains):
        # Initialize lists to store donor and acceptor information
        donor_atoms = []
        donor_atom_names = []
        donor_radii = []
        donor_coords = []
        donor_res = []
        donor_chains = []
        donor_resnames = []
        acceptor_atoms = []
        acceptor_atom_names = []
        acceptor_radii = []
        acceptor_coords = []
        acceptor_res = []
        acceptor_chains = []
        acceptor_resnames = []

        # Iterate through all atoms in the structure
        for model in self.structure:
            for chain in model:
                if chain.get_id() not in target_chains:
                    continue
                for residue in chain:
                    res_name = residue.get_resname()
                    res_id = residue.get_id()[1]  # Get residue number
                    for atom in residue:
                        atom_name = atom.get_name()
                        atom_id = atom.get_serial_number()
                        coord = atom.get_coord()

                        if atom_name not in self.atom_radii:
                            continue  # Skip this atom if its radius is not known

                        # Check if atom is a potential donor
                        if res_name in self.hb_donors and atom_name in self.hb_donors[res_name]:
                            donor_atoms.append(atom_id)
                            donor_atom_names.append(atom_name)
                            donor_radii.append(self.atom_radii[atom_name])
                            donor_coords.append(coord)
                            donor_res.append(res_id)
                            donor_chains.append(chain.get_id())
                            donor_resnames.append(res_name)

                        # Check if atom is a potential acceptor
                        if res_name in self.hb_acceptors and atom_name in self.hb_acceptors[res_name]:
                            acceptor_atoms.append(atom_id)
                            acceptor_atom_names.append(atom_name)
                            acceptor_radii.append(self.atom_radii[atom_name])
                            acceptor_coords.append(coord)
                            acceptor_res.append(res_id)
                            acceptor_chains.append(chain.get_id())
                            acceptor_resnames.append(res_name)

        # Convert lists to numpy arrays
        donor_atoms = np.array(donor_atoms, dtype=int)
        donor_atom_names = np.array(donor_atom_names, dtype=str)
        donor_radii = np.array(donor_radii, dtype=float)
        donor_coords = np.array(donor_coords, dtype=float)
        donor_res = np.array(donor_res, dtype=int)
        donor_chains = np.array(donor_chains, dtype=str)
        donor_resnames = np.array(donor_resnames, dtype=str)
        acceptor_atoms = np.array(acceptor_atoms, dtype=int)
        acceptor_atom_names = np.array(acceptor_atom_names, dtype=str)
        acceptor_radii = np.array(acceptor_radii, dtype=float)
        acceptor_coords = np.array(acceptor_coords, dtype=float)
        acceptor_res = np.array(acceptor_res, dtype=int) 
        acceptor_chains = np.array(acceptor_chains, dtype=str)
        acceptor_resnames = np.array(acceptor_resnames, dtype=str)
        
        result_dict = {'donor_atoms': donor_atoms, 
                       'donor_atom_names': donor_atom_names, 
                       'donor_radii' : donor_radii, 
                       'donor_coords' : donor_coords,
                       'donor_res' : donor_res,
                       'donor_chains': donor_chains,
                       'donor_resnames': donor_resnames,
                       'acceptor_atoms' : acceptor_atoms,
                       'acceptor_atom_names' : acceptor_atom_names,
                       'acceptor_radii': acceptor_radii, 
                       'acceptor_coords': acceptor_coords,
                       'acceptor_res': acceptor_res,
                       'acceptor_chains': acceptor_chains,
                       'acceptor_resnames': acceptor_resnames}
        
        return result_dict
    '''
    def filter_distances(self, distances_1, data_dict_1, data_dict_2):
        
        # Convert inputs to numpy arrays if they aren't already
        atom_pairs = np.asarray(list(distances_1.keys()))
        atom_radii_1 = data_dict_1['donor_radii'][atom_pairs[:,0]]
        atom_radii_2 = data_dict_2['acceptor_radii'][atom_pairs[:,1]]
        donor_res = data_dict_1['donor_res'][atom_pairs[:,0]]
        acceptor_res = data_dict_2['acceptor_res'][atom_pairs[:,1]]
        distances = np.asarray(list(distances_1.values()))
            
        allowed_distances = atom_radii_1 + atom_radii_2 + self.distance_threshold
               
        dist_mask = (distances <= allowed_distances)
        valid_mask = dist_mask & (donor_res != acceptor_res)
        
        # Find indices of valid pairs
        valid_pairs = np.where(valid_mask)[0]
        valid_indices = atom_pairs[valid_pairs]

        # Update attributes
        data_dict_1['donor_atoms'] = data_dict_1['donor_atoms'][valid_indices[:,0]]
        data_dict_1['donor_atom_names'] = data_dict_1['donor_atom_names'][valid_indices[:,0]]
        data_dict_1['donor_radii'] = data_dict_1['donor_radii'][valid_indices[:,0]]
        data_dict_1['donor_coords'] = data_dict_1['donor_coords'][valid_indices[:,0]]
        data_dict_1['donor_res'] = data_dict_1['donor_res'][valid_indices[:,0]]

        data_dict_2['acceptor_atoms'] = data_dict_2['acceptor_atoms'][valid_indices[:,1]]
        data_dict_2['acceptor_atom_names'] = data_dict_2['acceptor_atom_names'][valid_indices[:,1]]
        data_dict_2['acceptor_radii'] = data_dict_2['acceptor_radii'][valid_indices[:,1]]
        data_dict_2['acceptor_coords'] = data_dict_2['acceptor_coords'][valid_indices[:,1]]
        data_dict_2['acceptor_res'] = data_dict_2['acceptor_res'][valid_indices[:,1]]

        #print("AFTER FILTERING DISTANCES", len(self.donor_atoms), len(self.donor_res), len(self.donor_atom_names), len(self.acceptor_atoms), len(self.acceptor_res), len(self.acceptor_atom_names))
        # Convert valid_indices to a set of tuples for faster lookup
        
        # Update distances
        distances_1 = {k: v for k, v in distances_1.items() if k in valid_indices}
        
        return distances_1, data_dict_1, data_dict_2
    '''
    def filter_distances(self, distances_1, data_dict_1, data_dict_2):
        # Convert dict keys/values to aligned NumPy arrays
        atom_pairs = np.asarray(list(distances_1.keys()))       # shape (n_pairs, 2) → (donor_idx, acceptor_idx)
        distances_arr = np.asarray(list(distances_1.values()))  # shape (n_pairs,)

        # Radii for both atoms in each pair
        atom_radii_1 = data_dict_1['donor_radii'][atom_pairs[:, 0]]
        atom_radii_2 = data_dict_2['acceptor_radii'][atom_pairs[:, 1]]

        # Residue numbers (to avoid same-residue bonds)
        donor_res    = data_dict_1['donor_res'][atom_pairs[:, 0]]
        acceptor_res = data_dict_2['acceptor_res'][atom_pairs[:, 1]]

        # Allowed max distance = sum of radii + extra threshold
        allowed_distances = atom_radii_1 + atom_radii_2 + self.distance_threshold

        # Mask: distances within allowed range
        dist_mask = distances_arr <= allowed_distances

        # Mask: donors and acceptors not from same residue
        valid_mask = dist_mask & (donor_res != acceptor_res)

        # Get index positions that are valid
        valid_indices = np.where(valid_mask)[0]

        # Slice arrays to only keep valid pairs
        atom_pairs    = atom_pairs[valid_indices]
        distances_arr = distances_arr[valid_indices]

        # Filter donor side arrays
        #data_dict_1['donor_atoms']      = data_dict_1['donor_atoms'][atom_pairs[:, 0]]
        #data_dict_1['donor_atom_names'] = data_dict_1['donor_atom_names'][atom_pairs[:, 0]]
        #data_dict_1['donor_radii']      = data_dict_1['donor_radii'][atom_pairs[:, 0]]
        #data_dict_1['donor_coords']     = data_dict_1['donor_coords'][atom_pairs[:, 0]]
        #data_dict_1['donor_res']        = data_dict_1['donor_res'][atom_pairs[:, 0]]
        #if 'donor_chains' in data_dict_1:     # in case chain IDs were added
        #    data_dict_1['donor_chains'] = data_dict_1['donor_chains'][atom_pairs[:, 0]]
        #if 'donor_resnames' in data_dict_1:   # in case residue names were added
        #    data_dict_1['donor_resnames'] = data_dict_1['donor_resnames'][atom_pairs[:, 0]]

        # Filter acceptor side arrays
        #data_dict_2['acceptor_atoms']      = data_dict_2['acceptor_atoms'][atom_pairs[:, 1]]
        #data_dict_2['acceptor_atom_names'] = data_dict_2['acceptor_atom_names'][atom_pairs[:, 1]]
        #data_dict_2['acceptor_radii']      = data_dict_2['acceptor_radii'][atom_pairs[:, 1]]
        #data_dict_2['acceptor_coords']     = data_dict_2['acceptor_coords'][atom_pairs[:, 1]]
        #data_dict_2['acceptor_res']        = data_dict_2['acceptor_res'][atom_pairs[:, 1]]
        #if 'acceptor_chains' in data_dict_2:
        #    data_dict_2['acceptor_chains'] = data_dict_2['acceptor_chains'][atom_pairs[:, 1]]
        #if 'acceptor_resnames' in data_dict_2:
        #    data_dict_2['acceptor_resnames'] = data_dict_2['acceptor_resnames'][atom_pairs[:, 1]]

        # Rebuild distances dict but now in the same order as the filtered donor/acceptor arrays
        distances_1 = dict(zip(map(tuple, atom_pairs), distances_arr))

        # Get the keys in final output order
        sorted_keys = list(distances_1.keys())
        donor_idxs = np.array([k[0] for k in sorted_keys])
        acceptor_idxs = np.array([k[1] for k in sorted_keys])

        for key in ['donor_atoms', 'donor_atom_names', 'donor_radii', 'donor_coords', 'donor_res', 'donor_chains', 'donor_resnames']:
            if key in data_dict_1:
                data_dict_1[key] = data_dict_1[key][donor_idxs]
        for key in ['acceptor_atoms', 'acceptor_atom_names', 'acceptor_radii', 'acceptor_coords', 'acceptor_res', 'acceptor_chains', 'acceptor_resnames']:
            if key in data_dict_2:
                data_dict_2[key] = data_dict_2[key][acceptor_idxs]
        # After reordering donor and acceptor arrays, reset indices to match new order:
        new_keys = list(zip(range(len(donor_idxs)), range(len(acceptor_idxs))))
        distances_1 = dict(zip(new_keys, distances_arr))

        return distances_1, data_dict_1, data_dict_2


    def find_acceptor_donor_roots(self, data_dict_1, data_dict_2):
        donor_roots = []
        acceptor_roots = []
        
        donor_root_rotable_bond_count = 0
        acceptor_root_rotable_bond_count = 0
        
        donors_included = []
        acceptors_included = []
        #print("BEFORE FINDING ROOTS", len(self.donor_atoms), len(self.donor_res), len(self.donor_atom_names), len(self.acceptor_atoms), len(self.acceptor_res), len(self.acceptor_atom_names))
        for donor_id, donor_res, donor_name, acceptor_id, acceptor_res, acceptor_name in zip(data_dict_1['donor_atoms'], data_dict_1['donor_res'], data_dict_1['donor_atom_names'],
                                                                  data_dict_2['acceptor_atoms'], data_dict_2['acceptor_res'], data_dict_2['acceptor_atom_names']):
            
            donor_root, donor_root_rotable_bonds = self.get_covalent_root(donor_id, donor_res, donor_name)
            donor_roots.append(donor_root)
            acceptor_root, acceptor_root_rotable_bonds = self.get_covalent_root(acceptor_id, acceptor_res, acceptor_name)
            acceptor_roots.append(acceptor_root)

            if donor_id not in donors_included:
                donor_root_rotable_bond_count += donor_root_rotable_bonds
                donors_included.append(donor_id)
            
            if acceptor_id not in acceptors_included:
                acceptor_root_rotable_bond_count += acceptor_root_rotable_bonds
                acceptors_included.append(acceptor_id)
                
            
        total_rotable_bond_count = donor_root_rotable_bond_count + acceptor_root_rotable_bond_count
        
        data_dict_1['donor_roots'] = donor_roots
        data_dict_2['acceptor_roots'] = acceptor_roots
        data_dict_1['donor_root_rotable_bond_count'] = donor_root_rotable_bond_count
        data_dict_2['acceptor_root_rotable_bond_count'] = acceptor_root_rotable_bond_count
        
        return total_rotable_bond_count, data_dict_1, data_dict_2


    def get_covalent_root(self, atom_id, res_id, target_atom_name):
        
        rotable_bond_count = 0
        roots = []
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[1] == res_id or (residue.get_id()[1] == (res_id+1) and target_atom_name == 'N'):
                        res_name = residue.get_resname()
                        if res_name in self.covalent_bonds:
                            if target_atom_name in self.covalent_bonds[res_name]:
                                for potential_root in residue:
                                    potential_root_name = potential_root.get_name()
                                    if residue.get_id()[1] == (res_id+1) and potential_root_name == 'C' and target_atom_name == 'N':
                                        roots.append(potential_root.get_coord())
                                    if residue.get_id()[1] == res_id and potential_root_name in self.covalent_bonds[res_name][target_atom_name]:
                                        rotable_bond_count += self.rotable_bonds[res_name][potential_root_name]
                                        roots.append(potential_root.get_coord())
                                if roots:
                                    return roots, rotable_bond_count
        print("ERROR! ROOT NOT FOUND", atom_id, res_id)
        return [None]  # If atom not found (shouldn't happen if data is consistent)ouldn't happen if data is consistent)

    
    def calculate_angles(self):

        def process_roots(roots):
            processed_roots = []
            for r in roots:
                r = np.asarray(r)
                if r.ndim == 1:
                    processed_roots.append(r)
                else:
                    processed_roots.append(np.mean(r, axis=0))
            return np.array(processed_roots)

        def calculate_angle(A, B, C):
            '''
            v1 = coords1 - coords2
            v2 = root_coords - coords2

            v1_norm = v1 / np.linalg.norm(v1, axis=1)[:, np.newaxis]
            v2_norm = v2 / np.linalg.norm(v2, axis=1)[:, np.newaxis]

            dot_product = np.sum(v1_norm * v2_norm, axis=1)
            dot_product = np.clip(dot_product, -1.0, 1.0)

            return np.arccos(dot_product) * 180 / np.pi
            '''
            # Calculate vectors
            BA = B - A
            BC = B - C

            # Calculate dot products
            dot_products = np.sum(BA * BC, axis=1)

            # Calculate magnitudes
            magnitudes_BA = np.linalg.norm(BA, axis=1)
            magnitudes_BC = np.linalg.norm(BC, axis=1)

            # Calculate cosine of the angles
            cos_angles = dot_products / (magnitudes_BA * magnitudes_BC)

            # Ensure the values are within [-1, 1] to avoid domain errors
            cos_angles = np.clip(cos_angles, -1.0, 1.0)

            # Calculate angles in radians and convert to degrees
            angles = np.arccos(cos_angles)
            angles_degrees = np.degrees(angles)

            return angles_degrees

        donor_roots_processed = process_roots(self.receptor_data['donor_roots'])
        acceptor_roots_processed = process_roots(self.ligand_data['acceptor_roots'])

        if len(self.ligand_data['acceptor_coords']) > 0 and len(self.receptor_data['donor_coords']) > 0:
            self.receptor_data['donor_root_angles'] = calculate_angle(self.ligand_data['acceptor_coords'], self.receptor_data['donor_coords'], donor_roots_processed)
            self.ligand_data['acceptor_root_angles'] = calculate_angle(self.receptor_data['donor_coords'], self.ligand_data['acceptor_coords'], acceptor_roots_processed)

        donor_roots_processed = process_roots(self.ligand_data['donor_roots'])
        acceptor_roots_processed = process_roots(self.receptor_data['acceptor_roots'])

        if len(self.receptor_data['acceptor_coords']) > 0 and len(self.ligand_data['donor_coords']) > 0:
            self.ligand_data['donor_root_angles'] = calculate_angle(self.receptor_data['acceptor_coords'], self.ligand_data['donor_coords'], donor_roots_processed)
            self.receptor_data['acceptor_root_angles'] = calculate_angle(self.ligand_data['donor_coords'], self.receptor_data['acceptor_coords'], acceptor_roots_processed)
    '''
    def filter_angles(self, distances, data_dict_1, data_dict_2):
        #Filters hydrogen bond data based on angle thresholds.
        atom_pairs = np.asarray(list(distances.keys()))

        # Create boolean masks based on the angle thresholds
        acceptor_mask = data_dict_2['acceptor_root_angles'] >= self.angle_threshold_acceptor
        donor_mask = data_dict_1['donor_root_angles'] >= self.angle_threshold_donor

        # Combine the masks:  We want to keep entries that satisfy BOTH angle criteria
        combined_mask = acceptor_mask & donor_mask
        valid_indices = np.where(combined_mask)[0]
        
        # Apply the combined mask to filter the data
        data_dict_1['donor_atoms'] = data_dict_1['donor_atoms'][valid_indices]
        data_dict_1['donor_radii'] = data_dict_1['donor_radii'][valid_indices]
        data_dict_1['donor_coords'] = data_dict_1['donor_coords'][valid_indices]
        data_dict_1['donor_res'] = data_dict_1['donor_res'][valid_indices]

        data_dict_2['acceptor_atoms'] = data_dict_2['acceptor_atoms'][valid_indices]
        data_dict_2['acceptor_radii'] = data_dict_2['acceptor_radii'][valid_indices]
        data_dict_2['acceptor_coords'] = data_dict_2['acceptor_coords'][valid_indices]
        data_dict_2['acceptor_res'] = data_dict_2['acceptor_res'][valid_indices]

        data_dict_2['acceptor_root_angles'] = data_dict_2['acceptor_root_angles'][valid_indices]
        data_dict_1['donor_root_angles'] =  data_dict_1['donor_root_angles'][valid_indices]
        
        # Convert valid_indices to a set of tuples for faster lookup
        valid_indices = atom_pairs[valid_indices]
        valid_pairs_set = set(map(tuple, valid_indices))
        
        # Update distances
        distances = {k: v for k, v in distances.items() if k in valid_pairs_set}
                    
        return distances, data_dict_1, data_dict_2

       '''

    def filter_angles(self, distances, data_dict_1, data_dict_2):
        atom_pairs = np.asarray(list(distances.keys()))
        distance_vals = np.array(list(distances.values()))  # <- capture values in parallel

        acceptor_mask = data_dict_2['acceptor_root_angles'] >= self.angle_threshold_acceptor
        donor_mask    = data_dict_1['donor_root_angles'] >= self.angle_threshold_donor
        combined_mask = acceptor_mask & donor_mask

        valid_indices = np.where(combined_mask)[0]

        #for key in ['donor_atoms','donor_radii','donor_coords','donor_res','donor_root_angles', 'donor_chains']:
        #    data_dict_1[key] = data_dict_1[key][valid_indices]
        #for key in ['acceptor_atoms','acceptor_radii','acceptor_coords','acceptor_res','acceptor_root_angles', 'acceptor_chains']:
        #    data_dict_2[key] = data_dict_2[key][valid_indices]

        # Also filter distances *by index* to preserve alignment
        atom_pairs = atom_pairs[valid_indices]
        distance_vals = distance_vals[valid_indices]
    
        # Rebuild dict with only kept pairs, keeping order matched to arrays
        distances_1 = dict(zip(map(tuple, atom_pairs), distance_vals))

        # Get the keys in final output order
        sorted_keys = list(distances_1.keys())
        donor_idxs = np.array([k[0] for k in sorted_keys])
        acceptor_idxs = np.array([k[1] for k in sorted_keys])

        for key in ['donor_atoms', 'donor_atom_names', 'donor_radii', 'donor_coords', 'donor_res', 'donor_root_angles','donor_chains', 'donor_resnames']:
            if key in data_dict_1:
                data_dict_1[key] = data_dict_1[key][donor_idxs]
        for key in ['acceptor_atoms', 'acceptor_atom_names', 'acceptor_radii', 'acceptor_coords', 'acceptor_res','acceptor_root_angles', 'acceptor_chains', 'acceptor_resnames']:
            if key in data_dict_2:
                data_dict_2[key] = data_dict_2[key][acceptor_idxs]

    
        return distances, data_dict_1, data_dict_2

    def count_hydrogen_bonds(self, data_dict_1, data_dict_2, distances_dict):
        #Counts hydrogen bonds, respecting maximum possible bonds per pair, and categorizes them.

        total_hb = 0
        category_counts = {
            'mainchain_mainchain': 0,
            'mainchain_sidechain': 0,
            'sidechain_mainchain': 0,
            'sidechain_sidechain': 0
        }

        #mapping
        hb_list = []
        mainchain_atoms = ['N', 'CA', 'C', 'O']
        distances_arr = np.array(list(distances_dict.values()))

        # Create a list of tuples representing each donor-acceptor pair (ATOM SERIAL NUMBERS)
        #donor_acceptor_pairs = list(zip(data_dict_1['donor_atoms'], data_dict_2['acceptor_atoms']))

        # Count occurrences of each donor-acceptor pair
        #pair_counts = Counter(donor_acceptor_pairs)
        '''
        # Function to get residue number and atom name
        def get_residue_number_and_atom_name(atom_id, is_donor):
            """Efficiently retrieves residue number and atom name from atom ID."""
            if is_donor:
                try:
                    index = np.where(data_dict_1['donor_atoms'] == atom_id)[0][0]
                    res_num = data_dict_1['donor_res'][index]
                except IndexError:
                    print(f"ERROR! Could not find donor residue for atom ID {atom_id}")
                    return None, None
            else:  # is_acceptor
                try:
                    index = np.where(data_dict_2['acceptor_atoms'] == atom_id)[0][0]
                    res_num = data_dict_2['acceptor_res'][index]
                except IndexError:
                    print(f"ERROR! Could not find acceptor residue for atom ID {atom_id}")
                    return None, None

            # Now find the atom name
            for model in self.structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[1] == res_num:
                            for atom in residue:
                                if atom.get_serial_number() == atom_id:
                                    return res_num, atom.get_name()

            print(f"ERROR! Atom ID {atom_id} not found in residue {res_num}")
            return None, None  # Not found
        '''
        for idx, (donor_atom_id, acceptor_atom_id) in enumerate(zip(data_dict_1['donor_atoms'], data_dict_2['acceptor_atoms'])):
            donor_res_num = data_dict_1['donor_res'][idx]
            donor_atom_name = data_dict_1['donor_atom_names'][idx]
            donor_chain = data_dict_1['donor_chains'][idx]
            donor_resname = data_dict_1['donor_resnames'][idx]
            acceptor_res_num = data_dict_2['acceptor_res'][idx]
            acceptor_atom_name = data_dict_2['acceptor_atom_names'][idx]
            acceptor_chain = data_dict_2['acceptor_chains'][idx]
            acceptor_resname = data_dict_2['acceptor_resnames'][idx]

            max_donor_hbonds = self.max_hbonds.get(donor_atom_name)
            max_acceptor_hbonds = self.max_hbonds.get(acceptor_atom_name)
            if max_donor_hbonds is None or max_acceptor_hbonds is None:
                continue

            valid_hbonds = min(1, min(max_donor_hbonds, max_acceptor_hbonds))  # after filtering, each row is one HB
            total_hb += valid_hbonds

            donor_mc = donor_atom_name in mainchain_atoms
            acceptor_mc = acceptor_atom_name in mainchain_atoms
            if donor_mc and acceptor_mc:
                cat = 'mainchain_mainchain'
            elif donor_mc:
                cat = 'mainchain_sidechain'
            elif acceptor_mc:
                cat = 'sidechain_mainchain'
            else:
                cat = 'sidechain_sidechain'
            category_counts[cat] += valid_hbonds

            dist = distances_arr[idx]   # <-- NEW: get distance by index directly
            ang_d = data_dict_1.get('donor_root_angles', [np.nan]*len(distances_arr))[idx]
            ang_a = data_dict_2.get('acceptor_root_angles', [np.nan]*len(distances_arr))[idx]

            hb_list.append({
                'Donor_Chain': donor_chain,
                'Donor_Residue_Name': donor_resname,
                'Donor_Residue_Num': donor_res_num,
                'Acceptor_Chain': acceptor_chain,
                'Acceptor_Residue_Name': acceptor_resname,
                'Acceptor_Residue_Num': acceptor_res_num,
                'Donor_Atom': donor_atom_name,
                'Acceptor_Atom': acceptor_atom_name,
                'Distance': dist,
                'Donor_Angle': ang_d,
                'Acceptor_Angle': ang_a,
                'Type': cat.replace('_','-').replace('mainchain','MC').replace('sidechain','SD')
            })

        return total_hb, category_counts, pd.DataFrame(hb_list)
        '''
        for pair, count in pair_counts.items():
            # Find the donor and acceptor atoms for this pair
            donor_atom_id, acceptor_atom_id = pair  # ATOM IDs

            # Get donor and acceptor info
            donor_res_num, donor_atom_name = get_residue_number_and_atom_name(donor_atom_id, True)
            acceptor_res_num, acceptor_atom_name = get_residue_number_and_atom_name(acceptor_atom_id, False)

            # Handle cases where not found
            if donor_res_num is None or acceptor_res_num is None:
                continue

            # Determine maximum H-bonds based on atom types
            max_donor_hbonds = self.max_hbonds.get(donor_atom_name)
            max_acceptor_hbonds = self.max_hbonds.get(acceptor_atom_name)

            max_hbonds_per_pair = min(max_donor_hbonds, max_acceptor_hbonds)

            # Count only up to the maximum allowed
            total_hb += min(count, max_hbonds_per_pair)


            # Determine if mainchain or sidechain
            is_donor_mainchain = donor_atom_name in mainchain_atoms
            is_acceptor_mainchain = acceptor_atom_name in mainchain_atoms


            # Categorize the H-bond
            if is_donor_mainchain and is_acceptor_mainchain:
                category_counts['mainchain_mainchain'] += min(count, max_hbonds_per_pair)
            elif is_donor_mainchain and not is_acceptor_mainchain:
                category_counts['mainchain_sidechain'] += min(count, max_hbonds_per_pair)
            elif not is_donor_mainchain and is_acceptor_mainchain:
                category_counts['sidechain_mainchain'] += min(count, max_hbonds_per_pair)
            else:
                category_counts['sidechain_sidechain'] += min(count, max_hbonds_per_pair)

        return total_hb, category_counts
        '''
    def save_hb_structure(self, output_pdb_file):
        """
        Saves the hydrogen-bonded atoms to a new PDB file.
        Args:
        output_pdb_file (str): The name of the output PDB file.
        """
        def extract_hbond_atoms(structure):
            for model in structure:
                for chain in model:
                    for residue in chain:
                        atoms_to_remove = [atom for atom in residue 
                                           if atom.get_serial_number() not in self.receptor_data['donor_atoms'] 
                                           and atom.get_serial_number() not in self.receptor_data['acceptor_atoms']
                                           and atom.get_serial_number() not in self.ligand_data['donor_atoms']
                                           and atom.get_serial_number() not in self.ligand_data['acceptor_atoms']]
                        for atom in atoms_to_remove:
                            residue.detach_child(atom.id)
            return structure

        # Assuming you have already parsed your structure
        new_structure = extract_hbond_atoms(self.structure)
                 
        # Extract the atoms involved in hydrogen bonds
        hbond_atoms = extract_hbond_atoms(self.structure)

         # Save the new structure to a PDB file
        io = PDBIO()
        io.set_structure(hbond_atoms)
        io.save(output_pdb_file)
        
    def get_results(self):
        """Return a dictionary of all calculated results."""
        return {
            "Total HB": self.total_hb, 
            "MC-MC HB": self.category_counts['mainchain_mainchain'],
            "MC-SD HB": self.category_counts['mainchain_sidechain'],
            "SD-MC HB": self.category_counts['sidechain_mainchain'],
            "SD-SD HB": self.category_counts['sidechain_sidechain'],
            "Total RT": self.total_rotable_bond_count,
            "LMC-MC HB": self.category_counts_1['mainchain_mainchain'],
            "LMC-SD HB": self.category_counts_1['mainchain_sidechain'],
            "LSD-MC HB": self.category_counts_1['sidechain_mainchain'],
            "LSD-SD HB": self.category_counts_1['sidechain_sidechain'],
            "RMC-MC HB": self.category_counts_2['mainchain_mainchain'],
            "RMC-SD HB": self.category_counts_2['mainchain_sidechain'],
            "RSD-MC HB": self.category_counts_2['sidechain_mainchain'],
            "RSD-SD HB": self.category_counts_2['sidechain_sidechain'],            
            "Ligand RT": self.total_rotable_bond_count_1,
            "Receptor RT": self.total_rotable_bond_count_2,
            "LDonor RT": self.ligand_data['donor_root_rotable_bond_count'],
            "RDonor RT": self.receptor_data['donor_root_rotable_bond_count'],
            "LAcceptor RT": self.ligand_data['acceptor_root_rotable_bond_count'],
            "RAcceptor RT": self.receptor_data['acceptor_root_rotable_bond_count']
            
        }
    
    def export_hb(self, outname):
        if not hasattr(self, "hb_df1") or not hasattr(self, "hb_df2"):
            #CONTINUE HERE CONVERT DICTIONARY INTO PANDAS DF AND EXPORT 
            #hb_df = self.empty_hb_list
            hb_df = pd.DataFrame([self.empty_hb_list])
            print("Hydrogen bond DataFrames not found. Creating empty result.")
        else:
            # Save both DataFrames
            hb_df = pd.concat([self.hb_df1, self.hb_df2])
        hb_df.to_csv(outname, index=False)


def extract_features(pdb_file, ligand_chains, receptor_chains, out_preffix):
    key_list = []
    value_list = []
    
    parser = PDBParser()
    structure = parser.get_structure("protein_name", pdb_file)

    hc_analyzer = HCAnalyzer(structure, ligand_chains, receptor_chains)
    hc_analyzer()
    results = hc_analyzer.get_results()
    for key,value in results.items():
        key_list.append(key)
        value_list.append(value)
        
    hat_analyzer = HATAnalyzer(structure, ligand_chains, receptor_chains)
    hat_analyzer()
    results = hat_analyzer.get_results()
    for key,value in results.items():
        key_list.append(key)
        value_list.append(value)

    hb_analyzer = HBAnalyzer(structure, ligand_chains, receptor_chains, distance_threshold=1.4)
    try:
        hb_analyzer() 
        #hb_analyzer.save_hb_structure("test.pdb")
        results = hb_analyzer.get_results()
    except Exception as e:
        results = hb_analyzer.empty_results
        pass
    for key,value in results.items():
        key_list.append('W'+key)
        value_list.append(value)     
    if out_preffix != None:
        hb_analyzer.export_hb(out_preffix+"_WHB.csv")

    hb_analyzer = HBAnalyzer(structure, ligand_chains, receptor_chains)
    try:
        hb_analyzer() 
        #hb_analyzer.save_hb_structure("test.pdb")
        results = hb_analyzer.get_results()
    except Exception as e:
        results = hb_analyzer.empty_results
        pass
    for key,value in results.items():
        key_list.append(key)
        value_list.append(value)
    if out_preffix != None:
        hb_analyzer.export_hb(out_preffix+"_MHB.csv")

    hb_analyzer = HBAnalyzer(structure, ligand_chains, receptor_chains, distance_threshold=0)
    try:
        hb_analyzer() 
        results = hb_analyzer.get_results()
    except Exception as e:
        results = hb_analyzer.empty_results
        pass
    for key,value in results.items():
        key_list.append('S'+key)
        value_list.append(value)        
    if out_preffix != None:
        hb_analyzer.export_hb(out_preffix+"_SHB.csv")

    vdw_analyzer = VDWAnalyzer(structure, ligand_chains, receptor_chains)
    vdw_analyzer()
    results = vdw_analyzer.get_results()
    for key,value in results.items():
        key_list.append(key)
        value_list.append(value)
    
    return np.array(key_list), np.array(value_list)


if __name__ == "__main__":
    import sys
    out_preffix = None
    if len(sys.argv) >= 5:
        out_preffix = sys.argv[4]
    if len(sys.argv) >= 4:
        keys, values = extract_features(sys.argv[1], sys.argv[2].split(','), sys.argv[3].split(','), out_preffix)
        for k, v in zip(keys,values):
            print(k,v)
    else:
        space="      "
        print(space+"Usage:")
        print(space+space+"python3 pptools.py <input> <protein1_chain> <protein2_chain> <optional_out_preffix>")
        print("")
        print(space+"Examples:")
        print(space+space+"python3 pptools.py proteins.pdb A B")
        print(space+space+"python3 pptools.py proteins.pdb A,B C,D,E my_output")
        sys.exit()

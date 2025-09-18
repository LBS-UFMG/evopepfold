import pymol
from pymol import stored
import pandas as pd
import pyrosetta
import numpy as np
from itertools import combinations
from time import time

# Avoiding a very cluttered output
pyrosetta.init("-mute all")

PROTOCOL = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="fa_standard" weights="ref2015.wts"/>
    </SCOREFXNS>
    <MOVERS>
        <FlexPepDock name="ppack"  ppk_only="true"/>
        <FlexPepDock name="lowres" lowres_abinitio="true" pep_refine="true"/>
        <FlexPepDock name="minimize"  min_only="true"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover="ppack"/>
        <Add mover="lowres"/>
        <Add mover="minimize"/>
    </PROTOCOLS>
    <OUTPUT/>
    </ROSETTASCRIPTS>""").get_mover("ParsedProtocol")

FULL_PROTOCOL = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="fa_standard" weights="ref2015.wts"/>
    </SCOREFXNS>
    <MOVERS>
        <FlexPepDock name="ppack"  ppk_only="true"/>
        <FlexPepDock name="lowres" lowres_abinitio="true" pep_refine="false"/>
        <FlexPepDock name="fpd" lowres_abinitio="false" pep_refine="true"/>
        <FlexPepDock name="minimize"  min_only="true"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover="ppack"/>
        <Add mover="lowres"/>
        <Add mover="fpd"/>
        <Add mover="minimize"/>
    </PROTOCOLS>
    <OUTPUT/>
    </ROSETTASCRIPTS>""").get_mover("ParsedProtocol")


def evaluateDiversity(sequences):
    """
    Computes the mean pairwise normalized Hamming distance 
    between all sequences in the list.
    
    Parameters:
    sequences (list of str): List of amino acid sequences (all same length).
    
    Returns:
    float: Average diversity (normalized Hamming distance).
    """
    if not sequences:
        return 0.0
    
    n = len(sequences)
    sequence_length = len(sequences[0])
    
    total_distance = 0
    pair_count = 0
    
    for seq1, seq2 in combinations(sequences, 2):
        # Ensure sequences have the same length
        if len(seq1) != sequence_length or len(seq2) != sequence_length:
            raise ValueError("All sequences must have the same length.")
        
        # Compute Hamming distance
        hamming_distance = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
        
        # Normalize by sequence length
        normalized_distance = hamming_distance / sequence_length
        
        total_distance += normalized_distance
        pair_count += 1
    
    if pair_count == 0:
        return 0.0  # Only one individual or empty list
    
    return total_distance / pair_count

def countContacts(candidate_file, target = None, cutoff = 5.0):
    if target == None:
        target = [('20', 'VAL'), ('21', 'THR'), ('23', 'GLY'), ('24', 'THR'), ('25', 'THR'), 
        ('26', 'THR'), ('27', 'LEU'), ('41', 'ALA'), ('42', 'VAL'), ('44', 'CYS'), ('45', 'THR'),
        ('46', 'ALA'), ('49', 'MET'), ('50', 'LEU'), ('52', 'PRO'), ('67', 'LEU'), ('69', 'GLN'), 
        ('118', 'TYR'), ('119', 'ASN'), ('139', 'SER'), ('140', 'PHE'), ('141', 'LEU'), ('142', 'ASN'), 
        ('143', 'GLY'), ('144', 'SER'), ('145', 'CYS'), ('146', 'GLY'), ('163', 'HIS'), ('164', 'HIS'), 
        ('165', 'MET'), ('166', 'GLU'), ('167', 'LEU'), ('168', 'PRO'), ('172', 'HIS'), ('173', 'ALA'), 
        ('186', 'VAL'), ('187', 'ASP'), ('188', 'ARG'), ('189', 'GLN'), ('190', 'THR'), ('191', 'ALA'), 
        ('192', 'GLN'), ('193', 'ALA')]
    pymol.cmd.reinitialize()
    stored.list = []
    pymol.cmd.load(candidate_file)
    pymol.cmd.indicate(''.join(['chain A within ', str(cutoff), ' of chain B']))
    pymol.cmd.iterate('indicate', 'stored.list.append((resi,resn))')
    pymol.cmd.delete('indicate')
    count = 0
    for aa in set(stored.list):
        if aa in target:
            count += 1
    occupancy = count/len(target)
    return occupancy


def calculateFitness_Parallel(candidate, protocol, scorefxn, gen, reps, queue):
    best_pose = None
    best_fitness = 999999999
    scores = list()
    time_start = time()
    for i in range(reps):
        # Load the poses
        comp = pyrosetta.pose_from_pdb("Results/PDB_complex/G" + str(gen) + "/comp_" + candidate + ".pdb")
        # Apply the docking protocol and calculate fitness
        protocol.apply(comp)
        fitness = scorefxn(comp)
        scores.append(fitness)
        if fitness < best_fitness:
            best_pose = comp.clone()
            best_fitness = fitness
    best_pose.dump_pdb("Results/PDB_docking/G" + str(gen) + "/dock_" + candidate + ".pdb")
    time_end = time()
    time_elapsed = time_end-time_start
    queue.put([candidate, scores, best_fitness, time_elapsed])
    


def calculateFitness_fake(structure_path, score_function, output_folder):
    scorefxn = "ref2015"
    return np.random.randint(-1000, 1000)

def countContacts_fake(candidate_file, target = None, cutoff = 5.0):
    return np.random.randint(0,50)
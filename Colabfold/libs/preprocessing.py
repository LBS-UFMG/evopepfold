#import pymol
import os
import numpy as np
#import pyrosetta
import pandas as pd

# Avoiding a very cluttered output
#pyrosetta.init("-mute all")

def makeFasta(receptor, peptide, fasta_folder):
    fasta_path = os.path.join(fasta_folder, f'{peptide}.fasta')
    with open(fasta_path, 'w+') as fasta:
        fasta.write(f'{receptor}:{peptide}')
    return fasta_path

def makePeptide(candidate, gen=0):
    peptide = ['abc']#pyrosetta.pose_from_sequence(candidate)
    angles = np.random.rand(3,len(candidate))*360 - 180
    for i in range(len(candidate)):
        peptide.set_phi(i+1, angles[0,i])
        peptide.set_psi(i+1, angles[1,i])
        #peptide.set_chi(1, i+1, angles[2,i])   
    src = 'Temp/' + candidate + ".pdb"
    peptide.dump_pdb(src)
    dst = "Results/PDB_predicted/G" + str(gen) + "/" + candidate + ".pdb"
    #pymol.cmd.load("template_pos.pdb")
    #pymol.cmd.load(src, "new")
    #pymol.cmd.align("new", "template_pos")
    #pymol.cmd.save(dst, selection = "new")
    #pymol.cmd.reinitialize()
    os.remove('Temp/' + candidate + ".pdb")
    makeComplex(candidate, gen)


def makeComplex(candidate, gen):
    # pymol.cmd.load("target.pdb")
    # pymol.cmd.load("Results/PDB_predicted/G" + str(gen) + "/" + candidate + ".pdb", "pep")
    # #pymol.cmd.alter('chain A', 'chain="B"')
    # pymol.cmd.alter('pep', 'chain="B"')
    # #pymol.cmd.alter('chain E', 'chain="A"') # The rename is necessary because Rosetta has a problem in making chain A move
    # pymol.cmd.save("Results/PDB_complex/G" + str(gen) + "/comp(" + candidate + ").pdb") # Join both chains in a single PDB as a complex
    # pymol.cmd.reinitialize()
    pass

def makeBatch(population, cores, add_fake = False):
    # Population is the list of strings representing the peptides
    # n is the number of peptides per batch
    batches = list()
    x = int(len(population)/cores)
    if len(population) < cores:
        return [population]
    
    for i in range(x):
        if add_fake == True:
            new_batch = ["BJOUXZ"]
        else:
            new_batch = list()
        for j in range(cores):
            pos = i*cores + j
            new_batch.append(population[pos])
        batches.append(new_batch)
        
    if len(population)%cores != 0: # If true, there is an excess of peptides not covered yet
        if add_fake == True:
            final_batch = ["BJOUXZ"]
        else:
            final_batch = list()
        
        for k in range(len(population)%cores):
            pos = x*cores + k
            final_batch.append(population[pos])
        batches.append(final_batch)
    return batches

def makeSummary(pr):
    generation = []
    best_ind = []
    best_fit = []
    avg_fit = []
    worst_fit = []
    total_time = []
    oc = []
    max_oc = []
    for gen in pr['generation'].unique(): 
        gen_pr = pr.loc[pr['generation']==gen]
        sorted_pr = gen_pr.sort_values(by='best', ascending=True).reset_index(drop=True)
        best = sorted_pr.iloc[0]
        generation.append(best['generation'])
        best_ind.append(best['genotype'])
        best_fit.append(best['best'])
        avg_fit.append(sorted_pr['best'].mean())
        worst_fit.append(sorted_pr['best'].max())
        total_time.append(99999)
        oc.append(best['occupancy'])
        max_oc.append(sorted_pr['occupancy'].max())

    sm = {'generation': generation,
    'best_ind': best_ind,
    'best_gen_fit': best_fit,
    'avg_gen_fit': avg_fit,
    'worst_gen_fit': worst_fit,
    'total_time': total_time,
    '%surface_occupied': oc,
    'maximum_occupancy': max_oc}

    return pd.DataFrame(sm)
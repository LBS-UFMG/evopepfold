#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 6 18:57:14 2021
Modified on Thu Feb 20 22:04:50 2025

@author: FCCarvalho
"""

import os
import time
import random as rd
import numpy as np
import joblib
from multiprocessing import Process, Queue
import shutil
import pandas as pd
from libs.operations import mutation, deletion, insertion, crossover, tournament
from libs.output import makeDirectories, makeOutputs, prepareDirectories, transferPrevious
from libs.evaluations import countContacts  # using the real evaluation function
from libs.preprocessing import makeFasta, makeSummary

# NEW: Import and initialize PyRosetta with the ref2015 scoring function
import pyrosetta
pyrosetta.init("-mute all")
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
scorefxn = ScoreFunctionFactory.create_score_function("ref2015")


COLABFOLD_PYTHON_PATH = "/hd/localcolabfold/colabfold-conda/bin/python3"
COLABFOLD_PATH = "/hd/localcolabfold/colabfold-conda/bin/colabfold_batch"

def GA(initial_population, tournament_size=3, max_generations=10, pop_size=None, nsims=50, mutation_rate=0.20, cross_rate=0.80, elite=1, elite_limit=5, workers=1):
    makeDirectories()
    previous_result = makeOutputs(nsims)
    removed = []
    initial = 0
    current_elite = ''
    current_elite_count = 0

    if pop_size is None:
        pop_size = len(initial_population)

    population = initial_population
    if previous_result:
        pr = pd.read_csv("Results/populations_GA.csv")
        fit_dict = {}
        ocup_dict = {}
        known_genotypes = pr["genotype"].tolist()
        known_scores = pr["best"].tolist()
        known_occupancies = pr["occupancy"].tolist()
        for i in range(len(known_genotypes)):
            fit_dict[known_genotypes[i]] = known_scores[i]
            ocup_dict[known_genotypes[i]] = known_occupancies[i]
        sm = pd.read_csv("Results/summary_GA.csv")
        if len(sm) == 0:
            sm = makeSummary(pr)
            sm.to_csv("Results/summary_GA.csv")
        elif max(sm['generation']) < max(pr['generation']): # Verificar isso
            sm = makeSummary(pr)
            sm.to_csv("Results/summary_GA.csv")
        best_genotypes = sm["best_ind"].tolist()
        best_fitness = sm["best_gen_fit"].tolist()
        average_fitness = sm["avg_gen_fit"].tolist()
        worst_fitness = sm["worst_gen_fit"].tolist()
        initial = sm["generation"].max()
        print("Imported", len(fit_dict), "known peptides")
        print("Starting from generation", initial+1)
        print("#" * 10, "\n" * 2)
        max_occupancy = []
        best_occupancy = sm['%surface_occupied'].to_list()
        population = pr.loc[pr['generation'] == initial, ['genotype']]['genotype'].to_list()
        occupancies = pr.loc[pr['generation'] == initial, ['occupancy']]['occupancy'].to_list()
        batch_pop = population.copy()
        logs = open('logs.txt', 'a')
        logs.write(f'\n\nSimulation resumed from generation {initial +1}\n!')
        logs.write(f'{len(fit_dict)} known peptides were imported.\n\n')
    else:
        best_fitness = []
        known_genotypes = []
        known_scores =[]
        known_occupancies = []
        best_genotypes = []
        worst_fitness = []
        average_fitness = []
        max_occupancy = []
        best_occupancy = []
        fit_dict = {}
        ocup_dict = {}
        occupancies = []
        logs = open('logs.txt', 'w+')
        logs.write(f'Simulation started!')
    z = 0
    
    for j in range(initial, initial + max_generations):
        if not previous_result:
            print("Beginning generation", j)
            prepareDirectories(j)
            start = time.time()

            # Remove repeated peptides from previous generations
            sm = pd.read_csv("Results/summary_GA.csv")
            best_peps = list(sm['best_ind'])
            for b in pd.Series(best_peps).unique():
                if best_peps.count(b) > 3 and pd.unique(sm[sm['best_ind'] == b]['best_gen_fit'])[0] < -100:
                    removed.append(b)

            batch_pop = pd.Series(population).unique().tolist()
            if not previous_result:
                for genotype in reversed(batch_pop):
                    if genotype in fit_dict.keys() and not previous_result and genotype not in removed:
                        with open("Results/populations_GA.csv", mode="a") as p_info:
                            p_info.write("\n" + str(j) + "," + genotype)
                            for r in range(nsims):
                                p_info.write("," + str(fit_dict[genotype]))
                            p_info.write("," + str(fit_dict[genotype]) + "," + str(best_occupancy[-1]))
                        #transferPrevious(genotype, j)
                        batch_pop.remove(genotype)
                

                # ===== NEW BLOCK: Write FASTA files, run ColabFold, and manage outputs =====

                # Create (or clear) the fasta directory
                if os.path.exists("fasta"):
                    shutil.rmtree("fasta")
                os.makedirs("fasta")

                # Write a FASTA file for each peptide (genotype)
                for genotype in batch_pop:
                    fasta_file = os.path.join("fasta", f"{genotype}.fasta")
                    with open(fasta_file, "w") as out_f:
                        #out_f.write(">peptideo\n")
                        #out_f.write(f"{genotype}\n")
                        out_f.write(f">{genotype}\n")
                        out_f.write(f"SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQ:{genotype}\n")

                # Ensure the new_output directory is clean
                if os.path.exists("new_output"):
                    shutil.rmtree("new_output")
                os.makedirs("new_output")

                # Execute the ColabFold batch command (os.system waits until completion—there is no built-in timeout)
                #os.system("colabfold_batch --num-relax 1 --use-gpu-relax fasta/ new_output/") #AQUI
                #for genotype in population: # DEBUG ONLY
                #    src = 'output/pep1.pdb'
                #    shutil.copy(src, f"new_output/{genotype}_relaxed_rank_001.pdb") # FAKE
                for genotype in batch_pop:
                    #os.system(f"colabfold_batch --num-relax 1 --use-gpu-relax fasta/{genotype}.fasta new_output/") 
                    if genotype != current_elite:
                        logs.write(f'\nModelling peptide {genotype}:')
                        before = time.time()
                        os.system(f"{COLABFOLD_PYTHON_PATH} {COLABFOLD_PATH} --num-relax 1 --num-recycle 3 --use-gpu-relax --msa-mode mmseqs2_uniref_env_envpair --pair-mode unpaired --model-type alphafold2_multimer_v3 fasta/{genotype}.fasta new_output/")
                        after = time.time()
                        logs.write(f'\nPeptide {genotype} took {int(after - before)} seconds to complete.')
                    
                # Make a complete copy of new_output for future research purposes in the 'output' directory
                dest_output = os.path.join("output", f"G{j}")
                if os.path.exists(dest_output):
                    shutil.rmtree(dest_output)
                shutil.copytree("new_output", dest_output)

                # Create destination folder for evaluation results in the Results folder
                dest_results = os.path.join("Results", f"G{j}")
                if not os.path.exists(dest_results):
                    os.makedirs(dest_results)

                if current_elite != '' and j >0:
                    src = os.path.join("Results", f"G{j-1}",f"dock_{current_elite}.pdb")
                    dst = os.path.join("Results", f"G{j}",f"dock_{current_elite}.pdb")
                    shutil.copy(src, dst)
                    
                # Process only the "relaxed_rank_001" PDB files:
                # For each file in new_output that matches the pattern, extract the peptide name and move it with the new naming format.
                for file in os.listdir("new_output"):
                    if file.endswith(".pdb") and "relaxed_rank_001" in file:
                        # Assuming the filename is formatted as {peptide}_relaxed_rank_001_...pdb
                        genotype = file.split("_relaxed_rank_001")[0]
                        if genotype in batch_pop:
                            new_name = f"dock_{genotype}.pdb"
                            src_file = os.path.join("new_output", file)
                            dest_file = os.path.join(dest_results, new_name)
                            shutil.move(src_file, dest_file)
                # =========================================================================

                # ===== Evaluate structures using PyRosetta and ref2015 scoring =====
                with open("Results/populations_GA.csv", mode="a") as p_info:
                    for pep in batch_pop:
                        if pep != current_elite:
                            pdb_file = os.path.join(dest_results, f"dock_{pep}.pdb")
                            # Load structure and compute score using ref2015
                            pose = pose_from_pdb(pdb_file)
                            s = scorefxn(pose)
                            oc = countContacts(pdb_file)
                            #oc = 0.9
                            scores_txt = ",".join([str(s)] * nsims)
                            p_info.write("\n" + str(j) + "," + pep + "," + scores_txt + "," + str(s) + "," + str(oc))
                            if oc < 0.3:
                                s = 9999999999  # Penalize peptides with low occupancy
                            fit_dict[pep] = s
                            ocup_dict[pep] = oc
                            occupancies.append(oc)
                        else:
                            occupancies.append(ocup_dict[pep])
                # =========================================================================

            fit_list = [fit_dict[genotype] for genotype in batch_pop]
            if len(batch_pop) > 0:
                best_fitness.append(min(fit_list))
                max_occupancy.append(max(occupancies))
                l = fit_list.index(min(fit_list))
                best_genotypes.append(batch_pop[l])
                best_occupancy.append(ocup_dict[best_genotypes[-1]])
                fit_list2 = [f for f in fit_list if f < 999999999]
                if len(fit_list2) > 0:
                    worst_fitness.append(max(fit_list2))
                    average = sum(fit_list2) / len(fit_list2)
                else:
                    worst_fitness.append(999999999)
                    average = np.mean(fit_list)
                average_fitness.append(average)
                known_genotypes.extend(batch_pop)
        else:
            start = time.time()
            fit_list = [fit_dict[genotype] for genotype in batch_pop]

        # ===== Create next generation =====
        spots_left = pop_size
        new_pop = []
        if elite > 0:
            ranking = [x for _, x in sorted(zip(fit_list, batch_pop))]
            new_elite = ranking[0] # O Rank 0 é o de menor score, logo, o de melhor score
            if new_elite == current_elite and current_elite_count > elite_limit:
                elite_range = 0
                current_elite_count = 0
            else:
                elite_range = elite
            for n in range(elite_range):
                new_pop.append(ranking[n])
            spots_left -= elite_range
            current_elite = new_elite
            current_elite_count +=1 
        operation_weights = [cross_rate, mutation_rate/3, mutation_rate/3, mutation_rate/3]
        while spots_left > 0:
            if spots_left == 1:
                operation = rd.choice(["insertion", "deletion", "mutation"])
            else:
                operation = rd.choices(["crossover", "insertion", "deletion", "mutation"], weights=operation_weights)[0]
            if operation == "crossover":
                parent1 = tournament(batch_pop, fit_list, tournament_size, reuse=True)
                parent2 = tournament(batch_pop, fit_list, tournament_size, reuse=True)
                child1, child2 = crossover(parent1, parent2)
                if child1 not in known_genotypes:
                    new_pop.append(child1)
                    spots_left -= 1
                if child2 not in known_genotypes:
                    new_pop.append(child2)
                    spots_left -= 1
            elif operation == "insertion":
                chosen = tournament(batch_pop, fit_list, tournament_size, reuse=False)
                mutant = insertion(chosen) if len(chosen) < 50 else deletion(chosen)
                if mutant not in known_genotypes:  
                    new_pop.append(mutant)
                    spots_left -= 1
            elif operation == "deletion":
                chosen = tournament(batch_pop, fit_list, tournament_size, reuse=False)
                mutant = deletion(chosen) if len(chosen) > 5 else insertion(chosen)
                if mutant not in known_genotypes:
                    new_pop.append(mutant)
                    spots_left -= 1
            else:
                chosen = tournament(batch_pop, fit_list, tournament_size, reuse=False)
                mutant = mutation(chosen)
                if mutant not in known_genotypes:
                    new_pop.append(mutant)
                    spots_left -= 1

        end = time.time()
        gen_time = end - start
        if len(batch_pop) > 0 and not previous_result:
            with open("Results/summary_GA.csv", mode="a") as summary:
                summary.write(str(j) + "," + str(best_genotypes[-1]) + "," + str(best_fitness[-1]) + "," +
                              str(average) + "," + str(worst_fitness[-1]) + "," + str(gen_time) + "," +
                              str(best_occupancy[-1]) + "," + str(max_occupancy[-1]) + "\n")
        population = new_pop.copy()
        if not previous_result:
            time.sleep(600)
        previous_result = False
    logs.close()
    return population, best_fitness, worst_fitness, average_fitness, best_genotypes, max_occupancy, best_occupancy

if __name__ == '__main__':
    s = time.time()
    #initial_population = ['SEQUENCIA']
    initial_population = joblib.load("best100.lst")
    p, bf, wf, af, bi, mo, bo = GA(initial_population, max_generations=100, nsims=5, pop_size=50, workers=25, tournament_size=2, mutation_rate=0.1, cross_rate=0.9, elite=1, elite_limit=3)
    e = time.time()
    joblib.dump(p, "Results/populations.lst")
    joblib.dump(bf, "Results/best_scores.lst")
    joblib.dump(wf, "Results/worst_scores.lst")
    joblib.dump(af, "Results/avg_scores.lst")
    joblib.dump(bi, "Results/best_genotypes.lst")
    joblib.dump(mo, "Results/max_occupancies.lst")
    joblib.dump(bo, "Results/best_occupancies.lst")
    print('\n\n\n')
    print('FINISHED IN', np.round((e-s)/60, 3), 'MINUTES')
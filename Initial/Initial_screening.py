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
from libs.evaluations import countContacts, PROTOCOL, calculateFitness_Parallel
from libs.preprocessing import makeComplex, makePeptide, makeBatch, makeSummary

# NEW: Import and initialize PyRosetta with the ref2015 scoring function
import pyrosetta
pyrosetta.init("-mute all")
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
scorefxn = ScoreFunctionFactory.create_score_function("ref2015")


COLABFOLD_PYTHON_PATH = "/hd/localcolabfold/colabfold-conda/bin/python3"
COLABFOLD_PATH = "/hd/localcolabfold/colabfold-conda/bin/colabfold_batch"

def GA(initial_population, tournament_size=3, max_generations=10, pop_size=None, nsims=50, mutation_rate=0.20, cross_rate=0.80, elite=1, workers=1):
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
    else:
        best_fitness = []
        best_genotypes = []
        worst_fitness = []
        average_fitness = []
        max_occupancy = []
        best_occupancy = []
        fit_dict = {}
        ocup_dict = {}
        occupancies = []
        known_genotypes = []
    
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

            # Start Queue 
            q = Queue()
        
            batch_pop = pd.Series(population).unique().tolist()
            for genotype in reversed(batch_pop):
                if genotype in fit_dict.keys() and not previous_result and genotype not in removed:
                    with open("Results/populations_GA.csv", mode="a") as p_info:
                        p_info.write("\n" + str(j) + "," + genotype)
                        for r in range(nsims):
                            p_info.write("," + str(fit_dict[genotype]))
                        p_info.write("," + str(fit_dict[genotype]) + "," + str(best_occupancy[-1]))
                    #transferPrevious(genotype, j)
                    batch_pop.remove(genotype)
            

            batches = makeBatch(batch_pop, cores = workers)
            occupancies = list()
            previous_result = False
            # Preparation step is sequential because of Modeller
            z = 0
            
            for batch in batches:
                print("######## BATCH", z, "########\n\n")
                for genotype in batch:
                    if genotype != "BJOUXZ": # The fake peptide speeds up simulation since the first process starts alone
                        #if not os.path.isdir("Results/PDB_docking/G0/" + genotype):
                        #    os.mkdir("Results/PDB_docking/G0/" + genotype)
                        #preparePeptide(candidate=genotype, gen=j)
                        makePeptide(candidate=genotype, gen=j)
                
                # Start the process list
                processes = list()
                for genotype in batch: # The docking step is parallel
                    p = Process(target=calculateFitness_Parallel, args=(genotype, PROTOCOL, scorefxn, j, nsims, q))
                    p.start()
                    processes.append(p)

                for p in processes: # We can only move forward after all the processes finish
                    p.join()
                
                print("##### END OF BATCH", z, "#####")
                print("\n", "#"*10, "Queue Size:", q.qsize(), "#"*10, "\n"*2)
                z += 1
                p_info = open("Results/populations_GA.csv", mode="a")
                for i in range(q.qsize()): # Fetching the scores of each genotype
                    [pep, s, bf] = q.get()
                    oc = countContacts("Results/PDB_docking/G" + str(j) + "/dock_" + pep + ".pdb")
                    scores_txt = [str(score) for score in s]
                    scores_txt = ",".join(scores_txt)
                    p_info.write("".join(["\n", str(j), ",", pep, ",", scores_txt, ",", str(bf), ",", str(oc)]))
                    if oc < 0.3:
                        s = [9999999999, 9999999999] # Make it so the peptide with less than 30% occupancy gets bad scores 
                    fit_dict[pep] = min(s) # The best of all the scores is selected as true score
                    ocup_dict[pep] = oc
                    if pep not in fit_dict.keys():
                        occupancies.append(oc)
                    else:
                        occupancies.append(ocup_dict[pep])
                p_info.close()
                time.sleep(15)
                
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
            if new_elite == current_elite and current_elite_count > 2:
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
        previous_result = False
        
    return population, best_fitness, worst_fitness, average_fitness, best_genotypes, max_occupancy, best_occupancy

if __name__ == '__main__':
    s = time.time()
    initial_population = joblib.load("initial.lst")
    p, bf, wf, af, bi, mo, bo = GA(initial_population, max_generations=1, nsims=5, pop_size=80, workers=40)
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
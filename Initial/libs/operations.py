import random as rd

# Implementation of the four chosen genetic operations
# First operation: Mutation
def mutation(peptide):
    aminoacids = ['A','C','D','E','F','G','H','I','L','M','N','P','Q','R','S','T','V','W','Y'] 
    pep = list(peptide)
    l = len(pep) - 1
    point = rd.randint(0,l)
    while True:
        new_aminoacid = rd.choice(aminoacids)
        if new_aminoacid != pep[point]:
            break
    pep[point] = new_aminoacid
    pep = "".join(pep)
    return pep

# Second operation: Deletion
def deletion(peptide):
    pep = list(peptide)
    l = len(pep) - 1
    point = rd.randint(0,l)
    pep.pop(point)
    pep = "".join(pep)
    return pep

# Third operation: Addition
def insertion(peptide):
    aminoacids = ['A','C','D','E','F','G','H','I','L','M','N','P','Q','R','S','T','V','W','Y']
    pep = list(peptide)
    l = len(pep)
    point = rd.randint(0,l)
    new_aminoacid = rd.choice(aminoacids)
    pep.insert(point, new_aminoacid)
    pep = "".join(pep)
    return pep

# Fourth operation: Crossover
def crossover(peptide1, peptide2, l_limit = 5, u_limit = 30):
    pep1 = list(peptide1)
    pep2 = list(peptide2)
    l1 = len(pep1) - 1
    l2 = len(pep2) - 1
    point1 = rd.randint(1,l1)
    point2 = rd.randint(1,l2)
    if point1 + (l2 - point2 + 1) > u_limit:
        point1 = u_limit - (l2 - point2 + 1)
    elif point1 + (l2 - point2 + 1) < l_limit:
        point1 = l_limit - (l2 - point2 + 1)
    if point2 + (l1 - point1 + 1) > u_limit:
        point2 = u_limit - (l1 - point1 + 1)
    elif point2 + (l1 - point1 + 1) < l_limit:
        point2 = l_limit - (l1 - point1 + 1)
    pep1_l = pep1[:point1]
    pep1_r = pep1[point1:]
    pep2_l = pep2[:point2]
    pep2_r = pep2[point2:]
    pep1 = pep1_l + pep2_r
    pep2 = pep2_l + pep1_r
    pep1 = "".join(pep1)
    pep2 = "".join(pep2)
    return pep1, pep2

# Fifth operation: tournament
def tournament(population, fitness, size = 5, reuse = False):
    # Randomly selects tournament participants among members of the population
    selected = []
    scores = []
    if size >= len(population):
        selected = population.copy()
        scores = fitness.copy()
    else:
        for i in range(size):
            new = rd.choice([p for p in population if p not in selected])
            selected.append(new)
            j = population.index(new)
            scores.append(fitness[j])
    ranking = [x for _,x in sorted(zip(scores,selected))]
    if not reuse:
        j = population.index(ranking[0])
        population.remove(ranking[0])
        fitness.pop(j)
    return ranking[0]
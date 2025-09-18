import os
import shutil

def makeDirectories():
# Create directories
    if not os.path.isdir("Results/"):
            os.mkdir("Results/")
    if not os.path.isdir("Temp/"):
            os.mkdir("Temp/")
    if not os.path.isdir("Results/PDB_docking/"):
            os.mkdir("Results/PDB_docking/")
    if not os.path.isdir("Results/PDB_predicted/"):
            os.mkdir("Results/PDB_predicted/")
    if not os.path.isdir("Results/PDB_complex/"):
            os.mkdir("Results/PDB_complex/")

def makeOutputs(nsims):
    previous_result = False
    if not os.path.isfile("Results/summary_GA.csv"):
        summary = open("Results/summary_GA.csv", mode="w+") # Create the summary csv file for later visualization of results
        summary.write("generation,best_ind,best_gen_fit,avg_gen_fit,worst_gen_fit,total_time,%surface_occupied,maximum_occupancy" + "\n") # The headers of summary
        summary.close()
    if not os.path.isfile("Results/populations_GA.csv"):
        p_info = open("Results/populations_GA.csv", mode="w+") # Create populations csv file for posterior checks
        p_info.write("generation,genotype,") # The headers of summary
        for i in range(nsims):
            p_info.write("f{},".format(i+1))
        p_info.write("best,occupancy") 
        p_info.close()   
    else:
        previous_result = True
        print("\n"*2, "PREVIOUS RESULTS FILE FOUND. CONTINUING FROM LAST RESULTS.", "\n"*2)
    return previous_result


def prepareDirectories(i):
    # Prepare the directories
    if not os.path.isdir("Results/PDB_docking/G" + str(i)):
        os.mkdir("Results/PDB_docking/G" + str(i))
    if not os.path.isdir("Results/PDB_predicted/G" + str(i)):
        os.mkdir("Results/PDB_predicted/G" + str(i))
    if not os.path.isdir("Results/PDB_complex/G" + str(i)):
        os.mkdir("Results/PDB_complex/G" + str(i))

def transferPrevious(genotype, j):
     for r in reversed(range(j)):
        if os.path.isfile("Results/PDB_docking/G" + str(r) + "/dock(" + genotype + ").pdb"):
            src1 = "Results/PDB_docking/G" + str(r) + "/dock(" + genotype + ").pdb"
            dst1 = "Results/PDB_docking/G" + str(j) + "/dock(" + genotype + ").pdb"
            src2 = "Results/PDB_predicted/G" + str(r) + "/" + genotype + ".pdb"
            dst2 = "Results/PDB_predicted/G" + str(j) + "/" + genotype + ".pdb"
            src3 = "Results/PDB_complex/G" + str(r) + "/comp(" + genotype + ").pdb"
            dst3 = "Results/PDB_complex/G" + str(j) + "/comp(" + genotype + ").pdb"
            shutil.copyfile(src1, dst1)
            shutil.copyfile(src2, dst2)
            shutil.copyfile(src3, dst3)
            break
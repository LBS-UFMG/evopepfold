<h1  align="center">
  EvoPepFold
</h1 >

## Overview
This repository contains code that depends on PyRosetta and Local ColabFold and is intended to run on a machine with substantial compute resources. Individual runs typically take **several days**. This README explains how to prepare the environment, install dependencies, and run the code safely in the background on a Linux server.

## Quick summary
- We recommend using **Miniconda** to manage environments.
- Install Python packages from `requirements.txt`.
- Install **PyRosetta** following the official instructions: https://www.pyrosetta.org/downloads
- Install **Local ColabFold** following the official instructions: https://github.com/YoshitakaMo/localcolabfold
- To run a version of the code, `cd` into that version's folder and run the Python file there.
- Runs usually take several days — run them with `nohup`, `tmux`, `screen`, or a job scheduler.

---

## Prerequisites
- Linux server with sufficient CPU/GPU, RAM, and disk.
- `git`, `wget`/`curl`, and a compiler toolchain (for some Python packages).
- Network access for model downloads and ColabFold server usage.
- Access to PyRosetta (registration/licensing required).

---

## Recommended environment (Miniconda)
```bash
# Download & install Miniconda (follow the installer's prompts)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# restart shell or source conda
source ~/.bashrc
```

# Create & activate environment (adjust Python version if needed)
```bash
conda create -n cf-env python=3.9 -y
conda activate cf-env
```

# Install project Python deps
```bash
pip install -r requirements.txt
```

---

## PyRosetta

PyRosetta must be installed manually. Visit the official site and follow the download/installation instructions:

[https://www.pyrosetta.org/downloads](https://www.pyrosetta.org/downloads)

Typical example (replace with the real filename you download):

```bash
pip install /path/to/PyRosetta-<version>-cp39-*.whl
```

> Do not commit PyRosetta binaries or license files to this repository.

---

## Local ColabFold

Install Local ColabFold per its repository instructions:

[https://github.com/YoshitakaMo/localcolabfold](https://github.com/YoshitakaMo/localcolabfold)

Follow that project's README to clone, install dependencies, and download required models.

---

## How to run

1. `cd` into the folder for the version of the code you want to run.
2. Run the Python script inside that folder.

Example:

```bash
cd versions/version_x
python run_version_x.py
```

### Long runs — run in background

Each run usually takes several days. Recommended approaches:

**nohup**

```bash
nohup python run_version_x.py > run_version_x.log 2>&1 &
# show PID of last background job (bash)
echo $!
# follow logs
tail -f run_version_x.log
```

---

## ColabFold server note

The current ColabFold code in this repository uses the **default server**. Be mindful of server rules, quotas, and fair usage. If your workload is heavy or frequent, consider:

* Hosting your own ColabFold server, or
* Using a dedicated compute cluster or cloud resources.

---

## Recommended .gitignore (suggested)

Add to your repo `.gitignore` to avoid committing large files, logs, or envs:

```
# Python
__pycache__/
*.pyc

# Logs
*.log
nohup.out

# PyRosetta, ColabFold, Gridsearch and Initial Screening results and cache
*.pdb
*.fasta
Results/
Temp/
```

---

## Troubleshooting (common issues)

* `pip install -r requirements.txt` fails: try `conda install -c conda-forge <package>` for heavy dependencies.
* PyRosetta import errors: confirm the wheel matches the Python version and OS; ensure you're in the conda env you installed into.
* Missing ColabFold models: follow Local ColabFold model download steps.
* Jobs killed unexpectedly: check system limits (`ulimit`), cluster scheduler policies, and system logs; prefer using a job scheduler for multi-day jobs.

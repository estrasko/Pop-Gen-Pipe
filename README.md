# Population Genomics Pipeline

Scripts useful for general population genomics workflows - Work in progress

## Project goals
In this project, we will write a script that streamlines our commonly used population genomics analyses into a clean 
pipeline. We know how to run the individual analyses. Our goal is to learn how to connect them together in a way 
that allows us to work faster using Python commands throughout to connect and guide the pipeline.

Goals:
- reproducible population genomics assembly and analysis
- use a Conda environment
- collaborate using git to create the pipeline
- test the pipeline using ~3 pop gen datasets

### End Goal: Produce quality figures and analyses results reproducable for publication

## Analysis using multiple SNPs per locus
- AMOVA
- DAPC
- IBD
- divmigrate

## Overview
This pipeline uses Python as the master program handling orchestration of runs, validation, and file handling. Every analysis
is able to be individually run as a function. The entire script is importable as a module.

All analyses are in the R programming language, where Python calls R for running the analyses and plotting the results.

## Workflow:

Input data
   ↓
Sample validation (Python)
   ↓
R analyses:
   ├── AMOVA
   ├── DAPC (with cross-validation)
   ├── IBD (Mantel + MRM)
   └── divMigrate (gene flow)
   ↓
Structured output directories + plots

## Installation
Clone the repository:

git clone https://github.com/estrasko/Pop-Gen-Pipe.git
cd Pop-Gen-Pipe

NOTE: you can also fork the repository and work that way.

## Create the Conda Environment
1) Put "environment.yml" in your designated working folder
2) Create conda environment under a new name: *conda env create -f environment.yml --name Pop-Gen-PipeEnv*
3) Activate your new environment: *conda activate Pop-Gen-PipeEnv*

## Required Input Files
*NOTE: ALL INPUT FILES MUST BE IN SAME ORDER (by population)!* critically, FST and geo matrices must have identical dimensions
and identical population order.

### Genepop files
| File           | Purpose           |
| -------------- | ----------------- |
| `haps.genepop` | AMOVA             |
| `snps.genepop` | DAPC + divMigrate |

### Popmap (popmap.csv)
CSV file with:
Sample,Population
LT-pop_01,Buxahatchee
LT-pop_02,Buxahatchee
...

the pattern is name of individual followed by population of origin. comma separated values.

### FST matrix (fst.csv)
square matrix:
0,0.24,0.23,0.18
0.24,0,0.19,0.13
0.23,0.19,0,0.13
0.18,0.13,0.13,0

### Geographic distance matrix (geo.csv)
square matrix:
0,170.41,138.18,80.14
170.41,0,77.68,90.25
138.18,77.68,0,57.79
80.14,90.25,57.79,0

## Run the Pipeline
python Pop_script_2.py \
  --haps-genepop populations.haps.genepop \
  --multi-snp-genepop populations.snps.genepop \
  --popmap popmap.csv \
  --fst-csv fst.csv \
  --geo-csv geo.csv \
  --outdir results \
  --scripts-dir . \
  --run-amova \
  --run-dapc \
  --run-ibd \
  --run-divmigrate

any of the above analyses can be removed or run individually, like so:

python Pop_script_2.py \
  --haps-genepop populations.haps.genepop \
  --multi-snp-genepop populations.snps.genepop \
  --popmap popmap.csv \
  --fst-csv fst.csv \
  --geo-csv geo.csv \
  --outdir results \
  --scripts-dir . \
  --run-amova \

## Optional: Multithreading for divmigrate
We created this pipeline to run on personal laptops, clusters, or whatever you have to work with. The only sometimes
computationally expensive program is divmigrate. If no threading option is specified, the default behavior is threads = 1.

Users can utilize more CPU resources by optionally flagging --threads <N>

Example: 
--run-divmigrate --threads 12

Note: threading is only available in divmigrate, not the other tests. The others don't need it.

# Feedback

You have a good start here.
I would recommend adding PCA to your list of analyses; this can be done using
all SNPs or one SNP per locus.
Creating the whole pipeline, from raw reads to results with figures might be a
bit ambitious for your class project.
I recommend you start by focusing on one part of the overall pipeline first.
For example, you could focus on starting with the assembly (or assemblies)
output by STACKS or other assemblers, and automating several of the analyses.

If you use Python, you will likely want to use the subprocess module for
running other tools outside of Python.
[Here is an intro to Python's subprocess module](https://www.geeksforgeeks.org/python/python-subprocess-module/).

You can use other languages for the project too (Python is not required).

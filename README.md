# Population Genomics Pipeline

Scripts useful for general population genomics workflows - Work in progress

## Project goals
In this project, we will write a script that streamlines our commonly used population genomics analyses into a clean 
pipeline. We know how to run the individual analyses. Our goal is to learn how to connect them together in a way 
that allows us to work faster using Python commands throughout to connect and guide the pipeline.

Goals:
- reproducible population genomics assembly and analysis
- use a container 
- collaborate using git to create the pipeline
- test the pipeline using ~3 pop gen datasets
### End Goal: Produce quality figures and analyses results reproducable for publication

### Assembly
- STACKS
- populations

### Analysis using multiple SNPs per locus
- Population summary statistics (from STACKS)
- AMOVA
- DAPC
- IBD
- divmigrate
- fineradstructure

### Analysis using single SNP per locus
- AdmixPipe

### Analysis using haplotypes from multiple SNPs per locus dataset
- allelic richness

## Conda Environment Installation
1) Put "environment.yml" in your designated working folder
2) Create conda environment under a new name: *conda env create -f environment.yml --name Pop-Gen-PipeEnv*
3) Activate your new environment: *conda activate Pop-Gen-PipeEnv*


## Pop_Script outline
- AMOVA
    - input file
        - columns
        - rows
    - functions
    - genclone object
    - AMOVA results
    - Randomization Test Results
- DAPC
- IBD


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

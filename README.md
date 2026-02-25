# Population Genomics Pipeline

Scripts useful for general population genomics workflows - Work in progress

## Project goals
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

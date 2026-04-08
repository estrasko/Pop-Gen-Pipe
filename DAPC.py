### DAPC ###

setwd("C:/R Folder")

library(adegenet)     # core DAPC tools
library(vcfR)         # quick VCF import
library(tidyverse)    # pipes + data wrangling
library(cowplot)      # to replicate the little inset barplot later


library(adegenet)
library(ape)
rm(list=ls())

data<-read.genepop(file="populations.snps.gen")

find.clusters(data,max.n.clust=25)

groups2<-find.clusters(data,max.n.clust=25)
# 200, 2
groups2
help(find.clusters)

dapc2<-dapc(data,groups2$grp)
# 4, 3
# 3 is number of eigenvalues


pdf("Leptoxis_coosaensis_DAPC.pdf")
scatter(dapc2,scree.da=TRUE,posi.da = "topleft", cell=0, cstar=0, clab=0, col=c("#FF6600","#0099E5"))  #can use grp=data$pop to color by collection locality
#scatter(dapc2,scree.da=TRUE,posi.da = "topleft", cell=0, cstar=0, clab=0, grp=data$pop)


## assumes you already have dapc2 from previous steps
scatter(dapc2,
        scree.da = TRUE,           # tiny bar-plot of DA eigenvalues
        posi.da  = "bottomleft",      # where that bar-plot goes
        cell     = 0,              # no inertia ellipses
        cstar    = 0,              # no star segments
        clab     = 0,              # hide point labels
        col      = c("#FF6600", "#0099E5", "#7C93C8", "#5FB58A")[1:nlevels(dapc2$grp)])


print(scatter)

dev.off()

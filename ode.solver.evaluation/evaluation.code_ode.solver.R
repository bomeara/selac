## header taken from selac.R
###LOAD REQUIRED PACKAGES -- eventually move to namespace:
library(ape)
library(expm)
library(nnet)
library(nloptr)
library(seqinr)
library(phangorn)
library(MASS)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(inline)
##need to load deSolve before dyn.load(*.so)
## or else you get errors like
## Error in checkDLL(func, jacfunc, dllname, initfunc, verbose, nout, outnames) : 
##  'initfunc' not loaded  initmod_selacHMM
library(deSolve)

## don't need to load selac package
## load compiled c code
## May need to recompile if source file has changed.
## To do so go to src/ and enter
##      R CMD SHLIB selacHMM.c
## only need to do this once
## We're doing this because we're not
## loading the selac package
dyn.load("../src/selacHMM.so") 


##sourcing selac.R
source("../R/selac.R")



##Load dataset
tree <- read.tree("rokasYeast.tre") ##in directory testthat/ which should be the working dir
phy <- drop.tip(tree, "Calb")
##reset branch length of phy$edge.length[1]=12 for testing long branch lengths.
yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
chars <- DNAbinToCodonNumeric(yeast.gene)
codon.data <- chars[phy$tip.label,]
aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=1)
aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
aa.optim.full.list <- aa.optim
codon.freq.by.aa <- GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=1)
codon.freq.by.gene <- GetCodonFreqsByGeneHMM(codon.data[,-1])
aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
colnames(aa.optim.frame.to.add) <- colnames(codon.data)
codon.data <- rbind(codon.data, aa.optim.frame.to.add)
codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
aa.optim = codon.data$optimal.aa
codon.index.matrix <- CreateCodonMutationMatrixIndexEvolveAA()


##source my modified files
source("./treeTraversalODETests.R")


runDiagnostics=TRUE; ## for ode solver output

###Calc likelihood using standard set up

## methods from ode help page
if(FALSE){
    odeMethodVec = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
                "euler", "rk4", "ode23", "ode45", "radau", 
                "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "iteration")
    }else{
        odeMethodVec = c("lsoda")
    }

for(odeMethod in odeMethodVec){

    ##WARNING: RUNNIGN THIS FROM WITHIN EMACS SEEMS TO NOT ALWAYS WORK PROPERLY
    ##TRY CUTTING AND PASTING OR SOURCING
    
    runTime = system.time(
        selac.unrest.evolveAA <- GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, rep(1,11), 0.01)), codon.data, phy, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="UNREST", codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=1) #n.cores.by.gene?
    )
    roundLLik <- round(selac.unrest.evolveAA, 3);

    print(paste("Run time: ", runTime, " LLik for ", odeMethod, " with settings" , ": ", roundLLik, " \t deviation:", -8677.442 - roundLLik , sep=""))
    comparison <- identical(round(selac.unrest.evolveAA, 3), -8677.442) ##Correct LLik value seen using lsoda and matrix exponentiation
    }

 


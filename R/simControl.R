
source("selacSim.R")

##

selacSimControl <- function(start.rep=1, end.rep=100, nuc.model, numcode=1, alpha=0.4, beta=0.1, Ne=5e6, numgenes=5, yeast.fit.output, aa.properties=NULL){
    phy = yeast.fit.output$phy
	rows.to.sample = sample(1:106, numgenes)
	write.table(rows.to.sample, file="GenesSampledForSim", quote=FALSE, sep="\t", row.names=FALSE)
    for(nrep.index in start.rep:end.rep){
        system(paste("mkdir", paste("fastaSet",nrep.index, sep="_"), sep=" "))
        for(gene.no.index in 1:numgenes){
            empirical.codon.freqs = yeast.fit.output$empirical.codon.freqs[[rows.to.sample[gene.no.index]]]
            par.vec = c(yeast.fit.output$mle.pars[rows.to.sample[gene.no.index],1], alpha, beta, Ne, yeast.fit.output$mle.pars[rows.to.sample[gene.no.index],5:12])
            optimal.aa = yeast.fit.output$aa.optim[[rows.to.sample[gene.no.index]]]
            tmp <- SelacSimulator(phy=phy, pars=par.vec, aa.optim_array=optimal.aa, root.codon.frequencies=empirical.codon.freqs, numcode=1, aa.properties=NULL, nuc.model=nuc.model, k.levels=0, diploid=FALSE)
			write.dna(tmp, file=paste(paste("fastaSet",nrep.index, sep="_"), "/gene",  gene.no.index, ".fasta", sep=""), format="fasta", colw=1000000)
        }
#      selac.fit <- EstimateParametersCodonGlobal(paste(paste("fastaSet", nrep.index, sep="_"),"/",sep=""), n.partitions=numgenes, nuc.model=nuc.model, diploid=FALSE, n.cores=numgenes)
#      save(selac.fit, file=paste("selac", nrep.index, sep="_"))
    }
}

load("yeast.finalSELAC106.Rdata")
print(ls())
selacSimControl(start.rep=1, end.rep=25, nuc.model="GTR", numgenes=5, yeast.fit.output=result)





source("selacSim.R")

##

selacSimControl <- function(start.rep=1, end.rep=100, nuc.model, numcode=1, alpha=0.4, beta=0.1, Ne=5e6, numgenes=5, yeast.fit.output, aa.properties=NULL){
    phy = yeast.fit.output$phy
#	rows.to.sample = sample(1:106, numgenes)
#	write.table(rows.to.sample, file="GenesSampledForSim", quote=FALSE, sep="\t", row.names=FALSE)
#	rows.to.sample <- c(66, 62, 93, 96,  3,  26,  81, 69, 105, 20)
	rows.to.sample <- c(106, 91, 81, 8, 24, 23, 53, 76, 12, 102)
    for(nrep.index in start.rep:end.rep){
        system(paste("mkdir", paste("fastaSet",nrep.index, sep="_"), sep=" "))
        for(gene.no.index in 1:numgenes){
            empirical.codon.freqs = yeast.fit.output$empirical.codon.freqs[[rows.to.sample[gene.no.index]]]
            par.vec = c(yeast.fit.output$mle.pars[rows.to.sample[gene.no.index],1], alpha, beta, Ne, yeast.fit.output$mle.pars[rows.to.sample[gene.no.index],5:12])
            optimal.aa = yeast.fit.output$aa.optim[[rows.to.sample[gene.no.index]]]
            tmp <- SelacSimulator(phy=phy, pars=par.vec, aa.optim_array=optimal.aa, root.codon.frequencies=empirical.codon.freqs, numcode=1, aa.properties=NULL, nuc.model=nuc.model, k.levels=0, diploid=FALSE)
			write.dna(tmp, file=paste(paste("fastaSet",nrep.index, sep="_"), "/gene",  gene.no.index, ".fasta", sep=""), format="fasta", colw=1000000)
        }
		selac.fit <- EstimateParametersCodonGlobal(paste(paste("fastaSet", nrep.index, sep="_"),"/",sep=""), n.partitions=numgenes, phy=phy, nuc.model=nuc.model, diploid=FALSE, n.cores=numgenes)
		save(selac.fit, file=paste("selac_", nrep.index, ".Rsave", sep=""))
    }
}

load("yeast.finalSELAC106.Rdata")
selacSimControl(start.rep=1, end.rep=50, nuc.model="GTR", numgenes=10, yeast.fit.output=result)

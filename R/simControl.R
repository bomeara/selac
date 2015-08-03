
source("selacSimulation.R")

selacSimControl <- function(nreps=50, nuc.model, numcode=1, sensitivity=2, alpha=0.3876658, beta=0.08483003, Ne=5e6, gene.no=5, yeast.fit.output, aa.properties=NULL, n.cores){
    phy = result$phy

    #Backtransform to Cphiqs#
    C = 2
    Phi <- 0.5
    q <- 4e-7
    C.Phi.q.s <- C * Phi * q * sensitivity
    ##########################
    
    for(nrep.index in 1:nreps){
        system(paste("mkdir", paste("fastaSet",nrep.index, sep="_"), sep=" "))
        rows.to.sample = sample(gene.no, 1:106)
        for(gene.no.index in 1:gene.no){
            empirical.codon.freqs = result$empirical.codon.freqs[[rows.to.sample[gene.no.index]]]
            par.vec = c(C.Phi.q.s, alpha, beta, Ne, result$mle.pars[rows.to.sample[gene.no.index],5:12])
            optimal.aa = result$aa.optimal[[rows.to.sample[gene.no.index]]]
            tmp <- SelacSimulator(nsites=result$nsites[rows.to.sample[gene.no.index]], phy=phy, pars=par.vec, aa.optim_array=aa.optim, empirical.codon.freqs=empirical.codon.freqs, numcode=1, aa.properties=NULL, nuc.model=nuc.model)
        }
        selac.fit <- EstimateParametersCodonGlobal(paste(paste("fastaSet", nrep.index, sep="_"),"/",sep=""), n.partitions=gene.no, nuc.model=nuc.model, n.cores=n.cores)
        save(selac.fit, file=paste("selac", nrep.index, sep="_"))
    }
}





###Modified version of function taken 4/25/17 from selac.R in git repository
## Changing to explore different ODE solvers in hopes of speeding things up.

TreeTraversalODE <- function(phy, Q_codon_array_vectored, liks.HMM, bad.likelihood=-100000, root.p) {
    if(runDiagnostics) print("site: Unknown")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode

    anc <- unique(phy$edge[,1])
    TIPS <- 1:nb.tip
    comp <- numeric(nb.tip + nb.node)
    #Start the postorder traversal indexing lists by node number:
    for (i in seq(from = 1, length.out = nb.node)) {
        focal <- anc[i]
        desRows <- which(phy$edge[,1]==focal)
        desNodes <- phy$edge[desRows,2]
        v = rep(1, dim(liks.HMM)[2])
        if(runDiagnostics) print(paste("Node: ", i, " out of ", nb.node, "\n")) 
        for (desIndex in sequence(length(desRows))){
            if(runDiagnostics) print(paste("descendent: ", desIndex, " out of ", length(desRows), "\n")) 
            yini <- liks.HMM[desNodes[desIndex],]
            times=c(0, phy$edge.length[desRows[desIndex]])
	    ##Change ODE solver here
            ## swapped lsoda() for ode() so we can define solver as an argument to function.
            ## Note to access all of the rk functions you can use
            ##         method = rkMethod("rk45dp7", densetype = NULL, nknots = 5)
            ##
            prob.subtree.cal.full <- ode(y=yini, times=times, func = "selacHMM", parms=Q_codon_array_vectored, initfunc="initmod_selacHMM", dllname = "selacHMM", method=odeMethod, hini=0)
            ## options for runge-kutta rk()
            ##rk(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
            ##verbose = FALSE, tcrit = NULL, hmin = 0, hmax = NULL,
            ##hini = hmax, ynames = TRUE, method = rkMethod("rk45dp7", ... ),
            ##maxsteps = 5000, dllname = NULL, initfunc = dllname,
            ##initpar = parms, rpar = NULL, ipar = NULL,
            ##nout = 0, outnames = NULL, forcings = NULL,
            ##initforc = NULL, fcontrol = NULL, events = NULL, ...)
            ##
            ## method : the integrator to use. This can either be a string constant
            ##naming one of the pre-defined methods or a call to function
            ##‘rkMethod’ specifying a user-defined method.  The most common
            ##methods are the fixed-step methods ‘"euler"’, second and
            ##fourth-order Runge Kutta (‘"rk2"’, ‘"rk4"’), or the variable
            ##step methods Bogacki-Shampine ‘"rk23bs"’,
            ##Runge-Kutta-Fehlberg ‘"rk34f"’, the fifth-order Cash-Karp
            ##method ‘"rk45ck"’ or the fifth-order Dormand-Prince method
            ##with seven stages ‘"rk45dp7"’.  As a suggestion, one may use
            ##‘"rk23bs"’ (alias ‘"ode23"’) for simple problems and
            ##‘"rk45dp7"’ (alias ‘"ode45"’) for rough problems.

            if(runDiagnostics){
                print("Running diagnostics\n")
                diagnostics(prob.subtree.cal.full)
                
            } 
            ######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
            if(attributes(prob.subtree.cal.full)$istate[1] < 0){
                return(bad.likelihood)
            }else{
                prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
            }
            ##############################################################################

            if(prob.subtree.cal[1]<0){
                return(bad.likelihood)
            }
            v <- v * prob.subtree.cal
        }
        comp[focal] <- sum(v)
        liks.HMM[focal,] <- v/comp[focal]
    }
    root.node <- nb.tip + 1L
    if (is.na(sum(log(liks.HMM[root.node,])))){
        return(bad.likelihood)
    }else{
        loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks.HMM[root.node,])))
    }
    return(loglik)
}

# HAS TO BE LOADED IN THIS ORDER OR IT MIGHT NOT WORK!!!
library(deSolve)
dyn.load("~/SELAC/selac/src/selacHMM.so")
load("~/SELAC/c_testing/Q_codon_vec.Rda")

n.times <- 20
run.times <- vector("numeric", n.times)

# Leaf initial conditions
#yini <- rep(0, 1344)
#yini[1] <- 1 

# internal node initial conditions
yini <- runif(1344)
yini <- yini / sum(yini)

for(t in 1:n.times)
{
  times <- c(0,t)
  
  start <- Sys.time()
  res <- euler(yini, times, func = "selacHMM", Q_codon_array_vectored, initfunc="initmod_selacHMM", dllname = "selacHMM")
  end <- Sys.time()
  
  cat(".")  
  run.times[t] <- end - start
}

pdf("~/SELAC/c_testing/ode_time_rand.pdf")
plot(cumsum(run.times), type="l", ylab = "Runtime [sec]", xlab = "ODE Time Units", main = "Random Initial Conditions")
lines(run.times, col="red")
legend("topleft", col=c("black", "red"), lty = 1, lwd = 2, legend = c("Cum Time", "Indiv Time"), bty="n")
abline(h=1, lty=2)
dev.off()

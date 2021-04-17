# Recreating figures 2.2 and 2.3

library(tidyverse)

###' @param N number of banks
###' @param a mean-reversion rate
###' @param sigma diffusion
###' @param X0 initial value for process X, assumed iid
###' @param end_time end of sim
###' @param dt Euler time step
simple_lending_model <- function(N,a,sigma,D=-0.7,rho=0,X0=0,end_time=1,dt=1e-4){
  ## set-up
  Nsteps <- end_time / dt + 1   # number of steps to take
  X <- matrix(nrow=N,ncol=Nsteps)   # empty matrix to store results
  ensemble_avg <- rep(NA,Nsteps)
  X[,1] <- rep(X0,N)   # X_0 = X0 for all banks
  ensemble_avg[1] <- X0
  defaults <- rep(0,N)
  
  for (n in 2:Nsteps){
    dW0t <- rnorm(1)   # common noise
    dWit <- rnorm(N)   # independent Brownian motion
    dWt <- rho * rep(dW0t,N) + sqrt(1-rho^2)*dWit   # total Brownian motion term
    X_bar <- sum(X[,n-1]) / N   # mean
    # Euler step:
    X[,n] = X[,n-1] + a * (rep(X_bar,N) - X[,n-1]) * dt + sigma * sqrt(dt) *  dWt
    ensemble_avg[n] <- ensemble_avg[n-1] + sqrt(dt) * sigma / N * sum(dWt)
    # record defaults
    defaults <- defaults + ifelse(X[,n]<= D,1,0)
  }
  # clean up defaults (1 if defaulted at least once)
  defaults <- ifelse(defaults>0,1,0)
  # clean up and return tibble of banks trajectories and times
  rownames(X) <- sapply(1:N, function(x) paste0("bank",x))
  times <- seq(0,end_time,by=dt)
  return(list(as_tibble(cbind(t(X),times,ensemble_avg)), defaults))
}



###' @param Nsims number of simulations for Monte Carlo
###' all other params as above
default_distribution <- function(Nsims,N,a,sigma=1,D=-0.7,rho=0,X0=0){
  # output is saved since it takes a long time to compute this
  # name file to save data:
  savefile <- sprintf("saved_sims/Nsims_%g_Nbanks_%g_a_%g_sigma_%g_rho_%g.Rdata",
                      Nsims,N,a,sigma,rho)
  if (!file.exists(savefile)){
    defaults <- rep(NA,Nsims)
    for (sim in 1:Nsims){
      # run the simple model Nsims times
      X <- simple_lending_model(N=N,a=a,sigma=sigma,D=D,rho=rho,X0=X0)
      # count number of defaults this run
      defaults[sim] <- sum(X[[2]])
    }
    save(defaults,file=savefile)
  } else{
    load(savefile)
  }
  return(defaults)
}

defaults1 <- default_distribution(Nsims=10000,N=10,a=0)
defaults2 <- default_distribution(Nsims=10000,N=10,a=10)
defaults3 <- default_distribution(Nsims=10000,N=10,a=100)

# figure 2.2 from paper
a1=as_tibble(defaults1) %>% count(value) %>% 
  mutate(n=n/10000) %>% rbind(c(value=10,n=0))
a2=as_tibble(defaults2) %>% count(value) %>% mutate(n=n/10000) 

fig_2_2a <- ggplot(data=a1,aes(value,n)) + 
  geom_line(colour='blue',lty=2) + 
  geom_line(data=a2, colour='blue') + 
  theme_classic() +
  xlab("# of default") + ylab("prob of # of default")

fig_2_2b <- ggplot(data=a1,aes(value,n)) + 
  geom_line(colour='blue',lty=2) + 
  geom_line(data=a2, colour='blue') + 
  theme_classic() + 
  coord_cartesian(xlim = c(6, 10),ylim=c(0,0.2)) +
  xlab("# of default") + ylab("prob of # of default")

gridExtra::grid.arrange(fig_2_2a,fig_2_2b,nrow=1)

# figure 2.2 with KDE instead
defaults12 <- rbind(mutate(as_tibble(defaults1),a = 0),
                    mutate(as_tibble(defaults2),a = 10)) 

plot_a <- defaults12 %>% mutate(a = as.factor(a)) %>%
  ggplot(aes(value, group=a,fill=a)) + 
  geom_density(adjust=2.5, alpha=0.5) + 
  theme_minimal() +
  xlab("# of default")

plot_b <- defaults12 %>% mutate(a = as.factor(a)) %>%
  ggplot(aes(value, group=a, fill=a)) + geom_density(adjust=2.5, alpha=0.5) + 
  theme_minimal() + coord_cartesian(xlim = c(6, 10),ylim=c(0,0.2)) +
  xlab("# of default")

gridExtra::grid.arrange(plot_a,plot_b,nrow=1)

#---
# figure 2.3 from paper
a3=as_tibble(defaults3) %>% count(value) %>% mutate(n=n/10000) 

fig_2_3a <- ggplot(data=a1,aes(value,n)) + 
  geom_line(colour='blue',lty=2) + 
  geom_line(data=a3, colour='blue') + 
  theme_classic() +
  xlab("# of default") + ylab("prob of # of default")

fig_2_3b <- ggplot(data=a1,aes(value,n)) + 
  geom_line(colour='blue',lty=2) + 
  geom_line(data=a3, colour='blue') + 
  theme_classic() + coord_cartesian(xlim = c(6, 10),ylim=c(0,0.2)) +
  xlab("# of default") + ylab("prob of # of default")

gridExtra::grid.arrange(fig_2_3a,fig_2_3b,nrow=1)

# figure 2.3 with KDE instead
defaults13 <- rbind(mutate(as_tibble(defaults1),a = 0),
                    mutate(as_tibble(defaults3),a = 100)) 

plot_c <- defaults13 %>% mutate(a = as.factor(a)) %>%
  ggplot(aes(value, group=a,fill=a)) + 
  geom_density(adjust=2.5, alpha=0.5) + 
  theme_minimal() +
  xlab("# of default")

plot_d <- defaults13 %>% mutate(a = as.factor(a)) %>%
  ggplot(aes(value, group=a, fill=a)) + geom_density(adjust=2.5, alpha=0.5) + 
  theme_minimal() + coord_cartesian(xlim = c(6, 10),ylim=c(0,0.2)) +
  xlab("# of default")

gridExtra::grid.arrange(plot_c,plot_d,nrow=1)



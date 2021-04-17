# making plots for use in the presentation and paper
library(tidyverse)
set.seed(720)

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

res1 = simple_lending_model(N=10,a=10,sigma=1,rho=0,X0=0)
plot1 <- res1[[1]] %>% 
  pivot_longer(cols = 'bank1':'bank10',values_to = 'value',names_to = 'bank') %>% 
  ggplot(aes(times,value,colour=bank)) + geom_line(alpha=0.5) +
  #geom_line(aes(times,ensemble_avg),colour='black') +
  geom_hline(yintercept = -0.7) + xlab('') + ylab('') +
  theme_bw(base_size = 18) + theme(legend.position = "none") + ylim(-2,2)
plot1
#ggsave('plot1a.png',height=5,width=7)

res2 = simple_lending_model(N=10,a=0,sigma=1,rho=0,X0=0)
plot2 <- res2[[1]] %>% 
  pivot_longer(cols = 'bank1':'bank10',values_to = 'value',names_to = 'bank') %>% 
  ggplot(aes(times,value,colour=bank))+ geom_line(alpha=0.5) +
  #geom_line(aes(times,ensemble_avg),colour='black') +
  geom_hline(yintercept = -0.7) + xlab('') + ylab('') +
  theme_bw(base_size = 18) + theme(legend.position = "none") + ylim(-2,2)
plot2
#ggsave('plot1b.png',height=5,width=7)

#-------
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


# figure 2.2 with KDE instead
defaults12 <- rbind(mutate(as_tibble(defaults1),a = 0),
                    mutate(as_tibble(defaults2),a = 10)) 

plot_a <- defaults12 %>% mutate(a = as.factor(a)) %>%
  ggplot(aes(value, group=a,fill=a)) + 
  geom_density(adjust=2.5, alpha=0.5) + 
  theme_minimal(base_size = 18) +
  xlab("Number of defaults") + ylab("Density")
plot_a
#ggsave('plot2a.png',height=5,width=7)

#---

# figure 2.3 with KDE instead
defaults13 <- rbind(mutate(as_tibble(defaults1),a = 0),
                    mutate(as_tibble(defaults3),a = 100)) 

plot_c <- defaults13 %>% mutate(a = as.factor(a)) %>%
  ggplot(aes(value, group=a,fill=a)) + 
  geom_density(adjust=2.5, alpha=0.5) + 
  theme_minimal(base_size = 18) +
  xlab("Number of defaults") + ylab("Density")
plot_c
#ggsave('plot2b.png',height=5,width=7)

#--------------
plot3a <- res1[[1]] %>% 
  pivot_longer(cols = 'bank1':'bank10',values_to = 'value',names_to = 'bank') %>% 
  ggplot(aes(times,value,colour=bank)) + geom_line(alpha=0.5) +
  geom_line(aes(times,ensemble_avg),colour='black') +
  geom_hline(yintercept = -0.7) + xlab('') + ylab('') +
  theme_bw(base_size = 18) + theme(legend.position = "none") + ylim(-1,0.8)
plot3a
#ggsave('plot3a.png',height=5,width=7)

res4 = simple_lending_model(N=10,a=100,sigma=1,rho=0,X0=0)
plot3b <- res4[[1]] %>% 
  pivot_longer(cols = 'bank1':'bank10',values_to = 'value',names_to = 'bank') %>% 
  ggplot(aes(times,value,colour=bank))+ geom_line(alpha=0.5) +
  geom_line(aes(times,ensemble_avg),colour='black') +
  geom_hline(yintercept = -0.7) + xlab('') + ylab('') +
  theme_bw(base_size = 18) + theme(legend.position = "none") + ylim(-1,0.8)
plot3b
#ggsave('plot3b.png',height=5,width=7)


#---------------
# Figure 2.4
res5 = simple_lending_model(N=10,a=10,sigma=1,rho=0.25,X0=0)
plot1 <- res5[[1]] %>% 
  pivot_longer(cols = 'bank1':'bank10',values_to = 'value',names_to = 'bank') %>% 
  ggplot(aes(times,value,colour=bank)) + geom_line(alpha=0.5) +
  geom_line(aes(times,ensemble_avg),colour='black') +
  geom_hline(yintercept = -0.7) + xlab('') + ylab('') +
  theme_bw(base_size = 18) + theme(legend.position = "none") + ylim(-2,2)
plot1
#ggsave('plot4a.png',height=5,width=7)

res6 = simple_lending_model(N=10,a=10,sigma=1,rho=0.75,X0=0)
plot2 <- res6[[1]] %>% 
  pivot_longer(cols = 'bank1':'bank10',values_to = 'value',names_to = 'bank') %>% 
  ggplot(aes(times,value,colour=bank))+ geom_line(alpha=0.5) +
  geom_line(aes(times,ensemble_avg),colour='black') +
  geom_hline(yintercept = -0.7) + xlab('') + ylab('') +
  theme_bw(base_size = 18) + theme(legend.position = "none") + ylim(-2,2)
plot2
#ggsave('plot4b.png',height=5,width=7)


#----------------
defaults4 <- default_distribution(Nsims=10000,N=10,a=10,rho=0.25)
defaults5 <- default_distribution(Nsims=10000,N=10,a=10,rho=0.5)
defaults6 <- default_distribution(Nsims=10000,N=10,a=10,rho=0.75)

# defaults2 is rho=0

defaults_combo <- rbind(mutate(as_tibble(defaults2),rho = 0),
                        mutate(as_tibble(defaults4),rho = 0.25),
                    mutate(as_tibble(defaults5),rho = 0.5),
                    mutate(as_tibble(defaults6),rho = 0.75)) %>% 
  mutate(rho = as.factor(rho))

plot <- ggplot(defaults_combo,aes(value, group=rho,fill=rho)) + 
  geom_density(adjust=2.5, alpha=0.5) + theme_minimal(base_size = 14) +
  xlab("Number of defaults") + ylab("Density")
plot
#ggsave('plot5.png',height=4,width=8)

# Comparing the open and closed loop solutions
eta <- function(t,c,N,a,q,epsilon,end_time){
  R <- (a+q)^2 + (1 - 1/(N^2))*(epsilon-q^2)
  delta_plus <- -(a+q) + sqrt(R)
  delta_minus <- -(a+q) - sqrt(R)
  exponent_term <- exp((delta_plus - delta_minus)*(end_time - t))
  eta_num <- -(epsilon-q^2) * (exponent_term - 1) - c *
    (delta_plus * exponent_term - delta_minus)
  eta_denom <- delta_minus * exponent_term - delta_plus - c *
    (1 - 1/(N^2)) * (exponent_term - 1)
  eta <- eta_num / eta_denom
  return(eta)
}

phi <- function(t,c,N,a,q,epsilon,end_time){
  R <- (a+q)^2 + (1 - 1/N)*(epsilon-q^2)
  delta_plus <- -(a+q) + sqrt(R)
  delta_minus <- -(a+q) - sqrt(R)
  exponent_term <- exp((delta_plus - delta_minus)*(end_time - t))
  phi_num <- -(epsilon-q^2) * (exponent_term - 1) - c *
    (delta_plus * exponent_term - delta_minus)
  phi_denom <- delta_minus * exponent_term - delta_plus - c *
    (1 - 1/N) * (exponent_term - 1)
  phi<- phi_num / phi_denom
  return(phi)
}

tseq <- seq(0,1,.01)
eta1 <- eta(t=tseq,c=0,N=10,a=1,q=1,epsilon=10,end_time=1)
eta2 <- eta(t=tseq,c=1,N=10,a=1,q=1,epsilon=10,end_time=1)
phi1 <- phi(t=tseq,c=0,N=10,a=1,q=1,epsilon=10,end_time=1)
phi2 <- phi(t=tseq,c=1,N=10,a=1,q=1,epsilon=10,end_time=1)

as_tibble(cbind(tseq,eta1,eta2,phi1,phi2)) %>% 
  pivot_longer(eta1:phi2, names_to= 'name', values_to = 'value') %>% 
  mutate(variable = fct_rev(as.factor(ifelse(name %in% c('phi1','phi2'),'phi','eta')))) %>% 
  mutate(c = as.factor(ifelse(name %in% c('phi1','eta1'),0,1))) %>% 
  ggplot() + theme_bw(base_size = 14) + xlab('t') +
  ylab(expression(paste(phi[t],", ",eta[t]))) +
  geom_line(aes(tseq,value,colour=c,linetype=variable), size=0.8) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c(expression(phi[t]),expression(eta[t])))+
  theme(legend.key.width = unit(2, "line"))


#ggsave('plot6.png',height=4,width=8)

# limit of the eta and phi as N to infty
eta_limit <- function(t,c,a,q,epsilon,end_time){
  R <- (a+q)^2 + (epsilon-q^2)
  delta_plus <- -(a+q) + sqrt(R)
  delta_minus <- -(a+q) - sqrt(R)
  exponent_term <- exp((delta_plus - delta_minus)*(end_time - t))
  phi_num <- -(epsilon-q^2) * (exponent_term - 1) - c *
    (delta_plus * exponent_term - delta_minus)
  phi_denom <- delta_minus * exponent_term - delta_plus - c *(exponent_term - 1)
  phi<- phi_num / phi_denom
  return(phi)
}


# Comparing the effective interbank borrowing rates
a <- 1
q <- 1
Nseq <- seq(10,100,0.1)
open_rate <- a + q + (1 - 1/Nseq) * phi(t=0,c=0,N=Nseq,a=a,q=q,epsilon=10,end_time=1)
closed_rate <- a + q + (1 - 1/Nseq) * eta(t=0,c=0,N=Nseq,a=a,q=q,epsilon=10,end_time=1)
limit <- a + q + eta_limit(t=0,c=0,a=a,q=q,epsilon=10,end_time=1)
limit <- rep(limit,length(Nseq))

as_tibble(cbind(Nseq,open_rate,closed_rate,limit)) %>%
  pivot_longer(open_rate:limit, names_to= 'name', values_to = 'value') %>%
  mutate(name = ifelse(name == "open_rate","open loop",name)) %>% 
  mutate(name = ifelse(name == "closed_rate","closed loop",name)) %>% 
  mutate(equilibrium = fct_relevel(as.factor(name),"open loop")) %>% 
  ggplot() + theme_bw(base_size = 14) + 
  geom_line(aes(Nseq,value,colour=equilibrium), size=0.8) +
  xlab('N') + ylab("Effective lending rate at equilibrium") +
  scale_colour_brewer(palette = 'Set2') 

#ggsave('plot7.png',height=4,width=8)



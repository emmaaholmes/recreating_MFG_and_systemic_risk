# Comparing the open and closed loop solutions
library(tidyverse)

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
  ggplot() + theme_bw() + xlab('t') +
  ylab(expression(paste(phi[t],", ",eta[t]))) +
  geom_line(aes(tseq,value,colour=c,linetype=variable)) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c(expression(phi[t]),expression(eta[t])))+
  theme(legend.key.width = unit(2, "line"))


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
  ggplot() + theme_bw() + 
  geom_line(aes(Nseq,value,colour=equilibrium)) +
  xlab('N') + ylab("Effective lending rate at equilibrium") +
  scale_colour_brewer(palette = 'Set2') 



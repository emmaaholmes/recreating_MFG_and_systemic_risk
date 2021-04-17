# Recreating figure 2.1

library(tidyverse)
set.seed(720)

###' @param N number of banks
###' @param a mean-reversion rate
###' @param sigma diffusion
###' @param X0 initial value for process X, assumed iid
###' @param end_time end of sim
###' @param dt Euler time step
simple_lending_model <- function(N,a,sigma,X0,end_time=1,dt=1e-4){
  ## set-up
  Nsims <- end_time / dt + 1   # number of steps to take
  X <- matrix(nrow=N,ncol=Nsims)   # empty matrix to store results
  X[,1] <- rep(X0,N)   # X_0 = X0 for all banks
  
  for (n in 2:Nsims){
    dWt <- sqrt(dt) * rnorm(N)   # Brownian motion
    X_bar <- sum(X[,n-1]) / N   # mean
    # Euler step:
    X[,n] = X[,n-1] + a * (rep(X_bar,N) - X[,n-1]) * dt + sigma * dWt
  }
  # clean up and return tibble of banks trajectories and times
  rownames(X) <- sapply(1:10, function(x) paste0("bank",x))
  times <- seq(0,end_time,by=dt)
  return(as_tibble(cbind(t(X),times)))
}

plot_banks <- function(res){
  as_tibble(cbind(t(res),times)) %>% 
    pivot_longer(cols = 'bank1':'bank10', values_to = 'value', names_to = 'bank') %>% 
    ggplot(aes(times,value,colour=bank))+ geom_line(alpha=0.5) +
    geom_hline(yintercept = -0.7) +
    theme_bw() + theme(legend.position = "none") + ylim(-2,2)
}

res1 = simple_lending_model(N=10,a=10,sigma=1,X0=0)
plot1 <- res1 %>% 
  pivot_longer(cols = 'bank1':'bank10',values_to = 'value',names_to = 'bank') %>% 
  ggplot(aes(times,value,colour=bank)) + geom_line(alpha=0.5) +
  geom_hline(yintercept = -0.7) + xlab('') + ylab('') +
  theme_bw() + theme(legend.position = "none") + ylim(-2,2) +
  ggtitle("a=10")

res2 = simple_lending_model(N=10,a=0,sigma=1,X0=0)
plot2 <- res2 %>% 
  pivot_longer(cols = 'bank1':'bank10',values_to = 'value',names_to = 'bank') %>% 
  ggplot(aes(times,value,colour=bank))+ geom_line(alpha=0.5) +
  geom_hline(yintercept = -0.7) + xlab('') + ylab('') +
  theme_bw() + theme(legend.position = "none") + ylim(-2,2) +
  ggtitle("a=0")

gridExtra::grid.arrange(plot1,plot2,nrow=1)



#simulates in R
sim_ARgamma = function(nu, rho_r, N){
  R = numeric(N)
  R[1] = rgamma(n = 1, shape = 1/nu, scale = nu)
  
  for(i in 2:N){
    m = (rho_r/nu)/(1-rho_r)*R[i-1]
    Npois = rpois(1, lambda = m)
    W = c(0, rexp(n = Npois, rate = 1))
    K = sum(W)/m
    R[i] = rho_r*K*R[i-1] + (1-rho_r)*rgamma(n = 1, shape = 1/nu, scale = nu)
  }
  
  return(R)
}
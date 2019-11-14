#params y = log returns, ysq = squared log returns
ARGVG_pred = function(y, ysq){
  
  #estimation of nu using the method of moments
  nu_hat = (moments::kurtosis(y)-3)/3
  
  #estimation of rho_r using regression method
  
  #calculates ACF for ysq
  rho_ysq_hat = acf(ysq, lag.max = length(ysq), plot=F)$acf[-1]
  #cuts off lags at (before) first insignificant lag in acf(Y^2)
  req_lags = which(rho_ysq_hat < 2/sqrt(length(y)))[1]
  rho_ysq_hat = rho_ysq_hat[1:(req_lags-1)]
  
  #estimates rho_r(k) hats
  rho_rk_hat = rho_ysq_hat*(2+3*nu_hat)/nu_hat
  #calculate linear fit of the log of these
  linreg = lm(log(rho_rk_hat)~seq(1:length(rho_rk_hat)))
  rho_r_hat = exp(linreg$coefficients[2])
  
  #predicts R_t+1 using Ysq_t according to developed method

  #calculate necessary bessel functions
  K_pp1 = besselK(sqrt(2*ysq/nu_hat), 1/nu_hat+1/2)
  K_p = besselK(sqrt(2*ysq/nu_hat), 1/nu_hat-1/2)
  
  #interpolates R_ts
  R_hat_int = sqrt(ysq*nu_hat/2)*K_pp1/K_p
  
  #using the expectation of R_t|R_(t-1): (sets R_1 = 1)
  R_hat_pred = c(1,rho_r_hat*R_hat_int + (1-rho_r_hat))
  
  #finally: if the return happened to be exactly zero for one day, Inf is produced in the bessel functions
  #primitive way to fix this is to set the NA/Inf hat values to the 1 (the expected value of R)
  #from tests this happens to less in <0.001 of days so it doesn't change much
  R_hat_pred[which(is.na(R_hat_pred))] <- 1
  
  
  return(list(nu_hat = nu_hat,
              rho_r_hat = rho_r_hat,
              R_hat_pred = R_hat_pred))
}

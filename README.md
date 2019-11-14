Some supplementary code in R for the ARGVG Model
================


Autoregressive Gamma variance Gaussian mixture model
----------------------------------------------------

This model aims to model volatile financial return series using the product *Y*<sub>*t*</sub> of two independent time series, namely a Gaussian (often limited to white noise in a financial return series context) *X*<sub>*t*</sub> and a Gamma component *R*<sub>*t*</sub> that modulates the variance. This repository provides R code to simulate *Y*<sub>*t*</sub>, estimate parameters on real data, and a primitive prediction method. For more on the model and its origins, see [Johannesson et al. (2016)](https://www.semanticscholar.org/paper/AR(1)-time-series-with-autoregressive-gamma-for-Johannesson-Podg%C3%B3rski/e42e7beb40de022361be7ce82a9d4f013ea8307a). For more on the methods to estimate parameters and predict *Y*<sub>*t*</sub>, see [[2]](http://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=8995981&fileOId=8996215).

In addition to the prediction function (`ARGVG_pred.R`) and the simulation function (`sim_arGamma.R`), two simple `rshiny` apps are provided in `shinyapp_argvg.R` and `shinyapp_est.R` to visualize the impact of the model's parameters, and to visualize the parameter estimation process respectively.



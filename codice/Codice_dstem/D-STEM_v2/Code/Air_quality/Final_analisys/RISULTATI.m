%% Confrontiamo i vari modelli
clear all
clc

load("result_data_multivariato_meteo_NOX.mat")
load("result_data_multivariato_meteo_PM10.mat")

tot_t_stat = (mean(result_data_multivariato_meteo_PM10{1,9},2) + mean(X{1,9},2))/2;

%{
ORDINE NELLA STRUCT
beta_cv;
theta_z_cv;
v_z_cv;
sigma_eta_cv;
G_cv;
sigma_eps_cv;
diag_varcov_cv;
log_likelihood_cv;
t_stat;
%}

load("result_data_pm10_univariate_selected.mat")
result_data_pm10_univariate_selected{1,7};

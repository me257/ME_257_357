%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ME357 Spring 2017
%  Final Project
%  1. Turbine Stator-Rotor Analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

cp = 1005;
gamma = 1.4;
T = 1600;
eta_st = 0.9;
u_x = 300;
U = 200;
n = 3;
alpha_b = 70 / 180 * pi;
beta_c  = 60 / 180 * pi;




[p_ratio, eta_t, T_out, Work_out] = turbfn(cp, gamma, T, eta_st, u_x, U, n, alpha_b, beta_c)
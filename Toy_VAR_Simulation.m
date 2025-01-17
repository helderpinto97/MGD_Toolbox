%% Test Detection of Synergy and Redundancy using Conditioning
clear all; close all; clc;

% Fix Seed for Reproducibility
rng('default');

% Noise Vector
x=randn(1000,4);

% Equations X_3 and X_4
for t=2:1000
    x(t,4)=0.9*x(t-1,3)+0.1*x(t,4);
    x(t,3)=0.5*x(t-1,1)+0.5*x(t-1,2)+0.1*x(t,3);
end

% Parameters
% Driver and Target 
i_driver=1;
i_target=4;
model_order=2;
q=22; % Auto correlation truncation lag
numsurr=100;
alpha= 0.05;

%% Estimation of VAR Model
[Am,Su,Yp,Up,Z,Yb]=mgd_idMVAR(x',model_order,0);

%% Pairwise MIR
ret = mgd_MIR_lin(Am,Su,q,i_driver,i_target);

%% Conditional MIR
Ixy_z=mgd_cMIR_lin(Am,Su,q,i_driver,i_target,2);

%% Greedy Search
ret_syn_red=mgd_mir_syn_red_est(x',model_order,i_driver,i_target,q,numsurr,alpha);
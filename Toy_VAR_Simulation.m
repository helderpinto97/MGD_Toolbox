%% Test Detection of Synergy and Redundancy using Conditioning
clear all; close all; clc;

% Fix Seed for Reproducibility
rng('default');
addpath('functions');

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
numsurr=10000;
alpha= 0.05;

%% Estimation of VAR Model
[Am,Su,Yp,Up,Z,Yb]=mgd_idMVAR(x',model_order,0);

%% Pairwise MIR
ret = mgd_MIR_lin(Am,Su,q,i_driver,i_target);

%% Conditional MIR
Ixy_z=mgd_cMIR_lin(Am,Su,q,i_driver,i_target,2);

%% Greedy Search
ret_syn_red=mgd_mir_syn_red_est(x',model_order,i_driver,i_target,q,numsurr,alpha);
eMIRm=ret_syn_red.MIRm;
eMIRM=ret_syn_red.MIRM;
ekkm=ret_syn_red.kkm;
ekkM=ret_syn_red.kkM;

% Display Estimated Results 
disp('%%%%% Estimated Results %%%%%%%%%%%%%%');
disp(['Pairwise MIR=' num2str(ret_syn_red.MIRp)])
disp(['Fully-conditioned MIR=' num2str(ret_syn_red.MIRf)])
disp(['Decomposition MIRM=S+R+U:' num2str(ret_syn_red.S+ret_syn_red.U+ret_syn_red.R)]);
disp(['MIRM = ' num2str(ret_syn_red.MIRM(end))]);
disp(['S = ' num2str(ret_syn_red.S)]);
disp(['R = ' num2str(ret_syn_red.R)]);
disp(['U = ' num2str(ret_syn_red.U)]);

epsi=10^-5;
ymin=0.9*min([eMIRm eMIRM]);
ymax=1.1*max([eMIRm eMIRM]);
if abs(ymin-ymax)<epsi, ymin=0; ymax=0.1; end

figure('WindowState','maximized');
subplot(1,2,1);
plot(eMIRm,'ob-','LineWidth',1.5); hold on
ylabel(['eMIR_{' int2str(i_driver) ';' int2str(i_target) '|m}'])
xmin=0.5; xmax=length(ekkm)+0.5; xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('k_m')
label_m{1}='[ ]';
for i=2:length(ekkm)
    label_m{i}=['X_' int2str(ekkm(i))];
    plot(i*ones(numsurr,1),ret_syn_red.MIRm_surr(:,i-1),'x','MarkerSize',10,'Color','k');
end
xticks(1:length(eMIRm));
xticklabels(label_m);
title('MIR minimization');
ax=gca;
ax.FontSize=18;

subplot(1,2,2);
plot(eMIRM,'or-','LineWidth',1.5); hold on;
ylabel(['eMIR_{' int2str(i_driver) ';' int2str(i_target) '|M}'])
xmin=0.5; xmax=length(ekkM)+0.5; xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('k_M');
label_M{1}='[ ]';
for i=2:length(ekkM)
    label_M{i}=['X_' int2str(ekkM(i))];
    plot(i*ones(numsurr,1),ret_syn_red.MIRM_surr(:,i-1),'x','MarkerSize',10,'Color','k')
end
xticks(1:length(eMIRM));
xticklabels(label_M);
title('MIR maximization');
ax=gca;
ax.FontSize=18;
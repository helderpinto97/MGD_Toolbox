%% MIR decomposition with greedy approach - estimation from data

% Y, M*N matrix of time series (each time series is in a row)
% p, model order
% q: truncantion lag for correlations
% ii: index of first process (X)
% jj: index of second process (Y)
% numsur: number of surrogates
% alpha: level of significance
% epsi: accepted error

function ret = mir_syn_red_est_old(Y,p,ii,jj,q,numsurr,alpha,epsi)

if nargin<8, epsi=10^-5; end

[eAm,eSu]=lrp_idMVAR(Y,p);

M=size(eAm,1);
p=size(eAm,2)/M;
[eAmt,eSut,ki]=VARreorder(eAm,eSu,ii,jj);

% X and Y in the first two entries

%%% pairwise and fully conditioned MIR
eMIRp=MIR_lin(eAmt,eSut,q,1,2); % pairwise MIR

%%% Surrogates -- Check pairwise  is significant
% for s=1:100
%     Y_surr=Y;
%     Y_j_shift= circshift(Y(jj,:),50);
%     Y_surr(jj,:)=Y_j_shift;
%     [eAm_surr,eSu_surr]=lrp_idMVAR(Y_surr,p);
%     M=size(eAm_surr,1);
%     p=size(eAm_surr,2)/M;
%     [eAmt_surr,eSut_surr,~]=VARreorder(eAm_surr,eSu_surr,ii,jj);
%     eMIRp_surr=MIR_lin(eAmt_surr,eSut_surr,q,1,2);
%     tmp_surr_pair(s)=eMIRp_surr;
% end

eMIRf=cMIR_lin(eAmt,eSut,q,1,2,3:M); % fully conditioned MIR

%verification of reordering
% GCp_ver=lrp_cGC(Am,Su,ii,jj,[],q); % pairwise TE
% GCf_ver=lrp_cGC(Am,Su,ii,jj,kk,q); % fully conditioned TE
% [GCp GCp_ver GCf GCf_ver]

%% greedy search for minimum
exitcrit=0;
cm=[];
cset=3:M; %initial conditioning set
eMIRm=eMIRp; %initial min = pairwise GC
MIRm_surr=[];
while exitcrit==0
    MIRtmp=nan*ones(length(cset),1);
    for i=1:length(cset)
        MIRtmp(i)=cMIR_lin(eAmt,eSut,q,1,2,[cm cset(i)]);
    end
    [MIRmin,imin]=min(MIRtmp);
    imin=imin(1);
    
    if MIRmin>eMIRp
        break
    end
    
    MIRtmps=nan*ones(numsurr,1);
    for is=1:numsurr
        y=Y(ki(cset(imin)),:)';
        ys=surr_iaafft(y);
        Ys=Y; Ys(ki(cset(imin)),:)=ys';
        [seAm,seSu]=lrp_idMVAR(Ys,p); 
        [seAmt,seSut,ki]=VARreorder(seAm,seSu,ii,jj);
        MIRtmps(is)=cMIR_lin(seAmt,seSut,q,1,2,[cm cset(imin)]);
    end
    MIRtmps_th=prctile(MIRtmps,100*alpha);
    MIRm_surr=[MIRm_surr MIRtmps];
    
    eMIRm=[eMIRm MIRmin];
    cm=[cm cset(imin)];
    if MIRtmps_th>MIRmin % MIR has decreased non-negligibly (surrogates)
        cset(imin)=[];
        if isempty(cset), exitcrit=1; end
    else %no decrease
        cm(end)=[];
        eMIRm(end)=[];
        exitcrit=1;
    end   
end
kkm=0;
for i=1:length(cm)
    kkm=[kkm ki(cm(i))];
end

%% greedy search for maximum
exitcrit=0;
cM=[];
cset=3:M; %initial conditioning set
eMIRM=eMIRp; %initial min = pairwise MIR
MIRM_surr=[];
while exitcrit==0
    MIRtmp=nan*ones(length(cset),1);
    for i=1:length(cset)
        MIRtmp(i)=cMIR_lin(eAmt,eSut,q,1,2,[cM cset(i)]);
    end
    [MIRmax,imax]=max(MIRtmp);
    imax=imax(1);

    if MIRmax<eMIRp
        break
    end
    
    MIRtmps=nan*ones(numsurr,1);
    for is=1:numsurr
        y=Y(ki(cset(imax)),:)';
        ys=surr_iaafft(y);
        Ys=Y; Ys(ki(cset(imax)),:)=ys';
        [seAm,seSu]=lrp_idMVAR(Ys,p); 
        [seAmt,seSut,ki]=VARreorder(seAm,seSu,ii,jj);
        MIRtmps(is)=cMIR_lin(seAmt,seSut,q,1,2,[cM cset(imax)]);
    end
    MIRtmps_th=prctile(MIRtmps,100*(1-alpha));
    MIRM_surr=[MIRM_surr MIRtmps];   

    eMIRM=[eMIRM MIRmax];
    cM=[cM cset(imax)];
    if MIRmax>MIRtmps_th % % MIR has increased non-negligibly (surrogates)
        cset(imax)=[];
        if isempty(cset), exitcrit=1; end
    else %no increase
        cM(end)=[];
        eMIRM(end)=[];
        exitcrit=1;
    end   
end
kkM=0;
for i=1:length(cM)
    kkM=[kkM ki(cM(i))];
end

R=eMIRp-eMIRm(end);
S=eMIRM(end)-eMIRp;
U=eMIRm(end);
GM=eMIRM(end);


assert(abs(eMIRM(end)-(S+R+U))<epsi);

%% output
ret.MIRp=eMIRp;
ret.MIRf=eMIRf;
ret.GM=GM;
ret.R=R;
ret.S=S;
ret.U=U;
ret.MIRm=eMIRm;
ret.kkm=kkm;
ret.MIRM=eMIRM;
ret.kkM=kkM;
ret.MIRm_surr=MIRm_surr;
ret.MIRM_surr=MIRM_surr;

end
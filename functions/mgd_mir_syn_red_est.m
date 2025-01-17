%% MIR decomposition with greedy approach - estimation from data

%% V2: Perform surrogates to assess the statistical significance of the pairwise measure:
% - If not significant dont search for minimum

% Y, M*N matrix of time series (each time series is in a row)
% p, model order
% q: truncantion lag for correlations
% ii: index of first process (X)
% jj: index of second process (Y)
% numsur: number of surrogates
% alpha: level of significance
% epsi: accepted error

function ret = mgd_mir_syn_red_est(Y,p,ii,jj,q,numsurr,alpha,epsi)

if nargin<8, epsi=10^-5; end

[eAm,eSu]=mgd_idMVAR(Y,p);

M=size(eAm,1);
p=size(eAm,2)/M;
[eAmt,eSut,ki]=mgd_VARreorder(eAm,eSu,ii,jj);

% X and Y in the first two entries

%%% pairwise and fully conditioned MIR
eMIRp=mgd_MIR_lin(eAmt,eSut,q,1,2); % pairwise MIR
eMIRf=mgd_cMIR_lin(eAmt,eSut,q,1,2,3:M); % fully conditioned MIR

%%% Surrogates -- Check pairwise  is significant
for s=1:numsurr
    Y_surr=Y;
    minshift = 40;
    maxshift=length(Y_surr)-minshift;
    lagshift=fix(rand*(maxshift-minshift+1)+minshift);
    Y_j_shift= circshift(Y(jj,:),lagshift);
    Y_surr(jj,:)=Y_j_shift;
    [eAm_surr,eSu_surr]=mgd_idMVAR(Y_surr,p);
    M=size(eAm_surr,1);
    p=size(eAm_surr,2)/M;
    [eAmt_surr,eSut_surr,~]=mgd_VARreorder(eAm_surr,eSu_surr,ii,jj);
    eMIRp_surr=mgd_MIR_lin(eAmt_surr,eSut_surr,q,1,2);
    tmp_surr_pair(s)=real(eMIRp_surr);
end

%% Greedy search for minimum
if eMIRp<=prctile(tmp_surr_pair,(1-alpha)*100)
    % So Redudancy is not statistical significant
    eMIRm=eMIRp;
    eMIRp_sig=0;
    eMIRm_sig=0;
    R_sig=0;
    
    % No process added to the conditioning vector
    kkm=0;
    MIRm_surr=0;

else
    eMIRp_sig=1; % Pairwise is significant so search for minimum
    exitcrit=0;
    cm=[];
    cset=3:M; %initial conditioning set
    eMIRm=eMIRp; %initial min = pairwise GC
    MIRm_surr=[];
    while exitcrit==0
        MIRtmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            MIRtmp(i)=mgd_cMIR_lin(eAmt,eSut,q,1,2,[cm cset(i)]);
        end
        [MIRmin,imin]=min(MIRtmp);
        imin=imin(1);

        if MIRmin>eMIRp
            break
        end

        MIRtmps=nan*ones(numsurr,1);
        for is=1:numsurr
            y=Y(ki(cset(imin)),:)';
            ys=mgd_surr_iaafft(y);
            Ys=Y; Ys(ki(cset(imin)),:)=ys';
            [seAm,seSu]=mgd_idMVAR(Ys,p);
            [seAmt,seSut,ki]=mgd_VARreorder(seAm,seSu,ii,jj);
            MIRtmps(is)=real(mgd_cMIR_lin(seAmt,seSut,q,1,2,[cm cset(imin)]));
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

    if length(kkm)>1  % at leat one process that decreases significantly the MIR
        eMIRm_sig=1;
        R_sig=1;
    else
        eMIRm_sig=0;
        R_sig=0;
    end
end


%% Greedy search for maximum
exitcrit=0;
cM=[];
cset=3:M; %initial conditioning set
eMIRM=eMIRp; %initial min = pairwise MIR
MIRM_surr=[];
while exitcrit==0
    MIRtmp=nan*ones(length(cset),1);
    for i=1:length(cset)
        MIRtmp(i)=mgd_cMIR_lin(eAmt,eSut,q,1,2,[cM cset(i)]);
    end
    [MIRmax,imax]=max(MIRtmp);
    imax=imax(1);

    if MIRmax<eMIRp
        break
    end
    
    MIRtmps=nan*ones(numsurr,1);
    for is=1:numsurr
        y=Y(ki(cset(imax)),:)';
        ys=mgd_surr_iaafft(y);
        Ys=Y; Ys(ki(cset(imax)),:)=ys';
        [seAm,seSu]=mgd_idMVAR(Ys,p); 
        [seAmt,seSut,ki]=mgd_VARreorder(seAm,seSu,ii,jj);
        MIRtmps(is)=real(mgd_cMIR_lin(seAmt,seSut,q,1,2,[cM cset(imax)]));
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

% Verify significance of synergy -- at least one process increase
% significantly the MIR

if length(kkM)>1
    eMIRM_sig=1;
    S_sig=1;
else
    eMIRM_sig=0;
    S_sig=0;
end


%% Decomposition
R=eMIRp-eMIRm(end);
S=eMIRM(end)-eMIRp;
U=eMIRm(end);


%% Verify if decomposition holds
assert(abs(eMIRM(end)-(S+R+U))<epsi);

%% output

% Pairwise MIR
ret.MIRp=eMIRp;
ret.MIRp_sig=eMIRp_sig;

% Fully Conditioned MIR
ret.MIRf=eMIRf;

% Synergy
ret.R=R;
ret.R_sig=R_sig;

ret.S=S;
ret.S_sig=S_sig;

% Unique Information
ret.U=U;
ret.U_sig=eMIRm_sig;

% Minimum MIR
ret.MIRm=eMIRm;
ret.MIRm_sig=eMIRm_sig;
ret.kkm=kkm;

% Maximum MIR
ret.MIRM=eMIRM;
ret.MIRM_sig=eMIRM_sig;
ret.kkM=kkM;

% Percentiles of surrogates used on the test to find redundancy and synergy
ret.MIRm_surr=MIRm_surr;
ret.MIRM_surr=MIRM_surr;

end
%% MIR mM decomposition with greedy approach - theoretical values

%%% input: VAR parameters (Am, Su)
%%% input: index of driver (ii) and target (jj) processes
%%% input: lag over which correlations are computed via YW eq. (q)

function ret = mgd_mir_syn_red_th(Am,Su,ii,jj,q,epsi)

if nargin<6, epsi=10^-5; end

M=size(Am,1);
p=size(Am,2)/M;
[Amt,Sut,ki]=mgd_VARreorder(Am,Su,ii,jj);

%%% pairwise and fully conditioned MIR
MIRp=mgd_MIR_lin(Amt,Sut,q,1,2); % pairwise MIR
MIRf=mgd_cMIR_lin(Amt,Sut,q,1,2,3:M); % fully conditioned MIR

%verification of reordering
% GCp_ver=lrp_cGC(Am,Su,ii,jj,[],q); % pairwise TE
% GCf_ver=lrp_cGC(Am,Su,ii,jj,kk,q); % fully conditioned TE
% [GCp GCp_ver GCf GCf_ver]

%% greedy search for minimum
exitcrit=0;
cm=[];
cset=3:M; %initial conditioning set
MIRm=MIRp; %initial min = pairwise MIR
while exitcrit==0
    MIRtmp=nan*ones(length(cset),1);
    for i=1:length(cset)
        MIRtmp(i)=mgd_cMIR_lin(Amt,Sut,q,1,2,[cm cset(i)]);
    end
    [MIRmin,imin]=min(MIRtmp);
    imin=imin(1);
    if MIRm(end)-MIRmin>epsi % MIR has decreased non-negligibly (surrogates for non theor)
        MIRm=[MIRm MIRmin];
        cm=[cm cset(imin)];
        cset(imin)=[];
        if isempty(cset), exitcrit=1; end
    else %no decrease
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
MIRM=MIRp; %initial min = pairwise GC
while exitcrit==0
    MIRtmp=nan*ones(length(cset),1);
    for i=1:length(cset)
        MIRtmp(i)=mgd_cMIR_lin(Amt,Sut,q,1,2,[cM cset(i)]);
    end
    [MIRmax,imax]=max(MIRtmp);
    imax=imax(1);
    if MIRmax-MIRM(end)>epsi % GC has increased non-negligibly (surrogates for non theor)
        MIRM=[MIRM MIRmax];
        cM=[cM cset(imax)];
        cset(imax)=[];
        if isempty(cset), exitcrit=1; end
    else %no increase
        exitcrit=1;
    end   
end
kkM=0;
for i=1:length(cM)
    kkM=[kkM ki(cM(i))];
end

R=MIRp-MIRm(end);
S=MIRM(end)-MIRp;
U=MIRm(end);
GM=MIRM(end);

assert(abs(MIRM(end)-(S+R+U))<epsi);

%% output
ret.MIRp=MIRp;
ret.MIRf=MIRf;
ret.GM=GM;
ret.R=R;
ret.S=S;
ret.U=U;
ret.MIRm=MIRm;
ret.kkm=kkm;
ret.MIRM=MIRM;
ret.kkM=kkM;

end

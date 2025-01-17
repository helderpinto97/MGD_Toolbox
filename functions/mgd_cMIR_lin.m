function Ixy_z=mgd_cMIR_lin(Am,Su,q,ix,iy,iz)
retx_z=mgd_MIR_lin(Am,Su,q,ix,iz); %retx_z.Ixy; %MIR btw processes ix and iz
retx_yz=mgd_MIR_lin(Am,Su,q,ix,[iy iz]); %retx_yz.Ixy; %MIR btw processes ix and [iy,iz]
Ixy_z=retx_yz-retx_z; %cMIR btw processes ix and iy given processes iz
end
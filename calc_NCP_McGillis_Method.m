function [BMS] = calc_NCP_McGillis_Method(BMS, ztop, zbtm, Sal_est)

% NCP Calculations using McGillis et al. 2011 method 

% created by Michelle Platz 
% University of South Floruda 
% 12/18/2019

% Variables:
%     Utop = Horizontal velocity at the top pump height(m/s)
%     Ubtm = Horizontal velocity at the bottom pump height (m/s)
%     ztop = height of the top pump (m)
%     zbtm = height of the bottom pump (m)
%     rho = density of water (kg/m3)
%     ustar = friction velocity (m/s) (calculated from ustar_McGillis_Method)
%     DOtop = Concentration at the top pump (umol/kg)
%     DObtm = Concentration at the bottom pump (umol/kg)
%     Sal_est = average salinity from SP ((g/kg), or ppt)

% Output = J_NCP (mmol/m2/hr)

%NCP_McGillis = (rho*ustar*0.41)*((DOtop-DObtm)/(ln(ztop/zbtm)))


%calculate density of seawater at each timestep 
BMS.DENS = gsw_rho(ones(size(BMS.TC)).*Sal_est, BMS.TC, BMS.Pres);

% presize matrix
BMS.JDO_WM = NaN(size(BMS.SDN));
BMS.NEP_WM = NaN(size(BMS.SDN));

% DO loop: loop through each gradient 
for i = 1:length(BMS.SDN)
    
    rho = BMS.DENS(i);
    ustar = BMS.ustar_WM(i);
    DOtop = BMS.DOXY(1,i);
    DObtm = BMS.DOXY(2,i);
    
    BMS.JDO_WM(i) = (rho*ustar*0.41)*((DOtop-DObtm)/(log(ztop/zbtm)));
    
    BMS.NEP_WM(i) = -BMS.JDO_WM(i) * 3600/1000; % Convert flux [umol/m2/sec] to NEP [mmol/m2/hr]

end 
return 

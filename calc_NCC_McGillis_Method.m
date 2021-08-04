function [BMS] = calc_NCC_McGillis_Method(BMS, ztop, zbtm, Sal_est)

% NCC Calculations using McGillis et al. 2011 method 

% created by Michelle Platz 
% University of South Floruda 
% 12/18/2019

% Variables:
%     Utop = Horizontal velocity at the top pump height(m/s)
%     Ubtm = Horizontal velocity at the bottom pump height (m/s)
%     ztop = height of the top pump (m)
%     zbtm = height of the bottom pump (m)
%     rho = density of water (kg/m3)
%     ustar = friction velocity (m/s)(calculated from ustar_McGillis_Method)
%     TAtop = Concentration at the top pump (umol/kg)
%     TAbtm = Concentration at the bottom pump (mmol/kg)
%     Sal_est = average salinity from SP ((g/kg), or ppt)

% Output = J_NCC (umol/m2/hr)

%NCC_McGillis = (rho*ustar*0.41)*((TAtop-TAbtm)/(ln(ztop/zbtm)))


%calculate density of seawater at each timestep 
BMS.DENS = gsw_rho(ones(size(BMS.TC)).*Sal_est, BMS.TC, BMS.Pres);

% presize matrix
BMS.JTA_WM = NaN(size(BMS.SDN));
BMS.NEC_WM = NaN(size(BMS.SDN));

% DO loop: loop through each gradient 
for i = 1:length(BMS.SDN)
    
    rho = BMS.DENS(i);
    ustar = BMS.ustar_WM(i);
    TAtop = BMS.TAtop(i);
    TAbtm = BMS.TAbtm(3,i); %index for Q = 1 
    
    BMS.JTA_WM(i) = (rho*ustar*0.41)*((TAtop-TAbtm)/(log(ztop/zbtm)));
    
    BMS.NEC_WM(i) = BMS.JTA_WM(i) * 3600/2/1000; % Convert flux [umol/m2/sec] to NEC [mmol/m2/hr] and account for 2meq

end 
return 

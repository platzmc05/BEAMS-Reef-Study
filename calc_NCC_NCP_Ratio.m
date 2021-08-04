function NCC_NCP_Ratio = calc_NCC_NCP_Ratio (ADavg)

% Calculates hourly NCC/NCP ratios  

% ADavg is the structure containing Aquadopp data, averaged over 60 minutes (created by average_aquadopp.m)

% Calculated values are saved into the ADavg structure.

% Created by: Michelle Platz 
% USF
% Version updated 8/27/2019


NCC_NCP_Ratio = NaN(size(ADavg.SDN));


for i = 1:length(ADavg.SDN)
    
    NCC_NCP_Ratio(i) = ADavg.NEC(i)./ADavg.NEP(i);
    
    
    
end 
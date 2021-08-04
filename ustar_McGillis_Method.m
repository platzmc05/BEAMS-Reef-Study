function [ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)

%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)

% Calculates friction velocity from current profile measured by HR Aquadopp. 

% ADavg is the structure containing Aquadopp data,
% ztop = height of top pump (m)
% zbtm = height of bottom pump (m)
% bintop = row index for the data corresponding to the top pump height - adjusted .25m for cinder block under ADCP
% binbtm = row index for the data corresponding to the bottom pump height - adjusted .25m for cinder block under ADCP

% Created by Michelle Platz 
% USF
% Version 12/12/2019


% presize matrix
ADavg.ustar_WM = NaN(size(ADavg.SDN));
ADavg.Utop = ADavg.uv(bintop,:); % 
ADavg.Ubtm = ADavg.uv(binbtm,:);
i1m = find(ADavg.bin_depth <= 1, 1, 'last'); % index for 1 m above benthos. This will be used for U0. 

% loop through each profile
for i = 1:length(ADavg.SDN)
      
    
    ADavg.ustar_WM(i) = 0.41*((ADavg.Utop(i)-ADavg.Ubtm(i))/(log(ztop/zbtm)));
    % calculate ustar from velocity and pump intake heights  
    
         
% running average for ustar; window is 2m+1 (m = 2nd input), so ~1 hour.
% This can be adjusted based on your sampling interval.
ADavg.ustar_WM_runmean = runmean(ADavg.ustar_WM,3,1,'edge');


end

ADavg.U0 = ADavg.uv(i1m,:); %i1m = index for 1 m above benthos
ADavg.Cd_WM = NaN(size(ADavg.ustar_WM)); %Presize Matrix 
ADavg.Cd_WM = (ADavg.ustar_WM.^2)./ADavg.uv(i1m,:).^2'; %use U0 from 1 m

return



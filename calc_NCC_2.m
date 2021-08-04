function BMS = calc_NCC_2(BMS,z)

% Calculates hourly flux by fitting flux and C1 to concentration gradient
% measured by SP using ititial guesses. 

% z = heights of pumps = [pump1(top) pump2(bottom)]

% BMS is the structure containing SeaphOx data and Aquadopp data, averaged over 15 minutes (created by average_aquadopp.m)

% Calculated values are saved into the BMS structure.

% Created by: Michelle Platz  
% USF
% Version updated 3/4/2021


JTAguess = 1;
C1guess = BMS.TAtop(40);

% presize matrix
BMS.JTA = NaN(size(BMS.SDN));
BMS.C1 = NaN(size(BMS.SDN));

% TA loop: loop through each gradient 
for i = 1:length(BMS.SDN)
    
        ustar = BMS.ustar(i); % single value per timestep 
        d = 0.1; %displacement height
        a = [JTAguess, C1guess];
        
        TAz = [BMS.TAtop(i), BMS.TAbtm(3,i)]; % 2 values per timestep - ydata in lsqcurvefit
        %TAbtm row 3 is for Q=1
       
        if(isempty(find(isnan(TAz),1)))
            
            options = optimset('Display','none', 'diagnostics', 'off');
           
            % col 1: lb/ub of J, col2: lb/ub of C1
            lb = [-2, C1guess-10]; % lower bound for J and C1. 
            ub = [2, C1guess+10];    % upper bound for J and C1.
            
            % use lsqcurvefit to fit flux to concentration gradient data (TAz).
            [c, resnorm, resid, ~]  = lsqcurvefit(@flux_fit, a, z, TAz, lb, ub, options);
 
            % save values % total alk flux
            BMS.JTA(i) = c(1);
            BMS.NEC(i) = BMS.JTA(i) * 3600 / 2; % Convert alk flux [mmol/m2/sec] to NEC [mmol/m2/hr]
            BMS.C1(i) = c(2);
        else
            BMS.JTA(i) = NaN;
            BMS.C1(i) = NaN;
        end      
end


function F = flux_fit(a,z)


% 
% Nested function for fitting flux to concentration gradient data.
% a(1) = J (flux) 
% a(2) = C1 (integration constant)

% z is the two heights of the pumps.
% d is the displacement height of the coral (25 or 30 CM) 
% ustar is friction velocity. 


% Created by: Michelle Platz 
% USF
% Version 8-7-2019

F = (a(1)/(ustar*0.41))*log(z/d)+a(2);
end


end  
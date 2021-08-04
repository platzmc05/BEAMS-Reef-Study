function BMS = calc_NCP_2(BMS,z, C1guess)

% Calculates hourly flux by fitting flux and C1 to concentration gradient
% measured by SP using ititial guesses. 

% z = heights of pumps = [pump1(top) pump2(bottom)]

% BMS is the structure containing SeaphOx data and Aquadopp data, averaged over 15 minutes (created by average_aquadopp.m)

% Calculated values are saved into the BMS structure.

% Created by: Michelle Platz 
% USF
% Version updated 3/4/2021

JDOguess = 0;
C1guess = BMS.DOXY(1,C1guess); %initial guess DOXY value



% presize matrix
BMS.JDO = NaN(size(BMS.SDN));
BMS.C1 = NaN(size(BMS.SDN));

% DO loop: loop through each gradient 
for i = 1:length(BMS.SDN)
    
        ustar = BMS.ustar(i); % single value per timestep 
        d = 0.1; %height of coral placement - measuring relative to this height 
        a = [JDOguess, C1guess];
        
        DOz = [BMS.DOXY(1,i), BMS.DOXY(2,i)]; % 2 values per timestep, top to bottom - ydata in lsqcurvefit 
       
        if(isempty(find(isnan(DOz),1)))
            
            options = optimset('Display','none', 'diagnostics', 'off');
           
            % col 1: lb/ub of J, col2: lb/ub of C1
            lb = [-2, C1guess-10]; % lower bound for J and C1. 
            ub = [2, C1guess+10];    % upper bound for J and C1.
            
            % use lsqcurvefit to fit flux to concentration gradient data (DOz).
            [c, resnorm, resid, ~]  = lsqcurvefit(@flux_fit, a, z, DOz, lb, ub, options);
 
            % save values % total alk flux
            BMS.JDO(i) = c(1);
            BMS.NEP(i) = BMS.JDO(i) * 3600; % Convert alk flux [mmol/m2/sec] to NEP [mmol/m2/hr]
            BMS.C1(i) = c(2);
        else
            BMS.JDO(i) = NaN;
            BMS.C1(i) = NaN;
        end      
end


function F = flux_fit(a,z)


% 
% Nested function for fitting flux to concentration gradient data.
% a(1) = J (flux) 
% a(2) = C1 (integration constant)

% z is the two heights of the pumps. top to bottom
% d is the displacement height of the coral (25 or 30 CM) 
% ustar is friction velocity. 


% Created by: Michelle Platz 
% USF
% Version 8-7-2019

F =-(a(1)/(ustar*0.41))*log(z/d)+a(2);
end


end  
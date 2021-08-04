function [ADavg] = ustar_from_aquadopp2(ADavg, zrange, Sal_est)

% [ADavg] = ustar_from_aquadopp(ADavg, zrange, Sal_est);
%
% Calculates friction velocity from current profile measured by HR
% Aquadopp. ADavg is the structure containing Aquadopp data, averaged over
% X minutes (created by average_aquadopp.m), zrange is the depth range to
% use for the ustar calculation [zmin zmax], and Sal_est is the estimate
% for salinity (used for viscosity calculation). 
%
% Calculates drag coefficient based on U0 at 1 m, running average of ustar,
% viscosity, and Reynolds number. Calculated values are saved
% into the ADavg structure.
%
% Created by: Yui Takeshita
% MBARI
% Original Dec-6-2018
% Edited by Michelle Platz 
% USF
% Version 7/29/2019

% calculate ustar by fitting ustar, d, and z0 for 

% inital guesses for fit. ustar is estimated from an drag coefficient from
% Takeshita et al. 2016. It may need to be adjusted for less turbulent
% systems. 
dguess = 0.1;
z0guess = 0.001;
% a = [ustarguess, dguess, z0guess];
lb = [0, 0, 0.0000001]; % lower bound for ustar, d, and z0. 
ub = [1, 0.8, 0.03];    % upper bound for ustar, d, and z0.
iz = inrange(ADavg.bin_depth, zrange); %logical output: 1 if in range, 0 if not 
% index for depth range to use for fit
i1m = find(ADavg.bin_depth <= 1, 1, 'last'); % index for 1 m above benthos. This will be used for U0. 

% presize matrix
ADavg.ustar = NaN(size(ADavg.SDN));
ADavg.d = ADavg.ustar;
ADavg.z0 = ADavg.ustar;
ADavg.vel_fit = NaN(size(ADavg.uv)); % fitted velocity profile 

% loop through each profile
for i = 1:length(ADavg.SDN)
    
     % use data points that are not NaNs.
    inonan = ~isnan(ADavg.uv(iz,i)); % returns 42x1 logical of vel values in range
        %1 if not NaN
        %0 if NaN 


    % if current profile has more than 3 points in the desired range, do fit:
    if(length(find(inonan))>3)
        % ommit lowest point; large error
        inonan(1) = false; %set first point to zero 
                %find index when inonan is NaN - 42 points long
%         iz2 = find(inonan==false);
%         % indicate this value is NaN in iz - 108 points long
%         iz(iz2)=false;% iz now indicates only points that are both nonNaN and in range 
         
        vel = ADavg.uv(inonan,i); %ydata in lsqcurvefit - velocity 
        z = ADavg.bin_depth(inonan); %variable xdata - depths of bins in vel profile

        %ustarguess = sqrt(0.019);%ADavg.uv(i1m,i)*0.13; % , which is drag coefficient from Takeshita et al. 2016
        if (isnan(ADavg.uv(i1m,i))==1)% if obs at 1m is NaN
            ustarguess = sqrt(0.019);% drag coefficient from Takeshita et al. 2016
        else 
            ustarguess = ADavg.uv(i1m,i)*0.13; 
        end        
        
        a = [ustarguess, dguess, z0guess];
        % use lsqcurvefit to fit current profile data to law of the wall.
        options = optimset('Display','none', 'diagnostics', 'off');
                                % lsqcurvefit(@func, in.vals, xdata, ydata, lb, ub, options);
        [c, resnorm, resid, ~]  = lsqcurvefit(@lawofwall, a, z, vel, lb, ub, options);
             
%       % save values
        ADavg.ustar(i) = c(1);
        ADavg.d(i) = c(2);
        ADavg.z0(i) = c(3);
        ADavg.vel_fit(inonan,i) = lawofwall(c,z);       
   else
        % if there wasn't enough data in the profile, make it all NaN.
        ADavg.ustar(i) = NaN;
        ADavg.d(i) = NaN;
        ADavg.z0(i) = NaN;
        ADavg.vel_fit(:,i) = NaN(size(ADavg.vel_fit(:,i)));
    end
end

% running average for ustar; window is 2m+1 (m = 2nd input), so ~1 hour.
% This can be adjusted based on your sampling interval.
ADavg.ustar_runmean = runmean(ADavg.ustar,3,1,'edge');

% Some bonus calculations. Not sure how accurate or useful nu and Re are.
% Calcualte drag coefficient (at 1m), viscocity, and reynolds number

ADavg.U0 = ADavg.uv(i1m,:); %i1m = index for 1 m above benthos
ADavg.direction1m = ADavg.direction(i1m,:); %i1m = index for 1 m above benthos

ADavg.Cd = NaN(size(ADavg.ustar)); %Presize Matrix 
ADavg.nu = ADavg.Cd;%Presize Matrix
ADavg.Re = NaN(length(ADavg.bin_depth), length(ADavg.Cd)); %Presize Matrix

ADavg.Cd = (ADavg.ustar.^2)./ADavg.uv(i1m,:).^2'; %use U0 from 1 m
ADavg.nu = SW_Viscosity(ADavg.TC, 'C', ones(size(ADavg.TC)).*Sal_est, 'ppt')./(gsw_rho(ones(size(ADavg.TC)).*Sal_est, ADavg.TC, ADavg.Pres));
ADavg.Re = (ADavg.bin_depth*ADavg.ustar(:)')./(ones(length(ADavg.bin_depth),1)*ADavg.nu');


return



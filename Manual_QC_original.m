


% manual_QC_mote


clear
close all


% load('Mote_Aug2018_BEAMS_adp_01m_05m_pump_45_20.mat', 'BMS');
load('Combined_BMS_Aug2018_mote.mat', 'BMS');

% get DOXY std
idoxystd = find(BMS.DOXYstd > 0.8);
ihighdoxystd = [];
for i = 1:length(idoxystd)
    
    ihighdoxystd = vertcat(ihighdoxystd,[idoxystd(i)-1:1:idoxystd(i)+1]');
end
% get unique IDs
ihighdoxystd = unique(ihighdoxystd);
% remove 0's and out of index values
ihighdoxystd(ihighdoxystd==0) = [];
ihighdoxystd(ihighdoxystd> length(BMS.SDN)) = [];

% make it into index
trex = false(size(BMS.SDN));
trex(ihighdoxystd) = true;
ihighdoxystd = trex;
clear trex;

% plot to see what got removed
figure
hold on; box on;
plot(BMS.SDN, BMS.NEP, 'k');
plot(BMS.SDN(ihighdoxystd), BMS.NEP(ihighdoxystd), 'ro');
datetickzoom;

BMS.NEP(ihighdoxystd) = NaN;
BMS.NEC(:,ihighdoxystd) = NaN;
% Bin data

X = floor(nanmin(BMS.SDN)):1/24:ceil(nanmax(BMS.SDN));

BMSbin.SDN = X;
% bin data hourly. Vector variables first
BMSbin.PRES_AD  = bin_data_to_X_GF(BMS.SDN, BMS.PRES_AD, X);
BMSbin.U0       = bin_data_to_X_GF(BMS.SDN, BMS.U0, X);
BMSbin.DIR      = bin_data_to_X_GF(BMS.SDN, BMS.DIR, X);
BMSbin.ustar    = bin_data_to_X_GF(BMS.SDN, BMS.ustar, X);
BMSbin.ustar_rm = bin_data_to_X_GF(BMS.SDN, BMS.ustar_rm, X);
BMSbin.NEP      = bin_data_to_X_GF(BMS.SDN, BMS.NEP, X);
BMSbin.NEP_WM   = bin_data_to_X_GF(BMS.SDN, BMS.NEP_WM, X);

% variables from 3 differnet heights
BMSbin.TC(1,:)       = bin_data_to_X_GF(BMS.SDN, BMS.TC(1,:), X);
BMSbin.TC(2,:)       = bin_data_to_X_GF(BMS.SDN, BMS.TC(2,:), X);
% BMSbin.TC(3,:)       = bin_data_to_X_GF(BMS.SDN, BMS.TC(3,:), X);
BMSbin.DOXY(1,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DOXY(1,:), X);
BMSbin.DOXY(2,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DOXY(2,:), X);
% BMSbin.DOXY(3,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DOXY(3,:), X);
BMSbin.pH(1,:)       = bin_data_to_X_GF(BMS.SDN, BMS.pH(1,:), X);
BMSbin.pH(2,:)       = bin_data_to_X_GF(BMS.SDN, BMS.pH(2,:), X);
% BMSbin.pH(3,:)       = bin_data_to_X_GF(BMS.SDN, BMS.pH(3,:), X);
BMSbin.PSAL(1,:)     = bin_data_to_X_GF(BMS.SDN, BMS.PSAL(1,:), X);
BMSbin.PSAL(2,:)     = bin_data_to_X_GF(BMS.SDN, BMS.PSAL(2,:), X);
% BMSbin.PSAL(3,:)     = bin_data_to_X_GF(BMS.SDN, BMS.PSAL(3,:), X);
BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(BMS.SDN, BMS.O2SATPER(1,:), X);
BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(BMS.SDN, BMS.O2SATPER(2,:), X);
% BMSbin.O2SATPER(3,:) = bin_data_to_X_GF(BMS.SDN, BMS.O2SATPER(3,:), X);
BMSbin.Pres(1,:)     = bin_data_to_X_GF(BMS.SDN, BMS.Pres(1,:), X);
BMSbin.Pres(2,:)     = bin_data_to_X_GF(BMS.SDN, BMS.Pres(2,:), X);
% BMSbin.Pres(3,:)     = bin_data_to_X_GF(BMS.SDN, BMS.Pres(3,:), X);  
BMSbin.DENS(1,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DENS(1,:), X);
BMSbin.DENS(2,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DENS(2,:), X);
% BMSbin.DENS(3,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DENS(3,:), X);  
BMSbin.PAR(1,:)      = bin_data_to_X_GF(BMS.SDN, BMS.PAR(1,:), X);
BMSbin.PAR(2,:)      = bin_data_to_X_GF(BMS.SDN, BMS.PAR(2,:), X);
% BMSbin.PAR(3,:)      = bin_data_to_X_GF(BMS.SDN, BMS.PAR(3,:), X);  

% Now TA and NEC
for ii = 1:size(BMS.TA,1)
   BMSbin.TA(ii,1,:)    = bin_data_to_X_GF(BMS.SDN, squeeze(BMS.TA(ii,1,:)), X); 
   BMSbin.TA(ii,2,:)    = bin_data_to_X_GF(BMS.SDN, squeeze(BMS.TA(ii,2,:)), X); 
%    BMSbin.TA(ii,3,:)    = bin_data_to_X_GF(BMS.SDN, squeeze(BMS.TA(ii,3,:)), X);      
   BMSbin.NEC(ii,:)     = bin_data_to_X_GF(BMS.SDN, BMS.NEC(ii,:), X);
   BMSbin.NEC_WM(ii,:)  = bin_data_to_X_GF(BMS.SDN, BMS.NEC_WM(ii,:), X);
end
disp('Finished Metabolism Calculations');

save('Aug2018_mote_qcd.mat', 'BMS', 'BMSbin');



ibad = (BMSbin.DIR > 90 & BMSbin.DIR < 270) | BMSbin.U0 < 0.03;
nepuse = BMSbin.NEP;
nepuse(ibad) = NaN;
nepuse(nepuse < -7 | nepuse > 15) = NaN; % manual QC

necuse = BMSbin.NEC(3,:);
necuse(ibad) = NaN;

nepdbin = parse_to_diel(BMSbin.SDN, nepuse, 24);
necdbin = parse_to_diel(BMSbin.SDN, necuse, 24);


figure
hold on; box on;
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, nepdbin, 'ko', 'markersize', 3);
plot(1:24, necdbin, 'ro', 'markersize', 3);
plot(1:24, nanmedian(nepdbin,1), 'ko-');
plot(1:24, nanmedian(necdbin,1), 'ro-')
ylabel(['NEP or \color{red}NEC']);
title('August 2018');
xlabel('hour of day');





% 
% % additional QC
% ibad = (BMS.DIR > 90 & BMS.DIR < 270) | BMS.U0 < 0.03;
% 
% % bin to daily data
% irealbad = ibad(:) | ihighdoxystd(:);
% 
nepdiel = parse_to_diel(BMS.SDN(~irealbad), BMS.NEP(~irealbad), 24);
% 
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, -nepdiel, 'ro', 'markersize', 3);
% plot(1:24, -nanmean(nepdiel,1), 'ko-');
% ylabel('NEP');
% title('August 2018');
% xlabel('hour of day');
% 
% 
% % 
% % figure
% hold on; box on;
% plot(BMS.SDN, BMS.DOXY(1,:) - BMS.DOXY(2,:), 'k');
% plot(BMS.SDN(ihighdoxystd), BMS.DOXY(1,ihighdoxystd) - BMS.DOXY(2,ihighdoxystd), 'ro');
% 
% 
% 
% 
% 
% 
% ibad = (BMSbin.DIR > 90 & BMSbin.DIR < 270) | BMSbin.U0 < 0.03;
% 









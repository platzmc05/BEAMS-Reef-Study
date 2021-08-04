% August manual_QC_mote

close all 
clc

load('Mar2019_mote_qcd.mat');
Mar_BMS = BMS;
Mar_BMSbin = BMSbin;

ibad =  Mar_BMSbin.U0 < 0.03;
Mar_nepuse = Mar_BMSbin.NEP;
Mar_nepuse(ibad) = NaN;
Mar_nepuse(Mar_nepuse < -7 | Mar_nepuse > 15) = NaN; % manual QC

Mar_necuse = Mar_BMSbin.NEC(3,:);
Mar_necuse(ibad) = NaN;

Mar_nepdbin = parse_to_diel(Mar_BMSbin.SDN, Mar_nepuse, 24);
Mar_necdbin = parse_to_diel(Mar_BMSbin.SDN, Mar_necuse, 24);


figure
hold on; box on;
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, Mar_nepdbin, 'ko', 'markersize', 3);
plot(1:24, Mar_necdbin, 'ro', 'markersize', 3);
plot(1:24, nanmedian(Mar_nepdbin,1), 'ko-');
plot(1:24, nanmedian(Mar_necdbin,1), 'ro-')
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
% nepdiel = parse_to_diel(BMS.SDN(~irealbad), BMS.NEP(~irealbad), 24);
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









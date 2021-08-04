% March manual_QC_mote

close all 
clc

load('Mar2019_mote_qcd.mat')
Mar_BMS = BMS;
Mar_BMSbin = BMSbin;

% indexing data that was too slow (not turbulent boundary layer) 
Mar_ibad = Mar_BMSbin.NEP < -20 | Mar_BMSbin.NEC(3,:) < -20 | Mar_BMSbin.U0 < 0.03 ;% 2 CLEAR OUTLIERS REMOVED FROM THE DATASET 
Mar_igood = ~Mar_ibad;

Mar_nepuse = Mar_BMSbin.NEP;% good NEP data
Mar_nepuse(Mar_ibad) = NaN;

Mar_necuse = Mar_BMSbin.NEC(3,:);% good NEC data
Mar_necuse(Mar_ibad) = NaN;

Mar_nepuse_nonan = rmmissing(Mar_nepuse);
Mar_necuse_nonan = rmmissing(Mar_necuse);
Number_Mar_nepuse_nonan = length(Mar_nepuse_nonan)
Number_Mar_necuse_nonan = length(Mar_necuse_nonan)

%% COMPOSITE DIEL METABOLISM PLOTS
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
title('March 2019');
xlabel('hour of day');





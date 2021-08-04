% Reef Study Analysis 
% Michelle Platz - USF 
% 3/10/2021
clear all
close all
clc
%% Create Master Datasets - load post QC datasets

%load U datasets 
load('UJuly20_2.mat', 'UJuly_BMS', 'UJuly_BMSbin');
load('UAug120_2.mat', 'UAug1_BMS', 'UAug1_BMSbin');
load('UAug220_2.mat', 'UAug2_BMS', 'UAug2_BMSbin');
load('USept20_2.mat', 'USept_BMS', 'USept_BMSbin');
load('UOct20_2.mat', 'UOct_BMS', 'UOct_BMSbin');

%load M datasets 
load('MJuly20_2.mat', 'MJuly_BMS', 'MJuly_BMSbin');
load('MAug120_2.mat', 'MAug1_BMS', 'MAug1_BMSbin');
load('MAug220_2.mat', 'MAug2_BMS', 'MAug2_BMSbin');
load('MSept20_2.mat', 'MSept_BMS', 'MSept_BMSbin');
load('MOct20_2.mat', 'MOct_BMS', 'MOct_BMSbin');

%% Make Master Data Structures

% create a datestring by week for the entire length of the monitoring
DateString = {'07/01/2020 12:00:00';'07/08/2020 12:00:00';'07/15/2020 12:00:00';'07/22/2020 12:00:00';'07/29/2020 12:00:00';...
              '08/05/2020 12:00:00';'08/12/2020 12:00:00';'08/19/2020 12:00:00';'08/26/2020 12:00:00';...
              '09/02/2020 12:00:00';'09/09/2020 12:00:00';'09/16/2020 12:00:00';'09/23/2020 12:00:00';'09/30/2020 12:00:00';...
              '10/07/2020 12:00:00';'10/14/2020 12:00:00';'10/21/2020 12:00:00';'10/28/2020 12:00:00';...
              '11/04/2020 12:00:00';'11/11/2020 12:00:00'};
    
% formatIn = 'mm/dd/yyyy HH:MM:SS';
tick = datenum(DateString,'mm/dd/yyyy HH:MM:SS');
Xrange = [datenum('06-29-2020 16:15:00'), datenum('14-Nov-2020 10:00:00')];
Xrange2 = [datenum('06-29-2020 16:15:00'), datenum('14-Oct-2020 10:00:00')];

% Cudjoe
UMaster_BMSbin.SDN = [UJuly_BMSbin.SDN, UAug1_BMSbin.SDN, UAug2_BMSbin.SDN, USept_BMSbin.SDN, UOct_BMSbin.SDN];
UMaster_BMSbin.NEP_QC = [UJuly_BMSbin.NEP_QC, UAug1_BMSbin.NEP_QC, UAug2_BMSbin.NEP_QC, USept_BMSbin.NEP_QC, UOct_BMSbin.NEP_QC];
UMaster_BMSbin.NEC_QC = [UJuly_BMSbin.NEC_QC, UAug1_BMSbin.NEC_QC, UAug2_BMSbin.NEC_QC, USept_BMSbin.NEC_QC, UOct_BMSbin.NEC_QC];
UMaster_BMSbin.ustar = [UJuly_BMSbin.ustar, UAug1_BMSbin.ustar, UAug2_BMSbin.ustar, USept_BMSbin.ustar, UOct_BMSbin.ustar];
UMaster_BMSbin.U0 = [UJuly_BMSbin.U0, UAug1_BMSbin.U0, UAug2_BMSbin.U0, USept_BMSbin.U0, UOct_BMSbin.U0];

close all
figure
hold on; box on;
plot(UMaster_BMSbin.SDN, UMaster_BMSbin.NEP_QC, 'b-'); 
plot(UMaster_BMSbin.SDN, UMaster_BMSbin.NEC_QC, 'r-'); 
set(gca, 'xlim', Xrange2, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
xlabel('Date');
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
title('Cudjoe Hourly Binned Fluxes');


% M32 
MMaster_BMSbin.SDN = [MJuly_BMSbin.SDN, MAug1_BMSbin.SDN, MAug2_BMSbin.SDN, MSept_BMSbin.SDN, MOct_BMSbin.SDN];
MMaster_BMSbin.NEP_QC = [MJuly_BMSbin.NEP_QC, MAug1_BMSbin.NEP_QC, MAug2_BMSbin.NEP_QC, MSept_BMSbin.NEP_QC, MOct_BMSbin.NEP_QC];
MMaster_BMSbin.NEC_QC = [MJuly_BMSbin.NEC_QC, MAug1_BMSbin.NEC_QC, MAug2_BMSbin.NEC_QC, MSept_BMSbin.NEC_QC, MOct_BMSbin.NEC_QC];
MMaster_BMSbin.ustar = [MJuly_BMSbin.ustar, MAug1_BMSbin.ustar, MAug2_BMSbin.ustar, MSept_BMSbin.ustar, MOct_BMSbin.ustar];
MMaster_BMSbin.U0 = [MJuly_BMSbin.U0, MAug1_BMSbin.U0, MAug2_BMSbin.U0, MSept_BMSbin.U0, MOct_BMSbin.U0];

figure
hold on; box on;
M_NEPplot = plot(MMaster_BMSbin.SDN, MMaster_BMSbin.NEP_QC, 'b-'); 
M_NECplot = plot(MMaster_BMSbin.SDN, MMaster_BMSbin.NEC_QC, 'r-');
set(gca, 'xlim', Xrange2, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
xlabel('Date');
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
title('Marker 32 Hourly Binned Fluxes');

%% Cudjoe Site Characterizations Subplot

% Load full datasets
load('UJuly20_SiteChar.mat', 'UJuly_SiteChar')
load('UAug120_SiteChar.mat', 'UAug1_SiteChar' )
load('UAug220_SiteChar.mat', 'UAug2_SiteChar' )
load('USept20_SiteChar.mat', 'USept_SiteChar' )
load('UOct20_SiteChar.mat', 'UOct_SiteChar' )

% Cudjoe Master Site Char Datastructure
UMaster_SiteChar.SDN = [UJuly_SiteChar.SDN; UAug1_SiteChar.SDN; UAug2_SiteChar.SDN; USept_SiteChar.SDN; UOct_SiteChar.SDN];
UMaster_SiteChar.AD_SDN = [UJuly_SiteChar.AD_SDN; UAug1_SiteChar.AD_SDN; USept_SiteChar.AD_SDN];
UMaster_SiteChar.AD_TC = [UJuly_SiteChar.AD_TC; UAug1_SiteChar.AD_TC; USept_SiteChar.AD_TC];
UMaster_SiteChar.AD_Pres = [UJuly_SiteChar.AD_Pres; UAug1_SiteChar.AD_Pres; USept_SiteChar.AD_Pres];

% Cudjoe Subplot 
close all
figure
sgtitle('Cudjoe Ledge Reef Characterization Parameters')
subplot(5,1,1);%Temp
hold on; box on;
plot(UJuly_SiteChar.SDN, UJuly_SiteChar.TC, 'r')
plot(UAug1_SiteChar.AD_SDN, UAug1_SiteChar.AD_TC, 'r')
plot(USept_SiteChar.AD_SDN, USept_SiteChar.AD_TC, 'r')
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
% xlabel('Date');
ylabel('Temperature (C)');
title('Temperature');

UMaster_SiteChar.TC = [UJuly_SiteChar.TC; UAug1_SiteChar.AD_TC; USept_SiteChar.AD_TC];
UMaster_SiteChar.TC_median =  median(UMaster_SiteChar.TC)
    UMaster_SiteChar.TC_min = min(UMaster_SiteChar.TC)
    UMaster_SiteChar.TC_max = max(UMaster_SiteChar.TC)
UMaster_SiteChar.TC_mean = mean(UMaster_SiteChar.TC)
    UMaster_SiteChar.TC_std = std(UMaster_SiteChar.TC)

subplot(5,1,2);%PSAL
UMaster_SiteChar.PSAL = [UJuly_SiteChar.PSAL; UAug1_SiteChar.PSAL; UAug2_SiteChar.PSAL; USept_SiteChar.PSAL; UOct_SiteChar.PSAL];
hold on; box on;
plot(UMaster_SiteChar.SDN, UMaster_SiteChar.PSAL, 'b') 
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
% xlabel('Date');
ylabel('Salinity (ppt)');
title('Salinity');
UMaster_SiteChar.PSAL_median =  median(UMaster_SiteChar.PSAL)
    UMaster_SiteChar.PSAL_min = min(UMaster_SiteChar.PSAL)
    UMaster_SiteChar.PSAL_max = max(UMaster_SiteChar.PSAL)
UMaster_SiteChar.PSAL_mean = mean(UMaster_SiteChar.PSAL)
    UMaster_SiteChar.PSAL_std = std(UMaster_SiteChar.PSAL)

subplot(5,1,3);%PAR
UMaster_SiteChar.PAR = [UJuly_SiteChar.PAR; UAug1_SiteChar.PAR; UAug2_SiteChar.PAR; USept_SiteChar.PAR; UOct_SiteChar.PAR];
hold on; box on;
plot(UMaster_SiteChar.SDN, UMaster_SiteChar.PAR, 'color', [0.9290, 0.6940, 0.1250]) 
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
% xlabel('Date');
ylabel('PAR (umol photons m^-^2 s^-^1)');
title('Photosynthetically Active Light Radiation');
UMaster_SiteChar_PAR_day = find (UMaster_SiteChar.PAR>1)
UMaster_SiteChar.PAR_day = UMaster_SiteChar.PAR(UMaster_SiteChar_PAR_day);  
UMaster_SiteChar.PAR_day_median =  median(UMaster_SiteChar.PAR_day)
    UMaster_SiteChar.PAR_day_min = min(UMaster_SiteChar.PAR_day)
    UMaster_SiteChar.PAR_day_max = max(UMaster_SiteChar.PAR_day)
UMaster_SiteChar.PAR_day_mean = mean(UMaster_SiteChar.PAR_day)
    UMaster_SiteChar.PAR_day_std = std(UMaster_SiteChar.PAR_day)

subplot(5,1,4);%Pressure
UMaster_SiteChar.Pres = [UJuly_SiteChar.Pres; UAug1_SiteChar.Pres; UAug2_SiteChar.Pres; USept_SiteChar.Pres; UOct_SiteChar.Pres];
hold on; box on;
plot(UMaster_SiteChar.SDN, UMaster_SiteChar.Pres, 'color', [0.4940 0.1840 0.5560]) 
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
% xlabel('Date');
ylabel('Pressure (dBar)');
title('Pressure');
UMaster_SiteChar.Pres_median =  median(UMaster_SiteChar.Pres)
    UMaster_SiteChar.Pres_min = min(UMaster_SiteChar.Pres)
    UMaster_SiteChar.Pres_max = max(UMaster_SiteChar.Pres)
UMaster_SiteChar.Pres_mean = mean(UMaster_SiteChar.Pres)
    UMaster_SiteChar.Pres_std = std(UMaster_SiteChar.Pres)

subplot(5,1,5);%Velocity
UMaster_SiteChar.U0 = [UJuly_SiteChar.U0, UAug1_SiteChar.U0, USept_SiteChar.U0];
hold on; box on;
plot(UJuly_SiteChar.AD_SDN, UJuly_SiteChar.U0(1:24716), 'color', [0.4660 0.6740 0.1880])
plot(UAug1_SiteChar.AD_SDN, UAug1_SiteChar.U0, 'color', [0.4660 0.6740 0.1880])
plot(USept_SiteChar.AD_SDN, USept_SiteChar.U0, 'color', [0.4660 0.6740 0.1880])
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
xlabel('Date');
ylabel('Velocity (m s^-^1)');
title('Velocity at 1 m');
    U_U0_isnan = isnan(UMaster_SiteChar.U0);
    UMaster_SiteChar.U0_QC  = UMaster_SiteChar.U0;
    UMaster_SiteChar.U0_QC(U_U0_isnan)=[];
UMaster_SiteChar.U0_median = median (UMaster_SiteChar.U0_QC)
    UMaster_SiteChar.U0_min = min(UMaster_SiteChar.U0_QC)
    UMaster_SiteChar.U0_max = max(UMaster_SiteChar.U0_QC)
UMaster_SiteChar.U0_mean = mean (UMaster_SiteChar.U0_QC)
    UMaster_SiteChar.U0_std = std(UMaster_SiteChar.U0_QC)


%% Marker 32 Site Characterization Subplot
% Load full datasets
load('MJuly20_SiteChar.mat', 'MJuly_SiteChar')
load('MAug120_SiteChar.mat', 'MAug1_SiteChar' )
load('MAug220_SiteChar.mat', 'MAug2_SiteChar' )
load('MSept20_SiteChar.mat', 'MSept_SiteChar' )
load('MOct20_SiteChar.mat', 'MOct_SiteChar' )

% Marker 32 Master Datastructure  
MMaster_SiteChar.SDN = [MJuly_SiteChar.SDN; MAug1_SiteChar.SDN; MAug2_SiteChar.SDN; MSept_SiteChar.SDN; MOct_SiteChar.SDN];
MMaster_SiteChar.AD_SDN = [MJuly_SiteChar.AD_SDN; MAug1_SiteChar.AD_SDN; MSept_SiteChar.AD_SDN];
MMaster_SiteChar.AD_TC = [MJuly_SiteChar.AD_TC; MAug1_SiteChar.AD_TC; MSept_SiteChar.AD_TC];
MMaster_SiteChar.AD_Pres = [MJuly_SiteChar.AD_Pres; MAug1_SiteChar.AD_Pres; MSept_SiteChar.AD_Pres];

close all
figure
sgtitle('Marker 32 Reef Characterization Parameters')
subplot(5,1,1);%Temp
hold on; box on;
plot(MJuly_SiteChar.SDN, MJuly_SiteChar.TC, 'r')
plot(MAug1_SiteChar.AD_SDN, MAug1_SiteChar.AD_TC, 'r')
plot(MSept_SiteChar.AD_SDN, MSept_SiteChar.AD_TC, 'r')
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
% xlabel('Date');
ylabel('Temperature (C)');
title('Temperature');

MMaster_SiteChar.TC = [MJuly_SiteChar.TC; MAug1_SiteChar.AD_TC; MSept_SiteChar.AD_TC];
MMaster_SiteChar.TC_median =  median(MMaster_SiteChar.TC)
    MMaster_SiteChar.TC_min = min(MMaster_SiteChar.TC)
    MMaster_SiteChar.TC_max = max(MMaster_SiteChar.TC)
MMaster_SiteChar.TC_mean = mean(MMaster_SiteChar.TC)
    MMaster_SiteChar.TC_std = std(MMaster_SiteChar.TC)

subplot(5,1,2);%PSAL
hold on; box on;
plot(MJuly_SiteChar.SDN, MJuly_SiteChar.PSAL, 'b')
plot(MAug1_SiteChar.SDN, MAug1_SiteChar.PSAL, 'b')
plot(MAug2_SiteChar.SDN, MAug2_SiteChar.PSAL, 'b')
plot(MSept_SiteChar.SDN, MSept_SiteChar.PSAL, 'b')
plot(MOct_SiteChar.SDN, MOct_SiteChar.PSAL, 'b')
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
% xlabel('Date');
ylabel('Salinity (ppt)');
title('Salinity');

MMaster_SiteChar.PSAL = [MJuly_SiteChar.PSAL; MAug1_SiteChar.PSAL; MAug2_SiteChar.PSAL; MSept_SiteChar.PSAL; MOct_SiteChar.PSAL];
MMaster_SiteChar.PSAL_median =  median(MMaster_SiteChar.PSAL)
    MMaster_SiteChar.PSAL_min = min(MMaster_SiteChar.PSAL)
    MMaster_SiteChar.PSAL_max = max(MMaster_SiteChar.PSAL)
MMaster_SiteChar.PSAL_mean = mean(MMaster_SiteChar.PSAL)
    MMaster_SiteChar.PSAL_std = std(MMaster_SiteChar.PSAL)

subplot(5,1,3);%PAR
hold on; box on;
plot(MJuly_SiteChar.SDN, MJuly_SiteChar.PAR, 'color', [0.9290, 0.6940, 0.1250])
plot(MAug1_SiteChar.SDN, MAug1_SiteChar.PAR, 'color', [0.9290, 0.6940, 0.1250])
plot(MAug2_SiteChar.SDN, MAug2_SiteChar.PAR, 'color', [0.9290, 0.6940, 0.1250])
plot(MSept_SiteChar.SDN, MSept_SiteChar.PAR, 'color', [0.9290, 0.6940, 0.1250])
plot(MOct_SiteChar.SDN, MOct_SiteChar.PAR, 'color', [0.9290, 0.6940, 0.1250])
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
% xlabel('Date');
ylabel('PAR (umol photons m^-^2 s^-^1)');
title('Photosynthetically Active Light Radiation');

MMaster_SiteChar.PAR = [MJuly_SiteChar.PAR; MAug1_SiteChar.PAR; MAug2_SiteChar.PAR; MSept_SiteChar.PAR; MOct_SiteChar.PAR];
    MMaster_SiteChar_PAR_day = find (MMaster_SiteChar.PAR>1)
MMaster_SiteChar.PAR_day = MMaster_SiteChar.PAR(MMaster_SiteChar_PAR_day);  
MMaster_SiteChar.PAR_day_median =  median(MMaster_SiteChar.PAR_day)
    MMaster_SiteChar.PAR_day_min = min(MMaster_SiteChar.PAR_day)
    MMaster_SiteChar.PAR_day_max = max(MMaster_SiteChar.PAR_day)
MMaster_SiteChar.PAR_day_mean = mean(MMaster_SiteChar.PAR_day)
    MMaster_SiteChar.PAR_day_std = std(MMaster_SiteChar.PAR_day)

subplot(5,1,4);%Pressure
hold on; box on;
plot(MJuly_SiteChar.AD_SDN, MJuly_SiteChar.AD_Pres, 'color', [0.4940 0.1840 0.5560])
plot(MAug1_SiteChar.AD_SDN, MAug1_SiteChar.AD_Pres, 'color', [0.4940 0.1840 0.5560])
plot(MSept_SiteChar.AD_SDN, MSept_SiteChar.AD_Pres, 'color', [0.4940 0.1840 0.5560])
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
% xlabel('Date');
ylabel('Pressure (dBar)');
title('Pressure');

MMaster_SiteChar.Pres = [MJuly_SiteChar.AD_Pres; MAug1_SiteChar.AD_Pres; MSept_SiteChar.AD_Pres];
MMaster_SiteChar.Pres_median =  median(MMaster_SiteChar.Pres)
    MMaster_SiteChar.Pres_min = min(MMaster_SiteChar.Pres)
    MMaster_SiteChar.Pres_max = max(MMaster_SiteChar.Pres)
MMaster_SiteChar.Pres_mean = mean(MMaster_SiteChar.Pres)
    MMaster_SiteChar.Pres_std = std(MMaster_SiteChar.Pres)
    

subplot(5,1,5);%Velocity
hold on; box on;
plot(MJuly_SiteChar.AD_SDN, MJuly_SiteChar.U0, 'color', [0.4660 0.6740 0.1880])
plot(MAug1_SiteChar.AD_SDN, MAug1_SiteChar.U0, 'color', [0.4660 0.6740 0.1880])
plot(MSept_SiteChar.AD_SDN, MSept_SiteChar.U0, 'color', [0.4660 0.6740 0.1880])
set(gca, 'xlim', Xrange, 'XTick', tick, 'xticklabel', tick, 'XGrid', 'on');
datetick('x', 'mm/dd', 'keeplimits', 'keepticks');
xlabel('Date');
ylabel('Velocity (m s^-^1)');
title('Velocity at 1 m');

MMaster_SiteChar.U0 = [MJuly_SiteChar.U0, MAug1_SiteChar.U0, MSept_SiteChar.U0];
    M_U0_isnan = isnan(MMaster_SiteChar.U0);
    MMaster_SiteChar.U0_QC  = MMaster_SiteChar.U0;
    MMaster_SiteChar.U0_QC(M_U0_isnan)=[];
MMaster_SiteChar.U0_median = median (MMaster_SiteChar.U0_QC)
    MMaster_SiteChar.U0_min = min(MMaster_SiteChar.U0_QC)
    MMaster_SiteChar.U0_max = max(MMaster_SiteChar.U0_QC)
MMaster_SiteChar.U0_mean = mean (MMaster_SiteChar.U0_QC)
    MMaster_SiteChar.U0_std = std(MMaster_SiteChar.U0_QC)



%% Cudjoe Pre/Post Restoration Analysis

% Split Aug_2 into pre- and post- restoration - Restoration occured Monday Aug 17th 
% need UAug2_BMSbin.NEP_day_QC, UAug2_BMSbin.NEC_day_QC, UAug2_BMSbin.dDOXY_day_QC, UAug2_BMSbin.dTA_day_QC,    and day gradients
% datestr(UAug2_BMSbin.SDN_day)

Aug2_pre_end = find(UAug2_BMSbin.SDN_day==datenum('17-Aug-2020 11:00:00')); % - break at 132
Aug2_post_start = Aug2_pre_end+1;
Aug2_break = UAug2_BMSbin.SDN_day(1:Aug2_pre_end); 

%create new arrays for daytime data
Pre_UAug2_BMSbin.SDN_day = UAug2_BMSbin.SDN;
Pre_UAug2_BMSbin.PAR_day = UAug2_BMSbin.PAR(1,:);
Pre_UAug2_BMSbin.NEP_day_QC = UAug2_BMSbin.NEP_QC;
Pre_UAug2_BMSbin.NEC_day_QC = UAug2_BMSbin.NEC_QC;
Pre_UAug2_BMSbin.dDOXY_day_QC = UAug2_BMSbin.dDOXY_QC;      %DO Gradient
Pre_UAug2_BMSbin.dTA_day_QC = UAug2_BMSbin.dTA_QC;          %TA Gradient
Pre_UAug2_BMSbin.NEP_night_QC = UAug2_BMSbin.NEP_QC;
Pre_UAug2_BMSbin.NEC_night_QC = UAug2_BMSbin.NEC_QC;

%find all nightime datapoints - night = 1 
Pre_UAug2_BMSbin.inight = UAug2_BMSbin.PAR(1,:) < 10; 

%set all nightime values to NaN but don't remove
Pre_UAug2_BMSbin.SDN_day(Pre_UAug2_BMSbin.inight) = NaN;
Pre_UAug2_BMSbin.PAR_day (Pre_UAug2_BMSbin.inight) = NaN;
Pre_UAug2_BMSbin.NEP_day_QC(Pre_UAug2_BMSbin.inight) = NaN;
Pre_UAug2_BMSbin.NEC_day_QC(Pre_UAug2_BMSbin.inight) = NaN;
Pre_UAug2_BMSbin.dDOXY_day_QC(Pre_UAug2_BMSbin.inight) = NaN;      %DO Gradient
Pre_UAug2_BMSbin.dTA_day_QC(Pre_UAug2_BMSbin.inight) = NaN;        %TA Gradient

%break into pre- and post- restoration
Pre_UAug2_BMSbin.SDN_day_pre = Pre_UAug2_BMSbin.SDN_day(1:132);
Pre_UAug2_BMSbin.NEP_day_QC_pre = Pre_UAug2_BMSbin.NEP_day_QC(1:132);
Pre_UAug2_BMSbin.NEC_day_QC_pre = Pre_UAug2_BMSbin.NEC_day_QC(1:132);
Pre_UAug2_BMSbin.dDOXY_day_QC_pre = Pre_UAug2_BMSbin.dDOXY_day_QC(1:132);
Pre_UAug2_BMSbin.dTA_day_QC_pre = Pre_UAug2_BMSbin.dTA_day_QC(1:132);

Pre_UAug2_BMSbin.SDN_day_post = Pre_UAug2_BMSbin.SDN_day(133:end);
Pre_UAug2_BMSbin.NEP_day_QC_post = Pre_UAug2_BMSbin.NEP_day_QC(133:end);
Pre_UAug2_BMSbin.NEC_day_QC_post = Pre_UAug2_BMSbin.NEC_day_QC(133:end);
Pre_UAug2_BMSbin.dDOXY_day_QC_post = Pre_UAug2_BMSbin.dDOXY_day_QC(133:end);
Pre_UAug2_BMSbin.dTA_day_QC_post = Pre_UAug2_BMSbin.dTA_day_QC(133:end);

% remove NaNs from both pre and post datasets (isnan(UAug2_BMSbin.NEP_day_QC))=[];
Pre_UAug2_BMSbin.SDN_day_pre(isnan(Pre_UAug2_BMSbin.SDN_day_pre))=[];
Pre_UAug2_BMSbin.NEP_day_QC_pre(isnan(Pre_UAug2_BMSbin.NEP_day_QC_pre))=[];
Pre_UAug2_BMSbin.NEC_day_QC_pre(isnan(Pre_UAug2_BMSbin.NEC_day_QC_pre))=[];
Pre_UAug2_BMSbin.dDOXY_day_QC_pre(isnan(Pre_UAug2_BMSbin.dDOXY_day_QC_pre))=[];
Pre_UAug2_BMSbin.dTA_day_QC_pre(isnan(Pre_UAug2_BMSbin.dTA_day_QC_pre))=[];

Pre_UAug2_BMSbin.SDN_day_post(isnan(Pre_UAug2_BMSbin.SDN_day_post))=[];
Pre_UAug2_BMSbin.NEP_day_QC_post(isnan(Pre_UAug2_BMSbin.NEP_day_QC_post))=[];
Pre_UAug2_BMSbin.NEC_day_QC_post(isnan(Pre_UAug2_BMSbin.NEC_day_QC_post))=[];
Pre_UAug2_BMSbin.dDOXY_day_QC_post(isnan(Pre_UAug2_BMSbin.dDOXY_day_QC_post))=[];
Pre_UAug2_BMSbin.dTA_day_QC_post(isnan(Pre_UAug2_BMSbin.dTA_day_QC_post))=[];

% Make Master Pre and Post Restoration Datasets for ratios 
U_Pre.SDN  = [UJuly_BMSbin.SDN, UAug1_BMSbin.SDN, UAug2_BMSbin.SDN(1:Aug2_pre_end)];
U_Pre.NEP  = [UJuly_BMSbin.NEP_QC, UAug1_BMSbin.NEP_QC, UAug2_BMSbin.NEP_QC(1:Aug2_pre_end)];
U_Pre.NEC  = [UJuly_BMSbin.NEC_QC, UAug1_BMSbin.NEC_QC, UAug2_BMSbin.NEC_QC(1:Aug2_pre_end)];
U_Pre.dDOXY= [UJuly_BMSbin.dDOXY_QC, UAug1_BMSbin.dDOXY_QC, UAug2_BMSbin.dDOXY_QC(1:Aug2_pre_end)];
U_Pre.dTA  = [UJuly_BMSbin.dTA_QC, UAug1_BMSbin.dTA_QC, UAug2_BMSbin.dTA_QC(1:Aug2_pre_end)];

U_Post.SDN  = [UAug2_BMSbin.SDN(Aug2_post_start:end), USept_BMSbin.SDN, UOct_BMSbin.SDN];
U_Post.NEP  = [UAug2_BMSbin.NEP_QC(Aug2_post_start:end),USept_BMSbin.NEP_QC, UOct_BMSbin.NEP_QC];
U_Post.NEC  = [UAug2_BMSbin.NEC_QC(Aug2_post_start:end),USept_BMSbin.NEC_QC, UOct_BMSbin.NEC_QC];
U_Post.dDOXY= [UAug2_BMSbin.dDOXY_QC(Aug2_post_start:end),USept_BMSbin.dDOXY_QC, UOct_BMSbin.dDOXY_QC];
U_Post.dTA  = [UAug2_BMSbin.dTA_QC(Aug2_post_start:end),USept_BMSbin.dTA_QC, UOct_BMSbin.dTA_QC];

U_Pre.NEP_day = [UJuly_BMSbin.NEP_day_QC, UAug1_BMSbin.NEP_day_QC, Pre_UAug2_BMSbin.NEP_day_QC_pre];
U_Pre.NEC_day = [UJuly_BMSbin.NEC_day_QC, UAug1_BMSbin.NEC_day_QC, Pre_UAug2_BMSbin.NEC_day_QC_pre];
U_Pre.dDOXY_day = [UJuly_BMSbin.dDOXY_day_QC, UAug1_BMSbin.dDOXY_day_QC, Pre_UAug2_BMSbin.dDOXY_day_QC_pre];
U_Pre.dTA_day = [UJuly_BMSbin.dTA_day_QC, UAug1_BMSbin.dTA_day_QC, Pre_UAug2_BMSbin.dTA_day_QC_pre];

U_Post.NEP_day = [Pre_UAug2_BMSbin.NEP_day_QC_post,USept_BMSbin.NEP_day_QC, UOct_BMSbin.NEP_day_QC];
U_Post.NEC_day = [Pre_UAug2_BMSbin.NEC_day_QC_post,USept_BMSbin.NEC_day_QC, UOct_BMSbin.NEC_day_QC];
U_Post.dDOXY_day = [Pre_UAug2_BMSbin.dDOXY_day_QC_post,USept_BMSbin.dDOXY_day_QC, UOct_BMSbin.dDOXY_day_QC];
U_Post.dTA_day = [Pre_UAug2_BMSbin.dTA_day_QC_post,USept_BMSbin.dTA_day_QC, UOct_BMSbin.dTA_day_QC];

% make nighttime pre/post datasets
datestr(UAug2_BMSbin.SDN_night) 
find(UAug2_BMSbin.SDN_night==datenum('17-Aug-2020 06:00:00')) % - break at 55
% Aug2_night_pre = (1:55)
% Aug2_night_post = (56:end) 
U_Pre.SDN_night = [UJuly_BMSbin.SDN_night, UAug1_BMSbin.SDN_night, UAug2_BMSbin.SDN_night(1:55)];
U_Pre.NEP_night = [UJuly_BMSbin.NEP_night_QC, UAug1_BMSbin.NEP_night_QC, UAug2_BMSbin.NEP_night_QC(1:55)];
U_Pre.NEC_night = [UJuly_BMSbin.NEC_night_QC, UAug1_BMSbin.NEC_night_QC, UAug2_BMSbin.NEC_night_QC(1:55)];

U_Post.SDN_night = [UAug2_BMSbin.SDN_night(56:end),USept_BMSbin.SDN_night, UOct_BMSbin.SDN_night];
U_Post.NEP_night = [UAug2_BMSbin.NEP_night_QC(56:end),USept_BMSbin.NEP_night_QC, UOct_BMSbin.NEP_night_QC];
U_Post.NEC_night = [UAug2_BMSbin.NEC_night_QC(56:end),USept_BMSbin.NEC_night_QC, UOct_BMSbin.NEC_night_QC];

save('U_Pre_Post_ratios.mat', 'U_Pre', 'U_Post');

% load('U_Pre_Post_ratios.mat', 'U_Pre', 'U_Post');


%% Cudjoe Pre/Post GM Regressions
close all 
clc

Reg_X_range = (-30:30);

% pre_flux
[U_Pre.m_flux,U_Pre.b_flux,U_Pre.r_flux,U_Pre.sm_flux,U_Pre.sb_flux]=lsqfitgm(U_Pre.NEP_day,U_Pre.NEC_day);
U_Pre.Reg_Line_flux = U_Pre.m_flux*Reg_X_range + U_Pre.b_flux;
U_Pre.Ratio_flux = U_Pre.m_flux;
U_Pre.R2_flux = U_Pre.r_flux;
% plot
figure
hold on; box on;
plot(U_Pre.NEP_day,U_Pre.NEC_day,'o')
plot(Reg_X_range,U_Pre.Reg_Line_flux,'r')
xlabel('NEP');
ylabel('NEC');
title('Cudjoe Pre-Restoration NEC:NEP Ratio from Fluxes');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + U_Pre.Ratio_flux)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + U_Pre.R2_flux)

U_Pre.sm_flux

% % % pre_gradients 
% % % multiply o2 gradient by -1 for O2 production
% % U_Pre.dDOXY_day_Reg = -1.*U_Pre.dDOXY_day;
% % % divide TA data by 2 for alkalinity anomaly 
% % U_Pre.dTA_day_Reg = 0.5.*U_Pre.dTA_day;
% % 
% % [U_Pre.m_G,U_Pre.b_G,U_Pre.r_G,U_Pre.sm_G,U_Pre.sb_G]=lsqfitgm(U_Pre.dDOXY_day_Reg,U_Pre.dTA_day_Reg);
% % U_Pre.Reg_Line_G = U_Pre.m_G*Reg_X_range + U_Pre.b_G;
% % U_Pre.Ratio_G = U_Pre.m_G;
% % U_Pre.R2_G = U_Pre.r_G;
% % % plot
% % figure
% % hold on; box on;
% % plot(U_Pre.dDOXY_day_Reg,U_Pre.dTA_day_Reg,'o')
% % plot(Reg_X_range,U_Pre.Reg_Line_G,'r')
% % xlabel('NEP');
% % ylabel('NEC');
% % title('Cudjoe Pre-Restoration NEC:NEP Ratio from Gradients');
% % annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + U_Pre.Ratio_G)
% % annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + U_Pre.R2_G)

% post_flux
[U_Post.m_flux,U_Post.b_flux,U_Post.r_flux,U_Post.sm_flux,U_Post.sb_flux]=lsqfitgm(U_Post.NEP_day,U_Post.NEC_day);
U_Post.Reg_Line_flux = U_Post.m_flux*Reg_X_range + U_Post.b_flux;
U_Post.Ratio_flux = U_Post.m_flux;
U_Post.R2_flux = U_Post.r_flux;
% plot
figure
hold on; box on;
plot(U_Post.NEP_day,U_Post.NEC_day,'o')
plot(Reg_X_range,U_Post.Reg_Line_flux,'r')
xlabel('NEP');
ylabel('NEC');
title('Cudjoe Post-Restoration NEC:NEP Ratio from Fluxes');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + U_Post.Ratio_flux)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + U_Post.R2_flux)

% % % post_gradient
% % % multiply o2 gradient by -1 for O2 production
% % U_Post.dDOXY_day_Reg = -1.*U_Post.dDOXY_day;
% % % divide TA data by 2 for alkalinity anomaly 
% % U_Post.dTA_day_Reg = 0.5.*U_Post.dTA_day;
% % 
% % [U_Post.m_G,U_Post.b_G,U_Post.r_G,U_Post.sm_G,U_Post.sb_G]=lsqfitgm(U_Post.dDOXY_day_Reg,U_Post.dTA_day_Reg);
% % U_Post.Reg_Line_G = U_Post.m_G*Reg_X_range + U_Post.b_G;
% % U_Post.Ratio_G = U_Post.m_G;
% % U_Post.R2_G = U_Post.r_G;
% % % plot
% % figure
% % hold on; box on;
% % plot(U_Post.dDOXY_day_Reg,U_Post.dTA_day_Reg,'o')
% % plot(Reg_X_range,U_Post.Reg_Line_G,'r')
% % xlabel('NEP');
% % ylabel('NEC');
% % xlim([-3.5 3.5])
% % ylim([-3.5 3.5])
% % title('Cudjoe Post-Restoration NEC:NEP Ratio from Gradients');
% % annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + U_Post.Ratio_G)
% % annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + U_Post.R2_G)

%% Cudjoe overlapping ratios

% Ratios from Fluxes
close all
figure; % pre_flux
hold on; box on; 
plot(U_Pre.NEP_day,U_Pre.NEC_day,'r.', 'MarkerSize', 17) %pre = red
Pre_Reg = plot(Reg_X_range,U_Pre.Reg_Line_flux,'r');
plot(U_Post.NEP_day,U_Post.NEC_day,'bo', 'MarkerSize', 7, 'linewidth', 1.2)%post = blue
Post_Reg = plot(Reg_X_range,U_Post.Reg_Line_flux,'b');
xlabel('NEP');
ylabel('NEC');
legend([Pre_Reg Post_Reg], {'Pre-Restoration','Post-Restoration'}, 'location', 'northeast');
% set(gca,'Color',[0.9843    1.0000    0.9608])
title('Cudjoe NEC:NEP Ratios from Fluxes');
str1 = num2str(U_Pre.Ratio_flux,2);
str2 = num2str(U_Pre.R2_flux,2);
annotation('textbox', [0.131249999999996,0.892690677966101,0.166672568460813,0.046324529061482], 'String', "\color{red}Pre-NEC:NEP =" + str1, 'HorizontalAlignment', 'right', 'edgecolor', 'none', 'FitBoxToText', 'on')
annotation('textbox', [0.13138020833333,0.862023305084746,0.166672568460813,0.046324529061482], 'String', "\color{red}R^2 =" + str2, 'HorizontalAlignment', 'right', 'edgecolor', 'none', 'FitBoxToText', 'on')
str3 = num2str(U_Post.Ratio_flux,2);
str4 = num2str(U_Post.R2_flux,2);
annotation('textbox', [0.131249999999999,0.807245386192755,0.166672568460813,0.046324529061482], 'String', "\color{blue}Post-NEC:NEP =" + str3, 'HorizontalAlignment', 'right', 'edgecolor', 'none', 'FitBoxToText', 'on')
annotation('textbox', [0.131249999999999,0.793432203389831,0.166672568460813,0.046324529061482], 'String', "\color{blue}R^2 =" + str4, 'HorizontalAlignment', 'right', 'edgecolor', 'none', 'FitBoxToText', 'on')
%Calculate % increase in ratio
U_Ratio_increase_flux = ((U_Post.Ratio_flux-U_Pre.Ratio_flux)/U_Pre.Ratio_flux)*100
str5 = num2str(U_Ratio_increase_flux,2);
annotation('textbox', [0.451953124999998 0.890095338983051 0.130338541666668 0.0300000000000001], 'String', "Percent change in ratio=" + str5, 'FitBoxToText', 'on', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
ylim([-20 20])

% % Ratios from Gradients
% figure; 
% hold on; box on; 
% plot(U_Pre.dDOXY_day_Reg,U_Pre.dTA_day_Reg,'r.', 'MarkerSize', 17) %pre = red
% Pre_Reg = plot(Reg_X_range,U_Pre.Reg_Line_G,'r');
% plot(U_Post.dDOXY_day_Reg,U_Post.dTA_day_Reg,'bo', 'MarkerSize', 7, 'linewidth', 1.2)%post = blue
% Post_Reg = plot(Reg_X_range,U_Post.Reg_Line_G,'b');
% xlabel('dDO');
% ylabel('dTA');
% xlim([-4 4]);
% ylim([-4 4]);
% legend([Pre_Reg Post_Reg], {'Pre-Restoration Ratio','Post-Restoration Ratio'}, 'location', 'northeast');
% % set(gca,'Color',[0.9843    1.0000    0.9608])
% title('Cudjoe NEC:NEP Ratios from Gradients');
% str1 = num2str(U_Pre.Ratio_G,2);
% str2 = num2str(U_Pre.R2_G,2);
% annotation('textbox', [0.131249999999996,0.892690677966101,0.091927083333337,0.03], 'String', "\color{red}Pre-NEC:NEP =" + str1, 'HorizontalAlignment', 'right')
% annotation('textbox', [0.13138020833333,0.862023305084746,0.091927083333337,0.03], 'String', "\color{red}R^2 =" + str2, 'HorizontalAlignment', 'right')
% str3 = num2str(U_Post.Ratio_G,2);
% str4 = num2str(U_Post.R2_G,2);
% annotation('textbox', [0.131249999999999,0.823569915254237,0.091927083333337,0.03], 'String', "\color{blue}Post-NEC:NEP =" + str3, 'HorizontalAlignment', 'right')
% annotation('textbox', [0.131249999999999,0.793432203389831,0.091927083333337,0.03], 'String', "\color{blue}R^2 =" + str4, 'HorizontalAlignment', 'right')
% %Calculate % increase in ratio
% U_Ratio_increase_G = ((U_Post.Ratio_G-U_Pre.Ratio_G)/U_Pre.Ratio_G)*100
% % str5 = num2str(U_Ratio_increase_G,4);
% % annotation('textbox', [0.451953124999998 0.890095338983051 0.130338541666668 0.0300000000000001], 'String', "Percent change in ratio=" + str5, 'HorizontalAlignment', 'right')

%% M32 Pre/Post Restoration Analysis

% Make Master Pre and Post Restoration Datasets for ratios 
% M_Aug2 SP dies before restoration so don't have to split dataset

M_Pre.SDN  = [MJuly_BMSbin.SDN, MAug1_BMSbin.SDN, MAug2_BMSbin.SDN];
M_Pre.NEP  = [MJuly_BMSbin.NEP_QC, MAug1_BMSbin.NEP_QC, MAug2_BMSbin.NEP_QC];
M_Pre.NEC  = [MJuly_BMSbin.NEC_QC, MAug1_BMSbin.NEC_QC, MAug2_BMSbin.NEC_QC];
M_Pre.dDOXY= [MJuly_BMSbin.dDOXY_QC, MAug1_BMSbin.dDOXY_QC, MAug2_BMSbin.dDOXY_QC];
M_Pre.dTA  = [MJuly_BMSbin.dTA_QC, MAug1_BMSbin.dTA_QC, MAug2_BMSbin.dTA_QC];
M_Pre.SDN_day = [MJuly_BMSbin.SDN_day, MAug1_BMSbin.SDN_day, MAug2_BMSbin.SDN_day];
M_Pre.NEP_day = [MJuly_BMSbin.NEP_day_QC, MAug1_BMSbin.NEP_day_QC, MAug2_BMSbin.NEP_day_QC];
M_Pre.NEC_day = [MJuly_BMSbin.NEC_day_QC, MAug1_BMSbin.NEC_day_QC, MAug2_BMSbin.NEC_day_QC];
M_Pre.dDOXY_day = [MJuly_BMSbin.dDOXY_day_QC, MAug1_BMSbin.dDOXY_day_QC, MAug2_BMSbin.dDOXY_day_QC];
M_Pre.dTA_day = [MJuly_BMSbin.dTA_day_QC, MAug1_BMSbin.dTA_day_QC, MAug2_BMSbin.dTA_day_QC];
M_Pre.NEC_night = [MJuly_BMSbin.NEC_night_QC, MAug1_BMSbin.NEC_night_QC, MAug2_BMSbin.NEC_night_QC];
M_Pre.NEP_night = [MJuly_BMSbin.NEP_night_QC, MAug1_BMSbin.NEP_night_QC, MAug2_BMSbin.NEP_night_QC];
M_Pre.NEP_WM_day = [MJuly_BMSbin.NEP_WM_day_QC, MAug1_BMSbin.NEP_WM_day_QC, MAug2_BMSbin.NEP_WM_day_QC];
M_Pre.NEC_WM_day = [MJuly_BMSbin.NEC_WM_day_QC, MAug1_BMSbin.NEC_WM_day_QC, MAug2_BMSbin.NEC_WM_day_QC];

M_Post.SDN  = [MSept_BMSbin.SDN, MOct_BMSbin.SDN];
M_Post.NEP  = [MSept_BMSbin.NEP_QC, MOct_BMSbin.NEP_QC];
M_Post.NEC  = [MSept_BMSbin.NEC_QC, MOct_BMSbin.NEC_QC];
M_Post.dDOXY= [MSept_BMSbin.dDOXY_QC, MOct_BMSbin.dDOXY_QC];
M_Post.dTA  = [MSept_BMSbin.dTA_QC, MOct_BMSbin.dTA_QC];
M_Post.SDN_day = [MSept_BMSbin.SDN_day, MOct_BMSbin.SDN_day];
M_Post.NEP_day = [MSept_BMSbin.NEP_day_QC, MOct_BMSbin.NEP_day_QC];
M_Post.NEC_day = [MSept_BMSbin.NEC_day_QC, MOct_BMSbin.NEC_day_QC];
M_Post.dDOXY_day = [MSept_BMSbin.dDOXY_day_QC, MOct_BMSbin.dDOXY_day_QC];
M_Post.dTA_day = [MSept_BMSbin.dTA_day_QC, MOct_BMSbin.dTA_day_QC];
M_Post.NEC_night = [MSept_BMSbin.NEC_night_QC, MOct_BMSbin.NEC_night_QC];
M_Post.NEP_night = [MSept_BMSbin.NEP_night_QC, MOct_BMSbin.NEP_night_QC];
M_Post.NEP_WM_day = [MSept_BMSbin.NEP_WM_day_QC, MOct_BMSbin.NEP_WM_day_QC];
M_Post.NEC_WM_day = [MSept_BMSbin.NEC_WM_day_QC, MOct_BMSbin.NEC_WM_day_QC];

save('M_Pre_Post_ratios.mat', 'M_Pre', 'M_Post');

% load('M_Pre_Post_ratios.mat', 'M_Pre', 'M_Post');

%% M32 Pre/Post GM Regressions
close all 
clc

% pre_flux
[M_Pre.m_flux,M_Pre.b_flux,M_Pre.r_flux,M_Pre.sm_flux,M_Pre.sb_flux]=lsqfitgm(M_Pre.NEP_day,M_Pre.NEC_day);
M_Pre.Reg_Line_flux = M_Pre.m_flux*Reg_X_range + M_Pre.b_flux;
M_Pre.Ratio_flux = M_Pre.m_flux;
M_Pre.R2_flux = M_Pre.r_flux;
% plot
figure
hold on; box on;
plot(M_Pre.NEP_day,M_Pre.NEC_day,'o')
plot(Reg_X_range,M_Pre.Reg_Line_flux,'r')
xlabel('NEP');
ylabel('NEC');
title('Marker 32 Pre-Restoration NEC:NEP Ratio from Fluxes');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + M_Pre.Ratio_flux)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R^2 =" + M_Pre.R2_flux)

% % pre_gradients 
% % multiply o2 gradient by -1 for O2 production
% M_Pre.dDOXY_day_Reg = -1.*M_Pre.dDOXY_day;
% % divide TA data by 2 for alkalinity anomaly 
% M_Pre.dTA_day_Reg = 0.5.*M_Pre.dTA_day;
% 
% [M_Pre.m_G,M_Pre.b_G,M_Pre.r_G,M_Pre.sm_G,M_Pre.sb_G]=lsqfitgm(M_Pre.dDOXY_day_Reg,M_Pre.dTA_day_Reg);
% M_Pre.Reg_Line_G = M_Pre.m_G*Reg_X_range + M_Pre.b_G;
% M_Pre.Ratio_G = M_Pre.m_G;
% M_Pre.R2_G = M_Pre.r_G;
% % plot
% figure
% hold on; box on;
% plot(M_Pre.dDOXY_day_Reg,M_Pre.dTA_day_Reg,'o')
% plot(Reg_X_range,M_Pre.Reg_Line_G,'r')
% xlabel('dDO');
% ylabel('dTA');
% xlim([-5 5])
% ylim([-5 5])
% title('Marker 32 Pre-Restoration NEC:NEP Ratio from Gradients');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + M_Pre.Ratio_G)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R^2 =" + M_Pre.R2_G)

% post_flux
[M_Post.m_flux,M_Post.b_flux,M_Post.r_flux,M_Post.sm_flux,M_Post.sb_flux]=lsqfitgm(M_Post.NEP_day,M_Post.NEC_day);
M_Post.Reg_Line_flux = M_Post.m_flux*Reg_X_range + M_Post.b_flux;
M_Post.Ratio_flux = M_Post.m_flux;
M_Post.R2_flux = M_Post.r_flux;
% plot
figure
hold on; box on;
plot(M_Post.NEP_day,M_Post.NEC_day,'o')
plot(Reg_X_range,M_Post.Reg_Line_flux,'r')
xlabel('NEP');
ylabel('NEC');
title('Marker 32 Post-Restoration NEC:NEP Ratio from Fluxes');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + M_Post.Ratio_flux)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R^2 =" + M_Post.R2_flux)

% % post_gradient
% % multiply o2 gradient by -1 for O2 production
% M_Post.dDOXY_day_Reg = -1.*M_Post.dDOXY_day;
% % divide TA data by 2 for alkalinity anomaly 
% M_Post.dTA_day_Reg = 0.5.*M_Post.dTA_day;
% [M_Post.m_G,M_Post.b_G,M_Post.r_G,M_Post.sm_G,M_Post.sb_G]=lsqfitgm(M_Post.dDOXY_day_Reg,M_Post.dTA_day_Reg);
% M_Post.Reg_Line_G = M_Post.m_G*Reg_X_range + M_Post.b_G;
% M_Post.Ratio_G = M_Post.m_G;
% M_Post.R2_G = M_Post.r_G;
% % plot
% figure
% hold on; box on;
% plot(M_Post.dDOXY_day_Reg,M_Post.dTA_day_Reg,'o')
% plot(Reg_X_range,M_Post.Reg_Line_G,'r')
% xlabel('dDO');
% ylabel('dTA');
% title('Marker 32 Post-Restoration NEC:NEP Ratio from Gradients');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + M_Post.Ratio_G)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + M_Post.R2_G)

%% M32 overlapping ratios

% Ratios from Fluxes
% close all
figure; % pre_flux
hold on; box on; 
plot(M_Pre.NEP_day,M_Pre.NEC_day,'r.', 'MarkerSize', 17) %pre = red
Pre_Reg = plot(Reg_X_range,M_Pre.Reg_Line_flux,'r');
plot(M_Post.NEP_day,M_Post.NEC_day,'bo', 'MarkerSize', 7, 'linewidth', 1.2)%post = blue
Post_Reg = plot(Reg_X_range,M_Post.Reg_Line_flux,'b');
xlabel('NEP');
ylabel('NEC');
legend([Pre_Reg Post_Reg], {'Pre-Restoration Ratio','Post-Restoration Ratio'}, 'location', 'northeast');
% set(gca,'Color',[0.9843    1.0000    0.9608])
title('Marker 32 NEC:NEP Ratios from Fluxes');
str1 = num2str(M_Pre.Ratio_flux,2);
str2 = num2str(M_Pre.R2_flux,2);
annotation('textbox', [0.131249999999996,0.892690677966101,0.091927083333337,0.03], 'String', "\color{red}Pre-NEC:NEP =" + str1, 'HorizontalAlignment', 'right', 'edgecolor', 'none', 'FitBoxToText', 'on')
annotation('textbox', [0.13138020833333,0.862023305084746,0.091927083333337,0.03], 'String', "\color{red}R^2 =" + str2, 'HorizontalAlignment', 'right', 'edgecolor', 'none', 'FitBoxToText', 'on')
str3 = num2str(M_Post.Ratio_flux,2);
str4 = num2str(M_Post.R2_flux,2);
annotation('textbox', [0.131249999999999,0.823569915254237,0.091927083333337,0.03], 'String', "\color{blue}Post-NEC:NEP =" + str3, 'HorizontalAlignment', 'right', 'edgecolor', 'none', 'FitBoxToText', 'on')
annotation('textbox', [0.131249999999999,0.793432203389831,0.091927083333337,0.03], 'String', "\color{blue}R^2 =" + str4, 'HorizontalAlignment', 'right', 'edgecolor', 'none', 'FitBoxToText', 'on')
%Calculate % increase in ratio
M_Ratio_increase_flux = ((M_Post.Ratio_flux-M_Pre.Ratio_flux)/M_Pre.Ratio_flux)*100;
str5 = num2str(M_Ratio_increase_flux,2);
annotation('textbox', [0.451953124999998 0.890095338983051 0.130338541666668 0.0300000000000001], 'String', "Percent change in ratio=" + str5, 'FitBoxToText', 'on', 'HorizontalAlignment', 'right')

% % Ratios from Gradients
% figure; 
% hold on; box on; 
% plot(M_Pre.dDOXY_day_Reg,M_Pre.dTA_day_Reg,'r.', 'MarkerSize', 17) %pre = red
% Pre_Reg = plot(Reg_X_range,M_Pre.Reg_Line_G,'r');
% plot(M_Post.dDOXY_day_Reg,M_Post.dTA_day_Reg,'bo', 'MarkerSize', 7, 'linewidth', 1.2)%post = blue
% Post_Reg = plot(Reg_X_range,M_Post.Reg_Line_G,'b');
% xlabel('dDO');
% ylabel('dTA');
% xlim([-4 4]);
% ylim([-4 4]);
% legend([Pre_Reg Post_Reg], {'Pre-Restoration Ratio','Post-Restoration Ratio'}, 'location', 'northeast');
% % set(gca,'Color',[0.9843    1.0000    0.9608])
% title('Marker 32 NEC:NEP Ratios from Gradients');
% str1 = num2str(M_Pre.Ratio_G,2);
% str2 = num2str(M_Pre.R2_G,2);
% annotation('textbox', [0.131249999999996,0.892690677966101,0.091927083333337,0.03], 'String', "\color{red}Pre-NEC:NEP =" + str1, 'HorizontalAlignment', 'right')
% annotation('textbox', [0.13138020833333,0.862023305084746,0.091927083333337,0.03], 'String', "\color{red}R^2 =" + str2, 'HorizontalAlignment', 'right')
% str3 = num2str(M_Post.Ratio_G,2);
% str4 = num2str(M_Post.R2_G,2);
% annotation('textbox', [0.131249999999999,0.823569915254237,0.091927083333337,0.03], 'String', "\color{blue}Post-NEC:NEP =" + str3, 'HorizontalAlignment', 'right')
% annotation('textbox', [0.131249999999999,0.793432203389831,0.091927083333337,0.03], 'String', "\color{blue}R^2 =" + str4, 'HorizontalAlignment', 'right')
% M_Ratio_increase_G = ((M_Post.Ratio_G-M_Pre.Ratio_G)/M_Pre.Ratio_G)*100
% % str5 = num2str(M_Ratio_increase_G,4);
% % annotation('textbox', [0.451953124999998 0.890095338983051 0.130338541666668 0.0300000000000001], 'String', "Percent change in ratio=" + str5, 'HorizontalAlignment', 'right')
% M_Ratio_increase_G = ((M_Post.Ratio_G-M_Pre.Ratio_G)/M_Pre.Ratio_G)*100

%% Velocity Analysis
%Kruskal Wallis Test to see if velocities are significantly different
UMaster_SiteChar.U0 = [UJuly_SiteChar.U0, UAug1_SiteChar.U0, USept_SiteChar.U0];
MMaster_SiteChar.U0 = [MJuly_SiteChar.U0, MAug1_SiteChar.U0, MSept_SiteChar.U0];

%histograms
close all
figure
hold on; box on;
U_hist = histogram(UMaster_SiteChar.U0, 'BinWidth',0.002);
M_hist = histogram(MMaster_SiteChar.U0, 'BinWidth',0.002);
ylabel('Count');
xlabel('Velocity at 1m above bottom');
legend([U_hist M_hist], {'Cudjoe','Marker 32'}, 'location', 'northeast');
title('Velocity Histogram');

%Normplots 
figure
hold on; box on;
U_norm = normplot(UMaster_SiteChar.U0);
title('Cudjoe Normplot');

figure
hold on; box on;
M_norm = normplot(MMaster_SiteChar.U0);
title('M32 Normplot');

%Shapiro Wilk normality test: null hypothesis of composite normality
%       H = 0 => Do not reject the null hypothesis at significance level ALPHA.
%       H = 1 => Reject the null hypothesis at significance level ALPHA.
% clc
% [H, pValue, SWstatistic] = swtest(UMaster_SiteChar.U0, 0.05)
% [H, pValue, SWstatistic] = swtest(MMaster_SiteChar.U0, 0.05)

%Lilliefors Normal Distribution Test - returns h = 0, indicating that the null hypothesis 
%-- that the data is a sample from a normal distribution -- is not rejected.
% [h,p] = lillietest(UMaster_SiteChar.U0)
% [h,p] = lillietest(MMaster_SiteChar.U0)

% *************** both test results confirm non-normal datasets ************

%Kruskal-Wallis test - returns the p-value for the null hypothesis that the 
    %data in each column of the matrix x comes from the same distribution, using 
    %a Kruskal-Wallis test. The alternative hypothesis is that not all samples 
    %come from the same distribution.
        %an extension of the Wilcoxon rank sum test to more than two groups, the 
        %Kruskal-Wallis test is valid for data that has two or more groups. It 
        %compares the medians of the groups of data in x to determine if the 
        %samples come from the same population (or, equivalently, from different
        %populations with the same distribution).
%if p-value is less than 0.05, reject the null hypothesis that the data comes from the same distribution; 

U_U0_length = length (UMaster_SiteChar.U0);
M_U0_length = length (MMaster_SiteChar.U0);
close all 
clc
%create vector of all velocity values- append M32 to Cudjoe  
Master_Velocity = [UMaster_SiteChar.U0 MMaster_SiteChar.U0];
dataset1 = zeros(1,U_U0_length); %Zero indicates group U 
dataset2 = ones(1,M_U0_length); %one indicates group M32 
%create a second vector indicating the dataset from which the corresponding measurement is made
dataset = [dataset1 dataset2];
Group = cell(size(dataset));
Group(dataset==0)={'Cudjoe'};
Group(dataset==1)={'Marker 32'};
%p = kruskalwallis(x)
pVelocity = kruskalwallis(Master_Velocity,Group);


% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(UMaster_SiteChar.U0,MMaster_SiteChar.U0) 

% remove NaNs 
U_U0_isnan = isnan(UMaster_SiteChar.U0)
UMaster_SiteChar.U0_QC  = UMaster_SiteChar.U0
UMaster_SiteChar.U0_QC(U_U0_isnan)=[]

M_U0_isnan = isnan(MMaster_SiteChar.U0)
MMaster_SiteChar.U0_QC  = MMaster_SiteChar.U0
MMaster_SiteChar.U0_QC(M_U0_isnan)=[]

% calculate median velocities
UMaster_SiteChar.U0_median = median (UMaster_SiteChar.U0_QC)
    UMaster_SiteChar.U0_min = min(UMaster_SiteChar.U0_QC)
    UMaster_SiteChar.U0_max = max(UMaster_SiteChar.U0_QC)

MMaster_SiteChar.U0_median = median (MMaster_SiteChar.U0_QC)
    MMaster_SiteChar.U0_min = min(MMaster_SiteChar.U0_QC)
    MMaster_SiteChar.U0_max = max(MMaster_SiteChar.U0_QC)
    
% calculate mean velocities
UMaster_SiteChar.U0_mean = mean (UMaster_SiteChar.U0_QC)
    UMaster_SiteChar.U0_std = std(UMaster_SiteChar.U0_QC)
MMaster_SiteChar.U0_mean = mean (MMaster_SiteChar.U0_QC)
    MMaster_SiteChar.U0_std = std(MMaster_SiteChar.U0_QC)

%histograms
close all
figure
hold on; box on;
U_hist = histogram(UMaster_SiteChar.U0_QC);
M_hist = histogram(MMaster_SiteChar.U0_QC);
ylabel('Count');
xlabel('Velocity at 1m above bottom');
legend([U_hist M_hist], {'Cudjoe','Marker 32'}, 'location', 'northeast');
title('Velocity Histogram');

U_U0_length_QC = length (UMaster_SiteChar.U0_QC)
M_U0_length_QC = length (MMaster_SiteChar.U0_QC)
Master_Velocity_QC = [UMaster_SiteChar.U0_QC MMaster_SiteChar.U0_QC];

dataset1 = zeros(1,U_U0_length_QC); %Zero indicates group U 
dataset2 = ones(1,M_U0_length_QC); %one indicates group M32 
%create a second vector indicating the dataset from which the corresponding measurement is made
dataset = [dataset1 dataset2];
Group = cell(size(dataset));
Group(dataset==0)={'Cudjoe'};
Group(dataset==1)={'Marker 32'};

%p = kruskalwallis(x)
pVelocity = kruskalwallis(Master_Velocity_QC,Group)

%Boxplot
close all
figure
hold on; box on;
boxplot(Master_Velocity_QC,Group)
ylabel('Velocity at 1 m')
title('Velocity Boxplots Full Datasets')
annotation('textbox', [0.763151041666666,0.883474576271186,0.139583333333333,0.035699152542373], 'String', "Kruskal Wallis Test p-value =" + pVelocity)

%% Cudjoe Diel plots 

% Cudjoe
U_necdbin_Pre = parse_to_diel(U_Pre.SDN, U_Pre.NEC, 24);
U_necdbin_Post = parse_to_diel(U_Post.SDN, U_Post.NEC, 24);
U_necdbin_Pre_median = nanmedian(U_necdbin_Pre,1);
U_necdbin_Post_median = nanmedian(U_necdbin_Post,1);

U_nepdbin_Pre = parse_to_diel(U_Pre.SDN, U_Pre.NEP, 24);
U_nepdbin_Post = parse_to_diel(U_Post.SDN, U_Post.NEP, 24);
U_nepdbin_Pre_median = nanmedian(U_nepdbin_Pre,1);
U_nepdbin_Post_median = nanmedian(U_nepdbin_Post,1);

% U Diel Pre Plot
close all
figure
subplot(1,2,1)
hold on; box on;
% shade nighttime hours
    Night1x = [1; 7; 7; 1];
    Night1y = [-29.9;-29.9;30;30];
    patch(Night1x,Night1y,[1 1 1]*0.95, 'LineStyle','none')
    Night2x = [21; 24; 24; 21];
    Night2y = [-29.9;-29.9;30;30];
    patch(Night2x,Night2y,[1 1 1]*0.95, 'LineStyle','none')
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, U_nepdbin_Pre, 'bo', 'markersize', 3);
plot(1:24, U_necdbin_Pre, 'ro', 'markersize', 3);
plot(1:24, nanmedian(U_nepdbin_Pre,1), 'bo-');
plot(1:24, nanmedian(U_necdbin_Pre,1), 'ro-')
ylabel(['\color{blue}NEP \color{black}or \color{red}NEC']);
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
set(gca, 'Layer','top','XGrid', 'on', 'XColor', 'k');
xlim([1 24]);
ylim([-30 30]);
title('Cudjoe Pre-Restoration Diel Plot');
xlabel('hour of day');

% U Diel post plot 
subplot(1,2,2)
hold on; box on;
% shade nighttime hours
    Night1x = [1; 7; 7; 1];
    Night1y = [-29.9;-29.9;30;30];
    patch(Night1x,Night1y,[1 1 1]*0.95, 'LineStyle','none')
    Night2x = [21; 24; 24; 21];
    Night2y = [-29.9;-29.9;30;30];
    patch(Night2x,Night2y,[1 1 1]*0.95, 'LineStyle','none')
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, U_nepdbin_Post, 'bo', 'markersize', 3);
plot(1:24, U_necdbin_Post, 'ro', 'markersize', 3);
plot(1:24, nanmedian(U_nepdbin_Post,1), 'bo-');
plot(1:24, nanmedian(U_necdbin_Post,1), 'ro-')
ylabel(['\color{blue}NEP \color{black}or \color{red}NEC']);
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
set(gca,'Layer','top', 'XGrid', 'on');
xlim([1 24]);
ylim([-30 30]);
title('Cudjoe Post-Restoration Diel Plot');
xlabel('hour of day');

%integrate to find daily cumulative flux
U_pre_net_NEP = trapz(1:24,U_nepdbin_Pre_median)
U_pre_net_NEC = trapz(1:24,U_necdbin_Pre_median)
U_post_net_NEP = trapz(1:24,U_nepdbin_Post_median)
U_post_net_NEC = trapz(1:24,U_necdbin_Post_median)

% Change in cumulative NEP 
U_net_NEP_change = ((U_post_net_NEP - U_pre_net_NEP)/U_pre_net_NEP)*100

% Change in cumulative NEC
U_net_NEC_change = ((U_post_net_NEC - U_pre_net_NEC)/U_pre_net_NEC)*100

% kruskalwallis test for significance of net NEC change 
U_pre_NEC_length = length (U_Pre.NEC);
U_post_NEC_length = length (U_Post.NEC);
Master_U_NEC = [U_Pre.NEC U_Post.NEC];
dataset1 = zeros(1,U_pre_NEC_length); %Zero indicates group U 
dataset2 = ones(1,U_post_NEC_length); %one indicates group M32 
%create a second vector indicating the dataset from which the corresponding measurement is made
dataset = [dataset1 dataset2];
Group = cell(size(dataset));
Group(dataset==0)={'Pre-restoration'};
Group(dataset==1)={'Post-restoration'};
%p = kruskalwallis(x)
[U_NEC_p,tbl,stats] = kruskalwallis(Master_U_NEC,Group)

% kruskalwallis test for significance of net NEP change 
U_pre_NEP_length = length (U_Pre.NEP);
U_post_NEP_length = length (U_Post.NEP);
Master_U_NEP = [U_Pre.NEP U_Post.NEP];
dataset1 = zeros(1,U_pre_NEP_length); %Zero indicates group U 
dataset2 = ones(1,U_post_NEP_length); %one indicates group M32 
%create a second vector indicating the dataset from which the corresponding measurement is made
dataset = [dataset1 dataset2];
Group = cell(size(dataset));
Group(dataset==0)={'Pre-restoration'};
Group(dataset==1)={'Post-restoration'};
%p = kruskalwallis(x)
[U_NEP_p,tbl,stats] = kruskalwallis(Master_U_NEP,Group)

% calculate net nighttime calcification and percent change  
U_pre_net_NEC_night = trapz(1:6,U_necdbin_Pre_median(1:6))+trapz(21:24,U_necdbin_Pre_median(21:24))
U_post_net_NEC_night = trapz(1:6,U_necdbin_Post_median(1:6))+trapz(21:24,U_necdbin_Post_median(21:24))
U_percent_change_NEC_night = ((U_post_net_NEC_night-U_pre_net_NEC_night)/U_pre_net_NEC_night)*100

% calculate net daytime calcification and percent change  
U_pre_net_NEC_day = trapz(6:21,U_necdbin_Pre_median(6:21))
U_post_net_NEC_day = trapz(6:21,U_necdbin_Post_median(6:21))
U_percent_change_NEC_day = ((U_post_net_NEC_day-U_pre_net_NEC_day)/U_pre_net_NEC_day)*100

% calculate net daytime production and percent change 
U_pre_net_NEP_day = trapz(6:21,U_nepdbin_Pre_median(6:21))
U_post_net_NEP_day = trapz(6:21,U_nepdbin_Post_median(6:21))
U_percent_change_NEP_day = ((U_post_net_NEP_day-U_pre_net_NEP_day)/U_pre_net_NEP_day)*100

% calculate net nighttime production and percent change  
U_pre_net_NEP_night = trapz(1:6,U_nepdbin_Pre_median(1:6))+trapz(21:24,U_nepdbin_Pre_median(21:24))
U_post_net_NEP_night = trapz(1:6,U_nepdbin_Post_median(1:6))+trapz(21:24,U_nepdbin_Post_median(21:24))
U_percent_change_NEP_night = ((U_post_net_NEP_night-U_pre_net_NEP_night)/U_pre_net_NEP_night)*100

%% Marker 32 Diel Plots 
M_necdbin_Pre = parse_to_diel(M_Pre.SDN, M_Pre.NEC, 24);
M_necdbin_Post = parse_to_diel(M_Post.SDN, M_Post.NEC, 24);
M_necdbin_Pre_median = nanmedian(M_necdbin_Pre,1);
M_necdbin_Post_median = nanmedian(M_necdbin_Post,1);

M_nepdbin_Pre = parse_to_diel(M_Pre.SDN, M_Pre.NEP, 24);
M_nepdbin_Post = parse_to_diel(M_Post.SDN, M_Post.NEP, 24);
M_nepdbin_Pre_median = nanmedian(M_nepdbin_Pre,1);
M_nepdbin_Post_median = nanmedian(M_nepdbin_Post,1);

% M Diel Pre Plot  
%close all
figure
subplot(1,2,1)
hold on; box on;
% shade nighttime hours
    Night1x = [1; 7; 7; 1];
    Night1y = [-29.9;-29.9;30;30];
    patch(Night1x,Night1y,[1 1 1]*0.95, 'LineStyle','none')
    Night2x = [21; 24; 24; 21];
    Night2y = [-29.9;-29.9;30;30];
    patch(Night2x,Night2y,[1 1 1]*0.95, 'LineStyle','none')
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, M_nepdbin_Pre, 'bo', 'markersize', 3);
plot(1:24, M_necdbin_Pre, 'ro', 'markersize', 3);
plot(1:24, nanmedian(M_nepdbin_Pre,1), 'bo-');
plot(1:24, nanmedian(M_necdbin_Pre,1), 'ro-')
ylabel(['\color{blue}NEP \color{black}or \color{red}NEC']);
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
set(gca, 'Layer','top','XGrid', 'on', 'XColor', 'k');
xlim([1 24]);
ylim([-30 30]);
title('Marker 32 Pre-Restoration Diel Plot');
xlabel('hour of day');

% M Diel post plot 
subplot (1,2,2)
hold on; box on;
% shade nighttime hours
    Night1x = [1; 7; 7; 1];
    Night1y = [-29.9;-29.9;30;30];
    patch(Night1x,Night1y,[1 1 1]*0.95, 'LineStyle','none')
    Night2x = [21; 24; 24; 21];
    Night2y = [-29.9;-29.9;30;30];
    patch(Night2x,Night2y,[1 1 1]*0.95, 'LineStyle','none')
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, M_nepdbin_Post, 'bo', 'markersize', 3);
plot(1:24, M_necdbin_Post, 'ro', 'markersize', 3);
plot(1:24, nanmedian(M_nepdbin_Post,1), 'bo-');
plot(1:24, nanmedian(M_necdbin_Post,1), 'ro-')
ylabel(['\color{blue}NEP \color{black}or \color{red}NEC']);
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
set(gca,'Layer','top', 'XGrid', 'on');
xlim([1 24]);
ylim([-30 30]);
title('Marker 32 Post-Restoration Diel Plot');
xlabel('hour of day');

M_pre_net_NEP = trapz(1:24,M_nepdbin_Pre_median)
M_pre_net_NEC = trapz(1:24,M_necdbin_Pre_median)
M_post_net_NEP = trapz(1:24,M_nepdbin_Post_median)
M_post_net_NEC = trapz(1:24,M_necdbin_Post_median)

% Change in cumulative NEP 
M_net_NEP_change = ((M_post_net_NEP - M_pre_net_NEP)/M_pre_net_NEP)*100

% Change in cumulative NEC
M_net_NEC_change = ((M_post_net_NEC - M_pre_net_NEC)/M_pre_net_NEC)*100

% kruskalwallis Test significance of change in net NEC
M_pre_NEC_length = length (M_Pre.NEC)
M_post_NEC_length = length (M_Post.NEC)
Master_M_NEC = [M_Pre.NEC M_Post.NEC];
length(Master_M_NEC) 
dataset1 = zeros(1,M_pre_NEC_length); %Zero indicates group U 
dataset2 = ones(1,M_post_NEC_length); %one indicates group M32 
%create a second vector indicating the dataset from which the corresponding measurement is made
dataset = [dataset1 dataset2];
Group = cell(size(dataset));
Group(dataset==0)={'Pre-restoration'};
Group(dataset==1)={'Post-restoration'};
%p = kruskalwallis(x)
[M_NEC_p,tbl,stats] = kruskalwallis(Master_M_NEC,Group)

% kruskalwallis Test significance of change in net NEP
M_pre_NEP_length = length (M_Pre.NEP)
M_post_NEP_length = length (M_Post.NEP)
Master_M_NEP = [M_Pre.NEP M_Post.NEP];
length(Master_M_NEP) 
dataset1 = zeros(1,M_pre_NEP_length); %Zero indicates group U 
dataset2 = ones(1,M_post_NEP_length); %one indicates group M32 
%create a second vector indicating the dataset from which the corresponding measurement is made
dataset = [dataset1 dataset2];
Group = cell(size(dataset));
Group(dataset==0)={'Pre-restoration'};
Group(dataset==1)={'Post-restoration'};
%p = kruskalwallis(x)
[M_NEP_p,tbl,stats] = kruskalwallis(Master_M_NEP,Group)


% calculate net nighttime calcification and percent change  
M_pre_net_NEC_night = trapz(1:6,M_necdbin_Pre_median(1:6))+trapz(21:24,M_necdbin_Pre_median(21:24))
M_post_net_NEC_night = trapz(1:6,M_necdbin_Post_median(1:6))+trapz(21:24,M_necdbin_Post_median(21:24))
M_percent_change_NEC_night = ((M_post_net_NEC_night-M_pre_net_NEC_night)/M_pre_net_NEC_night)*100

% calculate net daytime calcification and percent change  
M_pre_net_NEC_day = trapz(6:21,M_necdbin_Pre_median(6:21))
M_post_net_NEC_day = trapz(6:21,M_necdbin_Post_median(6:21))
M_percent_change_NEC_day = ((M_post_net_NEC_day-M_pre_net_NEC_day)/M_pre_net_NEC_day)*100

% calculate net daytime production and percent change 
M_pre_net_NEP_day = trapz(6:21,M_nepdbin_Pre_median(6:21))
M_post_net_NEP_day = trapz(6:21,M_nepdbin_Post_median(6:21))
M_percent_change_NEP_day = ((M_post_net_NEP_day-M_pre_net_NEP_day)/M_pre_net_NEP_day)*100

% calculate net nighttime production and percent change  
M_pre_net_NEP_night = trapz(1:6,M_nepdbin_Pre_median(1:6))+trapz(21:24,M_nepdbin_Pre_median(21:24))
M_post_net_NEP_night = trapz(1:6,M_nepdbin_Post_median(1:6))+trapz(21:24,M_nepdbin_Post_median(21:24))
M_percent_change_NEP_night = ((M_post_net_NEP_night-M_pre_net_NEP_night)/M_pre_net_NEP_night)*100

%% Nighttime Bar Graphs
% Cudjoe
% want to consider night: 21:00 - 6:00am --> 10 hours for both pre and post
% U_necdbin_Post has 24 columns, each column = 1 hour 
U_necdbin_Pre = parse_to_diel(U_Pre.SDN, U_Pre.NEC, 24);
U_necdbin_Post = parse_to_diel(U_Post.SDN, U_Post.NEC, 24);
U_necdbin_Pre_median = nanmedian(U_necdbin_Pre,1)
U_necdbin_Post_median = nanmedian(U_necdbin_Post,1)

x = [-3:6];
U_val1 = [U_necdbin_Pre_median(21:24)'; U_necdbin_Pre_median(1:6)']
U_val2 = [U_necdbin_Post_median(21:24)'; U_necdbin_Post_median(1:6)']
U_vals = [U_val1 U_val2]

% Marker 32
% want to consider night: 21:00 - 6:00am --> 10 hours for both pre and post
% M_necdbin_Post has 24 columns, each column = 1 hour 
M_necdbin_Pre = parse_to_diel(M_Pre.SDN, M_Pre.NEC, 24);
M_necdbin_Post = parse_to_diel(M_Post.SDN, M_Post.NEC, 24);
M_necdbin_Pre_median = nanmedian(M_necdbin_Pre,1)
M_necdbin_Post_median = nanmedian(M_necdbin_Post,1)

x = [-3:6];
M_val1 = [M_necdbin_Pre_median(21:24)'; M_necdbin_Pre_median(1:6)']
M_val2 = [M_necdbin_Post_median(21:24)'; M_necdbin_Post_median(1:6)']
M_vals = [M_val1 M_val2]

% Plot Bar graph subplot
close all
subplot(1, 2, 1)
b = bar(x,U_vals);
xticks([-3:6])
xticklabels({'21:00','22:00','23:00','24:00','01:00','02:00','03:00','04:00','05:00','06:00',})
xlabel('Nighttime Hours');
ylabel('NEC (CaCO_3 m^-^2 hr^-^1)');
legend('Pre-Restoration', 'Post-Restoration', 'location', 'southeast')
title('Cudjoe Median Nighttime NEC');
ylim([-4.5 2.5]);

subplot(1, 2, 2)
b = bar(x,M_vals);
xticks([-3:6])
xticklabels({'21:00','22:00','23:00','24:00','01:00','02:00','03:00','04:00','05:00','06:00',})
xlabel('Nighttime Hours');
ylabel('NEC (CaCO_3 m^-^2 hr^-^1)');
legend('Pre-Restoration', 'Post-Restoration', 'location', 'southeast')
title('Marker 32 Median Nighttime NEC');
ylim([-4.5 2.5]);

%Marker 32
figure 
subplot(1, 2, 1)
boxplot(M_Pre.NEC_night)
ylim([-14 14])
title('Marker 32 Pre-Restoration Nighttime NEC')
subplot(1, 2, 2)
boxplot(M_Post.NEC_night)
ylim([-14 14])
title('Marker 32 Post-Restoration Nighttime NEC')


%% Kruzkal Wallis test between pre/post nighttime data to see if significantly different
% test all nighttime NEC values 
% test if data in each categorical group is from the same distribution
% Cudjoe night NEC
clc 
close all
U_dataset1 = zeros(1,length(U_Pre.NEC_night)); %Zero indicates group Pre
U_dataset2 = ones(1,length(U_Post.NEC_night)); %one indicates group Post
%create a second vector indicating the dataset from which the corresponding measurement is made
U_dataset = [U_dataset1 U_dataset2];
U_Group = cell(size(U_dataset)); % one row of data 
U_Group(U_dataset==0)={'Pre-Restoration'};
U_Group(U_dataset==1)={'Post-Restoration'};

U_length_Pre_NEC_night = length (U_Pre.NEC_night);
U_length_Post_NEC_night = length (U_Post.NEC_night);
U_Master_NEC_night = [U_Pre.NEC_night U_Post.NEC_night];% one row of data 
length(U_Master_NEC_night)==U_length_Pre_NEC_night+U_length_Post_NEC_night % 1= true

%p = kruskalwallis(x, group) - 
[U_NEC_night_p,tbl,stats] = kruskalwallis(U_Master_NEC_night, U_Group)

% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(U_Pre.NEC_night,U_Post.NEC_night) 

% Cudjoe night NEP
clc 
close all
U_dataset1 = zeros(1,length(U_Pre.NEP_night)); %Zero indicates group Pre
U_dataset2 = ones(1,length(U_Post.NEP_night)); %one indicates group Post
%create a second vector indicating the dataset from which the corresponding measurement is made
U_dataset = [U_dataset1 U_dataset2];
U_Group = cell(size(U_dataset)); % one row of data 
U_Group(U_dataset==0)={'Pre-Restoration'};
U_Group(U_dataset==1)={'Post-Restoration'};

U_length_Pre_NEP_night = length (U_Pre.NEP_night);
U_length_Post_NEP_night = length (U_Post.NEP_night);
U_Master_NEP_night = [U_Pre.NEP_night U_Post.NEP_night];% one row of data 
length(U_Master_NEP_night)==U_length_Pre_NEP_night+U_length_Post_NEP_night % 1= true

%p = kruskalwallis(x, group) - 
[U_NEP_night_p,tbl,stats] = kruskalwallis(U_Master_NEP_night, U_Group)

% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(U_Pre.NEP_night,U_Post.NEP_night) 



% Marker 32 NEC Night
clc 
close all
M_length_Pre_NEC_night = length (M_Pre.NEC_night);
M_length_Post_NEC_night = length (M_Post.NEC_night);
M_Master_NEC_night = [M_Pre.NEC_night M_Post.NEC_night];% one row of data 
length(M_Master_NEC_night)==M_length_Pre_NEC_night+M_length_Post_NEC_night % 1= true

M_dataset1 = zeros(1,length(M_Pre.NEC_night)); %Zero indicates group Pre
M_dataset2 = ones(1,length(M_Post.NEC_night)); %one indicates group Post
%create a second vector indicating the dataset from which the corresponding measurement is made
M_dataset = [M_dataset1 M_dataset2];
M_Group = cell(size(M_dataset)); % one row of data 
M_Group(M_dataset==0)={'Pre-Restoration'};
M_Group(M_dataset==1)={'Post-Restoration'};

%p = kruskalwallis(x, group) - test if data in each categorical group is from the
%same distribution
[M_NEC_night_p,tbl,stats] = kruskalwallis(M_Master_NEC_night, M_Group)

% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(M_Pre.NEC_night,M_Post.NEC_night)


% Marker 32 NEP Night
clc 
close all
M_length_Pre_NEP_night = length (M_Pre.NEP_night);
M_length_Post_NEP_night = length (M_Post.NEP_night);
M_Master_NEP_night = [M_Pre.NEP_night M_Post.NEP_night];% one row of data 
length(M_Master_NEP_night)==M_length_Pre_NEP_night+M_length_Post_NEP_night % 1= true

M_dataset1 = zeros(1,length(M_Pre.NEP_night)); %Zero indicates group Pre
M_dataset2 = ones(1,length(M_Post.NEP_night)); %one indicates group Post
%create a second vector indicating the dataset from which the corresponding measurement is made
M_dataset = [M_dataset1 M_dataset2];
M_Group = cell(size(M_dataset)); % one row of data 
M_Group(M_dataset==0)={'Pre-Restoration'};
M_Group(M_dataset==1)={'Post-Restoration'};

%p = kruskalwallis(x, group) - test if data in each categorical group is from the
%same distribution
[M_NEP_night_p,tbl,stats] = kruskalwallis(M_Master_NEP_night, M_Group)

% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(M_Pre.NEP_night,M_Post.NEP_night)

%% Kruzkal Wallis test between pre/post daytime data to see if significantly different
% test all daytime NEC values 
% test if data in each categorical group is from the same distribution
% Cudjoe NEC
clc 
close all
U_dataset1 = zeros(1,length(U_Pre.NEC_day)); %Zero indicates group Pre
U_dataset2 = ones(1,length(U_Post.NEC_day)); %one indicates group Post
%create a second vector indicating the dataset from which the corresponding measurement is made
U_dataset = [U_dataset1 U_dataset2];
U_Group = cell(size(U_dataset)); % one row of data 
U_Group(U_dataset==0)={'Pre-Restoration'};
U_Group(U_dataset==1)={'Post-Restoration'};

U_length_Pre_NEC_day = length (U_Pre.NEC_day);
U_length_Post_NEC_day = length (U_Post.NEC_day);
U_Master_NEC_day = [U_Pre.NEC_day U_Post.NEC_day];% one row of data 
length(U_Master_NEC_day)==U_length_Pre_NEC_day+U_length_Post_NEC_day % 1= true

%p = kruskalwallis(x, group) - 
[U_NEC_day_p,tbl,stats] = kruskalwallis(U_Master_NEC_day, U_Group)

% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(U_Pre.NEC_day,U_Post.NEC_day) 


% Cudjoe NEP
clc 
close all
U_dataset1 = zeros(1,length(U_Pre.NEP_day)); %Zero indicates group Pre
U_dataset2 = ones(1,length(U_Post.NEP_day)); %one indicates group Post
%create a second vector indicating the dataset from which the corresponding measurement is made
U_dataset = [U_dataset1 U_dataset2];
U_Group = cell(size(U_dataset)); % one row of data 
U_Group(U_dataset==0)={'Pre-Restoration'};
U_Group(U_dataset==1)={'Post-Restoration'};

U_length_Pre_NEP_day = length (U_Pre.NEP_day);
U_length_Post_NEP_day = length (U_Post.NEP_day);
U_Master_NEP_day = [U_Pre.NEP_day U_Post.NEP_day];% one row of data 
length(U_Master_NEP_day)==U_length_Pre_NEP_day+U_length_Post_NEP_day % 1= true


%p = kruskalwallis(x, group) - 
[U_NEP_day_p,tbl,stats] = kruskalwallis(U_Master_NEP_day, U_Group)

% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(U_Pre.NEP_day,U_Post.NEP_day) 


% Marker 32 NEC Day
clc 
close all
M_length_Pre_NEC_day = length (M_Pre.NEC_day);
M_length_Post_NEC_day = length (M_Post.NEC_day);
M_Master_NEC_day = [M_Pre.NEC_day M_Post.NEC_day];% one row of data 
length(M_Master_NEC_day)==M_length_Pre_NEC_day+M_length_Post_NEC_day % 1= true

M_dataset1 = zeros(1,length(M_Pre.NEC_day)); %Zero indicates group Pre
M_dataset2 = ones(1,length(M_Post.NEC_day)); %one indicates group Post
%create a second vector indicating the dataset from which the corresponding measurement is made
M_dataset = [M_dataset1 M_dataset2];
M_Group = cell(size(M_dataset)); % one row of data 
M_Group(M_dataset==0)={'Pre-Restoration'};
M_Group(M_dataset==1)={'Post-Restoration'};

%p = kruskalwallis(x, group) - test if data in each categorical group is from the
%same distribution
[M_NEC_day_p,tbl,stats] = kruskalwallis(M_Master_NEC_day, M_Group)

% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(M_Pre.NEC_day,M_Post.NEC_day)


% Marker 32 NEP Day
clc 
close all
M_length_Pre_NEP_day = length (M_Pre.NEP_day);
M_length_Post_NEP_day = length (M_Post.NEP_day);
M_Master_NEP_day = [M_Pre.NEP_day M_Post.NEP_day];% one row of data 
length(M_Master_NEP_day)==M_length_Pre_NEP_day+M_length_Post_NEP_day % 1= true

M_dataset1 = zeros(1,length(M_Pre.NEP_day)); %Zero indicates group Pre
M_dataset2 = ones(1,length(M_Post.NEP_day)); %one indicates group Post
%create a second vector indicating the dataset from which the corresponding measurement is made
M_dataset = [M_dataset1 M_dataset2];
M_Group = cell(size(M_dataset)); % one row of data 
M_Group(M_dataset==0)={'Pre-Restoration'};
M_Group(M_dataset==1)={'Post-Restoration'};

%p = kruskalwallis(x, group) - test if data in each categorical group is from the
%same distribution
[M_NEP_day_p,tbl,stats] = kruskalwallis(M_Master_NEP_day, M_Group)

% Mann-Whitney U-test == Wilcoxon rank sum test
% nonparametric test for equality of population medians of two independent samples X and Y
% null hypothesis that data in x and y are samples from continuous distributions with equal medians,
% h = 1 indicates a rejection of the null hypothesis, 
[p,h,stats] = ranksum(M_Pre.NEP_day,M_Post.NEP_day)


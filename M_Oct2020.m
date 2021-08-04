% Oct SeapHOx Data Analysis from Marker 32 Reef
% Michelle Platz - USF 
% 3/10/2021

% SeapHOx sensor deployed 9/2/2020 
    % SP datafile: 'M1001BMS.txt'
    % pump 1 height above benthos = 70 cm  
    % pump 2 height above benthos = 20 cm 
% ADP sensor deployed 9/02/2020 
    % ADCP datafile: 'M_9_02'
    % height from substrate to ADCP head = 18 cm
clear all
close all
clc

%% Initial look at data
% ***** create MOct_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
MOct_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('M1001BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
MOct_SPraw.DOXY = MOct_SPraw.O2SATPER.*calcO2sat(MOct_SPraw.MCAT_TC, MOct_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
MOct_SPraw.pHint_prelim = calc_dfet_pHint(MOct_SPraw.Vint, MOct_SPraw.DFET_TC, -0.4);

% ***** create MOct_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
MOct_SP = parse_to_pumpheights_ARM_2pump_Mar19(MOct_SPraw, [datenum('10-01-2020 08:30:00'), datenum('10-22-2020 10:59:30')]);

% Calculate Gradients 
MOct_SP = calc_TA_gradientV2(MOct_SP, 2368.31, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

MOct_SP.dDOXY = MOct_SP.DOXY(1,:) - MOct_SP.DOXY(2,:); %Oxygen Gradient
MOct_SP.dpH = MOct_SP.pH(1,:) - MOct_SP.pH(2,:); %pH Gradient 
MOct_SP.dTA = MOct_SP.TAtop - MOct_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
MOct_DateString = {'10/01/2020 12:00:00';'10/02/2020 12:00:00';'10/03/2020 12:00:00';...
    '10/04/2020 12:00:00';'10/05/2020 12:00:00';'10/06/2020 12:00:00';'10/07/2020 12:00:00';...
    '10/08/2020 12:00:00';'10/09/2020 12:00:00';'10/10/2020 12:00:00';'10/11/2020 12:00:00';...
    '10/12/2020 12:00:00';'10/13/2020 12:00:00';'10/14/2020 12:00:00';'10/15/2020 12:00:00';...
    '10/16/2020 12:00:00';'10/17/2020 12:00:00';'10/18/2020 12:00:00';'10/19/2020 12:00:00';...
    '10/20/2020 12:00:00';'10/21/2020 12:00:00';'10/22/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
MOct_tick = datenum(MOct_DateString,formatIn);

MOct_Xrange = [datenum('10-01-2020 10:30:00'), datenum('10-22-2020 10:59:30')];

figure
hold on; box on;
plot(MOct_SP.SDN, MOct_SP.dDOXY); %oxygen gradient 
plot(MOct_SP.SDN, MOct_SP.dTA); %TA gradient
plot(MOct_SP.SDN, zeros(size(MOct_SP.SDN))); %zero line
set(gca, 'xlim', MOct_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 Oct 2020 Unbinned DO and TA Gradients');

%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
MOct_SiteChar.SDN = MOct_SPraw.SDN;
MOct_SiteChar.TC = MOct_SPraw.OPT_TC;
MOct_SiteChar.PAR = MOct_SPraw.PAR;
MOct_SiteChar.PSAL = MOct_SPraw.PSAL;
MOct_SiteChar.Pres = MOct_SPraw.Pres;

%Plot pressure data to see when surface interval observations are
close all
figure 
hold on; 
plot(MOct_SiteChar.SDN, MOct_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
MOct_SiteChar.SDN = MOct_SPraw.SDN(139:end);
MOct_SiteChar.TC = MOct_SPraw.OPT_TC(139:end);
MOct_SiteChar.PAR = MOct_SPraw.PAR(139:end);
MOct_SiteChar.PSAL = MOct_SPraw.PSAL(139:end);
MOct_SiteChar.Pres = MOct_SPraw.Pres(139:end);

%extract full length of ADCP datafile  --> already done in Oct 
% % % MOct_ADfull=aquadoppraw2mat('U_9_02', 70, [datenum('09-02-2020 16:00:00'), datenum('10-30-2020 05:57:30')]);
% % % 
% % % % add AD variables to Site Char 
% % % MOct_SiteChar.AD_SDN = MOct_ADfull.SDN;
% % % MOct_SiteChar.AD_Pres = MOct_ADfull.Pres;
% % % MOct_SiteChar.AD_TC = MOct_ADfull.TC;
% % % MOct_SiteChar.bin_depth = MOct_ADfull.bin_depth;
% % % MOct_SiteChar.u = MOct_ADfull.u;
% % % MOct_SiteChar.v = MOct_ADfull.v;
% % % MOct_SiteChar.w = MOct_ADfull.w;
% % % MOct_SiteChar.uv = MOct_ADfull.uv;
% % % MOct_SiteChar.direction = MOct_ADfull.direction;
% % % 
% % % %Plot pressure data to see when surface interval observations are
% % % close all
% % % figure 
% % % hold on; 
% % % plot(MOct_SiteChar.AD_SDN, MOct_SiteChar.AD_Pres)
% % % 
% % % %clip ends of data to remove surfave interval observations 
% % % MOct_SiteChar.AD_SDN = MOct_ADfull.SDN(71:end);
% % % MOct_SiteChar.AD_Pres = MOct_ADfull.Pres(71:end);
% % % MOct_SiteChar.AD_TC = MOct_ADfull.TC(71:end);
% % % MOct_SiteChar.bin_depth = MOct_ADfull.bin_depth;
% % % MOct_SiteChar.u = MOct_ADfull.u(:,71:end);
% % % MOct_SiteChar.v = MOct_ADfull.v(:,71:end);
% % % MOct_SiteChar.w = MOct_ADfull.w(:,71:end);
% % % MOct_SiteChar.uv = MOct_ADfull.uv(:,71:end);
% % % MOct_SiteChar.direction = MOct_ADfull.direction(:,71:end);
% % % 
% % % % find U0 
% % % MOct_z1 = 0.70;
% % % MOct_z2 = 0.20;
% % % MOct_ADheight = 0.18;
% % % MOct_ADbin_depth_1m = 1-(MOct_ADheight);% = 0.82
% % % MOct_i1m = find(MOct_SiteChar.bin_depth==(0.82));
% % % MOct_SiteChar.U0 = MOct_SiteChar.uv(MOct_i1m,:);

% save data in separate datastructure
save('MOct20_SiteChar_2.mat', 'MOct_SiteChar' )


%% Constrain Xrange from graph results and extract good gradient data - 
close all 

MOct_good_Xrange = [datenum('10-01-2020 10:30:00'), datenum('10-05-2020 12:30:00')];

% plot to check range is correct
figure
hold on; box on;
plot(MOct_SP.SDN, MOct_SP.dDOXY); %oxygen gradient 
plot(MOct_SP.SDN, MOct_SP.dTA); %TA gradient
plot(MOct_SP.SDN, zeros(size(MOct_SP.SDN))); %zero line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 Oct 2020 Unbinned DO and TA Gradients');




%% Create good dataframe

close all


M_Oct_BMS_idx_start = find(MOct_SP.SDN==datenum('10-01-2020 10:30:00'))
M_Oct_BMS_idx_end = find(MOct_SP.SDN==datenum('10-05-2020 13:30:00'))

% Create new data vectors of just the good data
M_Oct_BMS_good_data = M_Oct_BMS_idx_start:M_Oct_BMS_idx_end;
Initial_data_points = length(M_Oct_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(MOct_SP);
for v = 1:length(vars)
    MOct_SP.(vars{v}) = (MOct_SP.(vars{v})(:,M_Oct_BMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: MOct_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
MOct_AD=aquadoppraw2mat('M_9_02', 70, [datenum('10-01-2020 10:30:00'), datenum('10-05-2020 13:45:00')]);%***********add 15 mins to SP end time

%averages data to the middle of the minute interval spacified 
MOct_ADavg = average_aquadopp(MOct_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.7m (pump 1) and 0.2m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.52m (bin 42) and 0.02m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
MOct_Sal_est = mean(MOct_SP.PSAL(1,3:end));

[MOct_ADavg] = ustar_from_aquadopp2(MOct_ADavg,[0.52 0.11], MOct_Sal_est); %bins adjusted 

clc
%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[MOct_ADavg] = ustar_McGillis_Method(MOct_ADavg, 0.70, 0.20, 42, 1);

%compare
close all
figure
hold on 
ustar_plot = plot(MOct_ADavg.SDN, MOct_ADavg.ustar, 'r');
ustar_WM_plot = plot(MOct_ADavg.SDN, MOct_ADavg.ustar_WM);
plot(MOct_ADavg.SDN, zeros(size(MOct_ADavg.SDN)),'k');
set(gca, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MOct_ADavg.SDN(1) MOct_ADavg.SDN(end)])
xlabel('Oct Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: MOct_BMS *****
% BMS = AD + SP datasets (full datasets)

ADavg_vars = fieldnames(MOct_ADavg);
for v = 1:length(ADavg_vars)
    MOct_BMS.(ADavg_vars{v}) = (MOct_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(MOct_SP);
for v = 1:length(SP_vars)
    MOct_BMS.(SP_vars{v}) = (MOct_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(MOct_BMS.SDN)
% min: 13-28-43-58 becuase SP was restarted in the field rather than on the
% minute
%% %% *************** Calculate Fluxes ****************

% actual pump heights  = 0.70m (pump 1) and 0.20m (pump 2) 
% 0.18m from substrate to ACDP head in Oct at U 
% adjusted height = 0.52 m (bin 42) and 0.02 m (too shallow) (bin 1)  
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[MOct_BMS] = calc_NCC_3(MOct_BMS,[0.52 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(MOct_BMS.DOXY(1,:));
[MOct_BMS] = calc_NCP_3(MOct_BMS, [0.52 0.11],C1guess); %estimate C1 guess using median DOXY(1,:) value  

% Plot NCP and NCC
%close all
figure
hold on; box on; 
NEPplot = plot(MOct_BMS.SDN, MOct_BMS.NEP);
NECplot = plot(MOct_BMS.SDN, MOct_BMS.NEC);
plot(MOct_BMS.SDN, zeros(size(MOct_BMS.SDN)),'k');
set(gca, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MOct_BMS.SDN(1) MOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Fluxes');

%McGillis method flux calculations 
[MOct_BMS] = calc_NCP_McGillis_Method(MOct_BMS, 0.70, 0.20, MOct_Sal_est);
[MOct_BMS] = calc_NCC_McGillis_Method(MOct_BMS, 0.70, 0.20, MOct_Sal_est);

% Plot NCP and NCC
% close all
% figure
% hold on; box on; 
% NEPplot = plot(MOct_BMS.SDN, MOct_BMS.NEP_WM);
% NECplot = plot(MOct_BMS.SDN, MOct_BMS.NEC_WM);
% plot(MOct_BMS.SDN, zeros(size(MOct_BMS.SDN)),'k');
% set(gca, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlim([MOct_BMS.SDN(1) MOct_BMS.SDN(end)])
% xlabel('Oct Days');
% ylabel('NCP or NCC [mmol/m2/hr]');
% legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
% title('Marker 32 Oct 2020 WM Fluxes');

%% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(MOct_BMS.SDN, MOct_BMS.NEP);
NEPplotWM = plot(MOct_BMS.SDN, MOct_BMS.NEP_WM);
plot(MOct_BMS.SDN, zeros(size(MOct_BMS.SDN)),'k');
set(gca, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MOct_BMS.SDN(1) MOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Marker 32 Oct 2020 NEP Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(MOct_BMS.SDN, MOct_BMS.NEC);
NECplotWM = plot(MOct_BMS.SDN, MOct_BMS.NEC_WM);
plot(MOct_BMS.SDN, zeros(size(MOct_BMS.SDN)),'k');
set(gca, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MOct_BMS.SDN(1) MOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Marker 32 Oct 2020 NEC Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
MOct_BMS.DOXYstd = std(MOct_BMS.DOXY);

% get DOXY std
MOct_idoxystd = find(MOct_BMS.DOXYstd > 2);
MOct_ihighdoxystd = [];
for i = 1:length(MOct_idoxystd)
    
    MOct_ihighdoxystd = vertcat(MOct_ihighdoxystd,[MOct_idoxystd(i)-1:1:MOct_idoxystd(i)+1]');
end
% get unique IDs
MOct_ihighdoxystd = unique(MOct_ihighdoxystd);
% remove 0's and out of index values
MOct_ihighdoxystd(MOct_ihighdoxystd==0) = [];
MOct_ihighdoxystd(MOct_ihighdoxystd> length(MOct_BMS.SDN)) = [];

% make it into index
trex = false(size(MOct_BMS.SDN));
trex(MOct_ihighdoxystd) = true;
MOct_ihighdoxystd = trex;
clear trex;

MOct_BMS.NEP_QC = MOct_BMS.NEP;
MOct_BMS.NEC_QC = MOct_BMS.NEC;
MOct_BMS.dDOXY_QC = MOct_BMS.dDOXY; %DO gradient
MOct_BMS.dTA_QC = MOct_BMS.dTA;     %TA gradient
MOct_BMS.NEP_WM_QC = MOct_BMS.NEP_WM;
MOct_BMS.NEC_WM_QC = MOct_BMS.NEC_WM;

% set observations when DOXYstd>0.8 to NaN
MOct_BMS.NEP_QC(MOct_ihighdoxystd) = NaN;  %NEP flux
MOct_BMS.NEC_QC(MOct_ihighdoxystd) = NaN;%NEC flux
MOct_BMS.dDOXY_QC(MOct_ihighdoxystd) = NaN; %DO gradient
MOct_BMS.dTA_QC(:,MOct_ihighdoxystd) = NaN;   %TA gradient
MOct_BMS.NEP_WM_QC(MOct_ihighdoxystd) = NaN;
MOct_BMS.NEC_WM_QC(:,MOct_ihighdoxystd) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(MOct_BMS.SDN, MOct_BMS.NEP, 'k');
NEPplotQC = plot(MOct_BMS.SDN, MOct_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(MOct_BMS.SDN, zeros(size(MOct_BMS.SDN)),'k');
set(gca, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MOct_BMS.SDN(1) MOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Fluxes');


%WM Plot
%close all
figure
hold on; box on;
NEPplot = plot(MOct_BMS.SDN, MOct_BMS.NEP_WM, 'k');
NEPplotQC = plot(MOct_BMS.SDN, MOct_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
plot(MOct_BMS.SDN, zeros(size(MOct_BMS.SDN)),'k');
set(gca, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MOct_BMS.SDN(1) MOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Fluxes WM');


%% Bin data to hourly intervals

X = floor(nanmin(MOct_BMS.SDN)):1/24:ceil(nanmax(MOct_BMS.SDN));

MOct_BMSbin.SDN = X;
% variables from 2 differnet heights
MOct_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.DOXY(1,:), X);
MOct_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.DOXY(2,:), X);

MOct_BMSbin.pH(1,:)       = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.pH(1,:), X);
MOct_BMSbin.pH(2,:)       = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.pH(2,:), X);

MOct_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.PSAL(1,:), X);
MOct_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.PSAL(2,:), X);

MOct_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.O2SATPER(1,:), X);
MOct_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.O2SATPER(2,:), X);

MOct_BMSbin.Pres(1,:)     = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.Pres(1,:), X);
MOct_BMSbin.Pres(2,:)     = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.Pres(2,:), X);

MOct_BMSbin.DENS(1,:)     = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.DENS(1,:), X);
MOct_BMSbin.DENS(2,:)     = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.DENS(2,:), X);

MOct_BMSbin.PAR(1,:)      = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.PAR(1,:), X);
MOct_BMSbin.PAR(2,:)      = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.PAR(2,:), X);

MOct_BMSbin.bin_depth     = MOct_BMS.bin_depth;

for i = 1:108
    MOct_BMSbin.uv(i,:)   = bin_data_to_X_GF(MOct_BMS.SDN,MOct_BMS.uv(i,:), X);
end
% bin data hourly. Vector variables 
MOct_BMSbin.PRES  = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.Pres, X);
MOct_BMSbin.U0       = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.U0, X);
MOct_BMSbin.DIR      = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.direction, X);
MOct_BMSbin.ustar    = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.ustar, X);
MOct_BMSbin.ustar_rm = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.ustar_runmean, X);
MOct_BMSbin.dTA      = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.dTA, X);
MOct_BMSbin.dTA_QC   = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.dTA_QC, X);
MOct_BMSbin.dpH      = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.dpH, X);
MOct_BMSbin.dDOXY    = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.dDOXY, X);
MOct_BMSbin.dDOXY_QC = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.dDOXY_QC, X);
MOct_BMSbin.NEP      = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.NEP, X);
MOct_BMSbin.NEP_QC   = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.NEP_QC, X);
MOct_BMSbin.NEC      = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.NEC, X);
MOct_BMSbin.NEC_QC   = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.NEC_QC, X);
MOct_BMSbin.NEP_WM   = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.NEP_WM, X);
MOct_BMSbin.NEP_WM_QC= bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.NEP_WM_QC, X);
MOct_BMSbin.NEC_WM   = bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.NEC_WM, X);
MOct_BMSbin.NEC_WM_QC= bin_data_to_X_GF(MOct_BMS.SDN, MOct_BMS.NEC_WM_QC, X);


%% Plot Binned Fluxes 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP, 'b'); 
% NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC, 'r-'); 
% plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Oct 2020 Hourly Binned Fluxes FluxFit');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC_QC, 'r-'); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned Fluxes FluxFit QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC_WM, 'r-'); 
% plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Oct 2020 WM Original Hourly Binned Fluxes');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC_WM_QC, 'r-'); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 WM Full QC Hourly Binned Fluxes');


%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
DOplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.dDOXY, 'c'); 
TAplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.dTA, 'k-'); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Binned Gradiets');
%% ***QC*** Remove sections when velocity is too slow
ibad = MOct_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow
sum(ibad)
% MOct_BMSbin.NEP(ibad) = NaN;
% MOct_BMSbin.NEC(ibad) = NaN;
MOct_BMSbin.NEP_QC(ibad) = NaN;
MOct_BMSbin.NEC_QC(ibad) = NaN;

MOct_BMSbin.dDOXY_QC(ibad) = NaN;
MOct_BMSbin.dTA_QC(ibad) = NaN;

% MOct_BMSbin.NEP_WM(ibad) = NaN;
% MOct_BMSbin.NEC_WM(ibad) = NaN;
MOct_BMSbin.NEP_WM_QC(ibad) = NaN;
MOct_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP, 'b'); 
% NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC, 'r-'); 
% plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Oct 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC_QC, 'r-'); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC_WM, 'r-'); 
% plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Oct 2020 Hourly Binned WM Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC_WM_QC, 'r-'); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned WM Fluxes Full QC');


%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
close all
figure
hold on; box on;
boxplot(MOct_BMSbin.NEP_QC)
ylabel('NEP')
title('Oct NEP Boxplots')

figure
hold on; box on;
boxplot(MOct_BMSbin.NEC_QC)
ylabel('NEC')
title('Oct NEC Boxplots')

figure
hold on; box on;
boxplot(MOct_BMSbin.dDOXY_QC)
ylabel('DO')
title('Oct dDO Boxplots')

figure
hold on; box on;
boxplot(MOct_BMSbin.dTA_QC)
ylabel('TA')
title('Oct dTA Boxplots')


% Outliers:
% 48:dTA outlier: 2.3
% 42: NEC outlier: -3.2
% 51: NEC outlier: 4.3
% 34 : NEC outlier: 6.2
    % bad profile
% 35: NEC outlier: 7.1
    % bad profile

datestr(MOct_BMSbin.SDN(42))


% Plot Profiles at outliers - 
close all 
for i = 35   %1:length(MAug2_BMSbin.SDN)
    figure (i)
    scatter(MOct_BMSbin.uv(1:108,i), MOct_BMSbin.bin_depth(1:108))
    title(['Cudjoe Oct Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end
%outliers to be removed: 34 35 
% MOct_BMSbin.NEP_QC(34) = NaN;
% MOct_BMSbin.NEC_QC(34) = NaN;
% MOct_BMSbin.dDOXY_QC(34) = NaN;
% MOct_BMSbin.dTA_QC(34) = NaN;
% 
% MOct_BMSbin.NEP_QC(35) = NaN;
% MOct_BMSbin.NEC_QC(35) = NaN;
% MOct_BMSbin.dDOXY_QC(35) = NaN;
% MOct_BMSbin.dTA_QC(35) = NaN;
close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC_QC, 'r-'); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned Fluxes Full QC');

%% Plot Diel Curves

% nepdbin = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEP, 24);
% necdbin = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin, 'bo', 'markersize', 3);
% plot(1:24, necdbin, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin,1), 'bo-');
% plot(1:24, nanmedian(necdbin,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('Oct Diel Plot 2020');
% xlabel('hour of day');

nepdbin_QC = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEC_QC, 24);
figure
hold on; box on;
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, nepdbin_QC, 'bo', 'markersize', 3);
plot(1:24, necdbin_QC, 'ro', 'markersize', 3);
plot(1:24, nanmedian(nepdbin_QC,1), 'bo-');
plot(1:24, nanmedian(necdbin_QC,1), 'ro-')
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
set(gca, 'XGrid', 'on');
xlim([0 24]);
ylabel(['NEP or \color{red}NEC']);
title('Oct Diel Plot 2020  QC');
xlabel('hour of day');

% WM data
% nepdbin_WM = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEC_WM, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('Oct 2020 WM');
% xlabel('hour of day');

nepdbin_WM_QC = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEC_WM_QC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM_QC, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM_QC, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM_QC,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM_QC,1), 'ro-')
% xlim([1 24]);
% ylabel(['NEP or \color{red}NEC']);
% title('Oct 2020 WM Full QC');
% xlabel('hour of day');



%% Extract Daytime data for Ratios
clc
% Extract daytime data using MOct_BMSbin.PAR
MOct_inight = MOct_BMSbin.PAR(1,:) < 5; %find all nightime datapoints 

%create new arrays for daytime data
MOct_BMSbin.SDN_day = MOct_BMSbin.SDN;
MOct_BMSbin.PAR_day = MOct_BMSbin.PAR(1,:);

MOct_BMSbin.NEP_day = MOct_BMSbin.NEP;
MOct_BMSbin.NEC_day = MOct_BMSbin.NEC;
MOct_BMSbin.NEP_day_QC = MOct_BMSbin.NEP_QC;
MOct_BMSbin.NEC_day_QC = MOct_BMSbin.NEC_QC;

MOct_BMSbin.dDOXY_day_QC = MOct_BMSbin.dDOXY_QC;      %DO Gradient
MOct_BMSbin.dTA_day_QC = MOct_BMSbin.dTA_QC;          %TA Gradient

MOct_BMSbin.NEP_WM_day = MOct_BMSbin.NEP_WM;
MOct_BMSbin.NEC_WM_day = MOct_BMSbin.NEC_WM;
MOct_BMSbin.NEP_WM_day_QC = MOct_BMSbin.NEP_WM_QC;
MOct_BMSbin.NEC_WM_day_QC = MOct_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
MOct_BMSbin.SDN_day(MOct_inight) = NaN;
MOct_BMSbin.PAR_day (MOct_inight) = NaN;

MOct_BMSbin.NEP_day(MOct_inight) = NaN;
MOct_BMSbin.NEC_day(MOct_inight) = NaN;
MOct_BMSbin.NEP_day_QC(MOct_inight) = NaN;
MOct_BMSbin.NEC_day_QC(MOct_inight) = NaN;

MOct_BMSbin.dDOXY_day_QC(MOct_inight) = NaN;      %DO Gradient
MOct_BMSbin.dTA_day_QC(MOct_inight) = NaN;        %TA Gradient

MOct_BMSbin.NEP_WM_day(MOct_inight) = NaN;
MOct_BMSbin.NEC_WM_day(MOct_inight) = NaN;
MOct_BMSbin.NEP_WM_day_QC(MOct_inight) = NaN;
MOct_BMSbin.NEC_WM_day_QC(MOct_inight) = NaN;

%Plot to check only nighttime points removed
figure 
hold on
scatter(MOct_BMSbin.SDN, MOct_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MOct_BMSbin.SDN_day, MOct_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
MOct_BMSbin.NEP_day(isnan(MOct_BMSbin.NEP_day))=[];
MOct_BMSbin.NEC_day(isnan(MOct_BMSbin.NEC_day))=[];
MOct_BMSbin.NEP_day_QC(isnan(MOct_BMSbin.NEP_day_QC))=[];
MOct_BMSbin.NEC_day_QC(isnan(MOct_BMSbin.NEC_day_QC))=[];

MOct_BMSbin.dDOXY_day_QC(isnan(MOct_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
MOct_BMSbin.dTA_day_QC(isnan(MOct_BMSbin.dTA_day_QC))=[];       %TA Gradient

MOct_BMSbin.NEP_WM_day(isnan(MOct_BMSbin.NEP_WM_day))=[];
MOct_BMSbin.NEC_WM_day(isnan(MOct_BMSbin.NEC_WM_day))=[];
MOct_BMSbin.NEP_WM_day_QC(isnan(MOct_BMSbin.NEP_WM_day_QC))=[];
MOct_BMSbin.NEC_WM_day_QC(isnan(MOct_BMSbin.NEC_WM_day_QC))=[];

% create nighttime hours datasets
MOct_BMSbin.SDN_night = MOct_BMSbin.SDN;
MOct_BMSbin.PAR_night = MOct_BMSbin.PAR(1,:);

MOct_BMSbin.NEP_night = MOct_BMSbin.NEP;
MOct_BMSbin.NEC_night = MOct_BMSbin.NEC;
MOct_BMSbin.NEP_night_QC = MOct_BMSbin.NEP_QC;
MOct_BMSbin.NEC_night_QC = MOct_BMSbin.NEC_QC;

MOct_BMSbin.dDOXY_night_QC = MOct_BMSbin.dDOXY_QC;      %DO Gradient
MOct_BMSbin.dTA_night_QC = MOct_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
MOct_BMSbin.SDN_night=MOct_BMSbin.SDN_night(MOct_inight);
MOct_BMSbin.PAR_night=MOct_BMSbin.PAR_night(MOct_inight);

MOct_BMSbin.NEP_night_QC=MOct_BMSbin.NEP_night_QC(MOct_inight);
MOct_BMSbin.NEC_night_QC=MOct_BMSbin.NEC_night_QC(MOct_inight);

MOct_BMSbin.dDOXY_night_QC=MOct_BMSbin.dDOXY_night_QC(MOct_inight);      %DO Gradient
MOct_BMSbin.dTA_night_QC=MOct_BMSbin.dTA_night_QC(MOct_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(MOct_BMSbin.SDN, MOct_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MOct_BMSbin.SDN_night, MOct_BMSbin.PAR_night, 'r.'); % night plot

%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(MOct_BMSbin.NEP_day,MOct_BMSbin.NEC_day);
% MOct_BMSbin.Reg_Line = m*MOct_BMSbin.NEP_day + b;
% MOct_BMSbin.Ratio = m;
% MOct_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(MOct_BMSbin.NEP_day,MOct_BMSbin.NEC_day,'o')
% plot(MOct_BMSbin.NEP_day,MOct_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Oct 2020 Pre-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MOct_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MOct_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(MOct_BMSbin.NEP_day_QC,MOct_BMSbin.NEC_day_QC);
MOct_BMSbin.Reg_Line_QC = m_QC*MOct_BMSbin.NEP_day_QC + b_QC;
MOct_BMSbin.Ratio_QC = m_QC;
MOct_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(MOct_BMSbin.NEP_day_QC,MOct_BMSbin.NEC_day_QC,'o')
plot(MOct_BMSbin.NEP_day_QC,MOct_BMSbin.Reg_Line_QC,'r')
xlabel('NCP');
ylabel('NCC');
title('Marker 32 Oct 2020 Post-Restoration NCC:NCP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MOct_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MOct_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(MOct_BMSbin.NEP_WM_day,MOct_BMSbin.NEC_WM_day);
% MOct_BMSbin.Reg_Line_WM = m_WM*MOct_BMSbin.NEP_WM_day + b_WM;
% MOct_BMSbin.Ratio_WM = m_WM;
% MOct_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(MOct_BMSbin.NEP_WM_day,MOct_BMSbin.NEC_WM_day,'o')
% plot(MOct_BMSbin.NEP_WM_day,MOct_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Oct 2020 Pre-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MOct_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MOct_BMSbin.R2_WM)


[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(MOct_BMSbin.NEP_WM_day_QC,MOct_BMSbin.NEC_WM_day_QC);
MOct_BMSbin.Reg_Line_WM_QC = m_WM_QC*MOct_BMSbin.NEP_WM_day_QC + b_WM_QC;
MOct_BMSbin.Ratio_WM_QC = m_WM_QC;
MOct_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(MOct_BMSbin.NEP_WM_day_QC,MOct_BMSbin.NEC_WM_day_QC,'o')
% plot(MOct_BMSbin.NEP_WM_day_QC,MOct_BMSbin.Reg_Line_WM_QC,'r')
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Oct 2020 Post-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MOct_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MOct_BMSbin.R2_WM_QC)
% 


%% For NEC:NEP Regressions Using Gradients
close all 
clc

% multiply o2 gradient by -1 for O2 production
MOct_BMSbin.dDOXY_Reg = -1.*MOct_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
MOct_BMSbin.dTA_Reg = 0.5.*MOct_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(MOct_BMSbin.dDOXY_day_QC);
figure 
hold on 
DOday = plot(1:Xlength, MOct_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, MOct_BMSbin.dDOXY_Reg);
xlabel('Oct Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, MOct_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, MOct_BMSbin.dTA_Reg);
xlabel('Oct Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(MOct_BMSbin.dDOXY_Reg, MOct_BMSbin.dTA_Reg);
MOct_BMSbin.Reg_Line_G = m_G*MOct_BMSbin.dDOXY_Reg + b_G;
MOct_BMSbin.Ratio_G = m_G;
MOct_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(MOct_BMSbin.dDOXY_Reg,MOct_BMSbin.dTA_Reg,'o')
plot(MOct_BMSbin.dDOXY_Reg,MOct_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Marker 32 Oct 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MOct_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MOct_BMSbin.R2_G)


clc
disp('Finished with Oct Metabolism Calculations');

save('MOct20_2.mat', 'MOct_BMS', 'MOct_BMSbin');

%% Plot Profiles 
% plot for profile within pump heights
% close all 
% for i =1:100 %length(MOct_ADavg.SDN)
%     figure (i)
%     scatter(MOct_ADavg.uv(1:108,i), MOct_ADavg.bin_depth(1:108))
%     title(['Marker 32 Oct Velocity Profile Number ',num2str(i),])
%     xlabel('Velocity (m/s)');
%     ylabel('Height (m)');
% end


%% Subplots 
close all
clc

sgtitle('Marker 32 October 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'10/01 12:00';'10/02 12:00';'10/03 12:00';...
    '10/04 12:00';'10/05 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
% legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(MOct_BMSbin.SDN, MOct_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(MOct_BMSbin.SDN, zeros(size(MOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
xticklabels({'10/01 12:00';'10/02 12:00';'10/03 12:00';...
    '10/04 12:00';'10/05 12:00'})
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MOct_BMSbin.SDN, MOct_BMSbin.NEC_QC, 24);
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, nepdbin_QC, 'bo', 'markersize', 3);
plot(1:24, necdbin_QC, 'ro', 'markersize', 3);
plot(1:24, nanmedian(nepdbin_QC,1), 'bo-');
plot(1:24, nanmedian(necdbin_QC,1), 'ro-')
ylabel(['\color{blue}NEP \color{black} or \color{red}NEC']);
title('Diel Plot');
xlabel('hour of day');

subplot(3,3,8); % Ratio Plot using fluxes 
hold on; box on;
plot(MOct_BMSbin.NEP_day_QC,MOct_BMSbin.NEC_day_QC,'o')
plot(MOct_BMSbin.NEP_day_QC,MOct_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(MOct_BMSbin.Ratio_QC,2);
str2 = num2str(MOct_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.075, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.075, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(MOct_BMSbin.dDOXY_Reg,MOct_BMSbin.dTA_Reg,'o')
plot(MOct_BMSbin.dDOXY_Reg,MOct_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
title('dTA:dDO Ratio from Gradients');
str1 = num2str(MOct_BMSbin.Ratio_G,2);
str2 = num2str(MOct_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.075, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.075, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')

%% Site Characterizations 
close all
figure
hold on; box on;
plot(MOct_BMS.SDN, MOct_BMS.U0); %Velocity 
% plot(MOct_BMS.SDN, MOct_SP.dTA); %TA gradient
plot(MOct_BMS.SDN, zeros(size(MOct_SP.SDN))); %zero line
set(gca, 'xlim', MOct_good_Xrange, 'XTick', MOct_tick, 'xticklabel', MOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('Velocity');
% legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 Oct 2020 Site Characterizations');



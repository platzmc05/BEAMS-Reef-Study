% Aug_1 SeapHOx Data Analysis from Marker 32 Ledge Reef
% Michelle Platz - USF 
% 3/10/2021

% SeapHOx sensor deployed 8/05/2020 at 11 am EST
    % SP datafile: 'M_805BMS.txt'
    % pump 1 height above benthos = 70 cm  
    % pump 2 height above benthos = 20 cm 
% ADP sensor deployed 8/05/2020 at 9:15 am EST
    % ADCP datafile: 'M_8_05'
    % height from substrate to ADCP head = 18 cm

close all
clc
clear all

%% Initial look at data
% ***** create MAug1_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
MAug1_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('M_805BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
MAug1_SPraw.DOXY = MAug1_SPraw.O2SATPER.*calcO2sat(MAug1_SPraw.MCAT_TC, MAug1_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
MAug1_SPraw.pHint_prelim = calc_dfet_pHint(MAug1_SPraw.Vint, MAug1_SPraw.DFET_TC, -0.4);

% ***** create MAug1_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
MAug1_SP = parse_to_pumpheights_ARM_2pump_Mar19(MAug1_SPraw, [datenum('08-05-2020 12:00:00'), datenum('08-12-2020 09:29:30')]);

% Calculate Gradients 
MAug1_SP = calc_TA_gradientV2(MAug1_SP, 2368.31, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

MAug1_SP.dDOXY = MAug1_SP.DOXY(1,:) - MAug1_SP.DOXY(2,:); %Oxygen Gradient
MAug1_SP.dpH = MAug1_SP.pH(1,:) - MAug1_SP.pH(2,:); %pH Gradient 
MAug1_SP.dTA = MAug1_SP.TAtop - MAug1_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
MAug1_DateString = {'08/05/2020 12:00:00';'08/06/2020 12:00:00';'08/07/2020 12:00:00';'08/08/2020 12:00:00';'08/09/2020 12:00:00';...
    '08/10/2020 12:00:00';'08/11/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
MAug1_tick = datenum(MAug1_DateString,formatIn);

MAug1_Xrange = [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 09:29:30')];

figure
hold on; box on;
plot(MAug1_SP.SDN, MAug1_SP.dDOXY); %oxygen gradient 
plot(MAug1_SP.SDN, MAug1_SP.dTA); %TA gradient
plot(MAug1_SP.SDN, zeros(size(MAug1_SP.SDN))); %zero line
set(gca, 'xlim', MAug1_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 Aug 2020 Unbinned DO and TA Gradients');

%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
MAug1_SiteChar.SDN = MAug1_SPraw.SDN;
MAug1_SiteChar.TC = MAug1_SPraw.OPT_TC;
MAug1_SiteChar.PAR = MAug1_SPraw.PAR;
MAug1_SiteChar.PSAL = MAug1_SPraw.PSAL;
MAug1_SiteChar.Pres = MAug1_SPraw.Pres;

close all
figure 
hold on; 
plot(MAug1_SiteChar.SDN, MAug1_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
MAug1_SiteChar.SDN = MAug1_SPraw.SDN(205:20207);
MAug1_SiteChar.TC = MAug1_SPraw.OPT_TC(205:20207);
MAug1_SiteChar.PAR = MAug1_SPraw.PAR(205:20207);
MAug1_SiteChar.PSAL = MAug1_SPraw.PSAL(205:20207);
MAug1_SiteChar.Pres = MAug1_SPraw.Pres(205:20207);

%extract full length of ADCP datafile  
MAug1_ADfull=aquadoppraw2mat('M_8_05', 70, [datenum('08-05-2020 09:00:00'), datenum('09-01-2020 09:00:00')]);

% add AD variables to Site Char 
MAug1_SiteChar.AD_SDN = MAug1_ADfull.SDN;
MAug1_SiteChar.AD_Pres = MAug1_ADfull.Pres;
MAug1_SiteChar.AD_TC = MAug1_ADfull.TC;
MAug1_SiteChar.bin_depth = MAug1_ADfull.bin_depth;
MAug1_SiteChar.u = MAug1_ADfull.u;
MAug1_SiteChar.v = MAug1_ADfull.v;
MAug1_SiteChar.w = MAug1_ADfull.w;
MAug1_SiteChar.uv = MAug1_ADfull.uv;
MAug1_SiteChar.direction = MAug1_ADfull.direction;

%Plot to see when surface interval observations are
close all
figure 
hold on; 
plot(MAug1_SiteChar.AD_SDN, MAug1_SiteChar.AD_Pres)

%clip ends of data to remove surfave interval observations 
MAug1_SiteChar.AD_SDN = MAug1_ADfull.SDN(210:end);
MAug1_SiteChar.AD_Pres = MAug1_ADfull.Pres(210:end);
MAug1_SiteChar.AD_TC = MAug1_ADfull.TC(210:end);
MAug1_SiteChar.bin_depth = MAug1_ADfull.bin_depth;
MAug1_SiteChar.u = MAug1_ADfull.u(:,210:end);
MAug1_SiteChar.v = MAug1_ADfull.v(:,210:end);
MAug1_SiteChar.w = MAug1_ADfull.w(:,210:end);
MAug1_SiteChar.uv = MAug1_ADfull.uv(:,210:end);
MAug1_SiteChar.direction = MAug1_ADfull.direction(:,210:end);

% find U0 
MAug1_z1 = 0.70;
MAug1_z2 = 0.20;
MAug1_ADheight = 0.18;
MAug1_ADbin_depth_1m = 1-0.18;% = 0.82
MAug1_i1m = find(MAug1_SiteChar.bin_depth==(0.82));
MAug1_SiteChar.U0 = MAug1_SiteChar.uv(MAug1_i1m,:);

% save data in separate datastructure
save('MAug120_SiteChar_2.mat', 'MAug1_SiteChar' )


%% Constrain Xrange from graph results and extract good gradient data - 
close all 

MAug1_good_Xrange = [datenum('08-05-2020 12:00:00'), datenum('08-12-2020 09:00:00')];

% plot to check range is correct
figure
hold on; box on;
plot(MAug1_SP.SDN, MAug1_SP.dDOXY); %oxygen gradient 
plot(MAug1_SP.SDN, MAug1_SP.dTA); %TA gradient
plot(MAug1_SP.SDN, zeros(size(MAug1_SP.SDN))); %zero line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 Aug 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all

M_Aug_1_BMS_idx_start = find(MAug1_SP.SDN==datenum('08-05-2020 12:00:00'))
M_Aug_1_BMS_idx_end = find(MAug1_SP.SDN==datenum('08-12-2020 09:00:00'))

% Create new data vectors of just the good data
M_Aug_1_BMS_good_data = M_Aug_1_BMS_idx_start:M_Aug_1_BMS_idx_end;
Initial_data_points = length(M_Aug_1_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(MAug1_SP);
for v = 1:length(vars)
    MAug1_SP.(vars{v}) = (MAug1_SP.(vars{v})(:,M_Aug_1_BMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: MAug1_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
MAug1_AD=aquadoppraw2mat('M_8_05', 70, [datenum('08-05-2020 12:00:00'), datenum('08-12-2020 09:15:00')]);

%averages data to the middle of the minute interval spacified 
MAug1_ADavg = average_aquadopp(MAug1_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.7m (pump 1) and 0.2m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.52m (bin 42) and 0.02m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
MAug1_Sal_est = mean(MAug1_SP.PSAL(1,3:end));

[MAug1_ADavg] = ustar_from_aquadopp2(MAug1_ADavg,[0.52 0.11], MAug1_Sal_est); %bins adjusted 

clc
%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[MAug1_ADavg] = ustar_McGillis_Method(MAug1_ADavg, 0.70, 0.20, 42, 1);

%compare
close all
figure
hold on 
ustar_plot = plot(MAug1_ADavg.SDN, MAug1_ADavg.ustar, 'r');
ustar_WM_plot = plot(MAug1_ADavg.SDN, MAug1_ADavg.ustar_WM);
plot(MAug1_ADavg.SDN, zeros(size(MAug1_ADavg.SDN)),'k');
set(gca, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug1_ADavg.SDN(1) MAug1_ADavg.SDN(end)])
xlabel('Aug Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: MAug1_BMS *****

ADavg_vars = fieldnames(MAug1_ADavg);
for v = 1:length(ADavg_vars)
    MAug1_BMS.(ADavg_vars{v}) = (MAug1_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(MAug1_SP);
for v = 1:length(SP_vars)
    MAug1_BMS.(SP_vars{v}) = (MAug1_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(MAug1_BMS.SDN)


%% %% *************** Calculate Fluxes ****************

% actual pump heights  = 0.70m (pump 1) and 0.20m (pump 2) 
% 0.18m from substrate to ACDP head in Aug at U 
% adjusted height = 0.52 m (bin 42) and 0.02 m (too shallow) (bin 1)  
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[MAug1_BMS] = calc_NCC_3(MAug1_BMS,[0.52 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(MAug1_BMS.DOXY(1,2:end));
[MAug1_BMS] = calc_NCP_3(MAug1_BMS, [0.52 0.11],C1guess); %C1 guess - DOXY(1,x) guess is current saved, 

% Plot NCP and NCC
%close all
figure
hold on; box on; 
NEPplot = plot(MAug1_BMS.SDN, MAug1_BMS.NEP);
NECplot = plot(MAug1_BMS.SDN, MAug1_BMS.NEC);
plot(MAug1_BMS.SDN, zeros(size(MAug1_BMS.SDN)),'k');
set(gca, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug1_BMS.SDN(1) MAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');

%McGillis method flux calculations 
[MAug1_BMS] = calc_NCP_McGillis_Method(MAug1_BMS, 0.70, 0.20, MAug1_Sal_est);
[MAug1_BMS] = calc_NCC_McGillis_Method(MAug1_BMS, 0.70, 0.20, MAug1_Sal_est);

% Plot NCP and NCC
% close all
figure
hold on; box on; 
NEPplot = plot(MAug1_BMS.SDN, MAug1_BMS.NEP_WM);
NECplot = plot(MAug1_BMS.SDN, MAug1_BMS.NEC_WM);
plot(MAug1_BMS.SDN, zeros(size(MAug1_BMS.SDN)),'k');
set(gca, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug1_BMS.SDN(1) MAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 WM Fluxes');

%% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(MAug1_BMS.SDN, MAug1_BMS.NEP);
NEPplotWM = plot(MAug1_BMS.SDN, MAug1_BMS.NEP_WM);
plot(MAug1_BMS.SDN, zeros(size(MAug1_BMS.SDN)),'k');
set(gca, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug1_BMS.SDN(1) MAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(MAug1_BMS.SDN, MAug1_BMS.NEC);
NECplotWM = plot(MAug1_BMS.SDN, MAug1_BMS.NEC_WM);
plot(MAug1_BMS.SDN, zeros(size(MAug1_BMS.SDN)),'k');
set(gca, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug1_BMS.SDN(1) MAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
MAug1_BMS.DOXYstd = std(MAug1_BMS.DOXY);

% get DOXY std
MAug1_idoxystd = find(MAug1_BMS.DOXYstd > 2);% 58.4% of data is greater than 0.8 stdev
MAug1_ihighdoxystd = [];
for i = 1:length(MAug1_idoxystd)
    
    MAug1_ihighdoxystd = vertcat(MAug1_ihighdoxystd,[MAug1_idoxystd(i)-1:1:MAug1_idoxystd(i)+1]');
end
% get unique IDs
MAug1_ihighdoxystd = unique(MAug1_ihighdoxystd);
% remove 0's and out of index values
MAug1_ihighdoxystd(MAug1_ihighdoxystd==0) = [];
MAug1_ihighdoxystd(MAug1_ihighdoxystd> length(MAug1_BMS.SDN)) = [];

% make it into index
trex = false(size(MAug1_BMS.SDN));
trex(MAug1_ihighdoxystd) = true;
MAug1_ihighdoxystd = trex;
clear trex;

MAug1_BMS.NEP_QC = MAug1_BMS.NEP;
MAug1_BMS.NEC_QC = MAug1_BMS.NEC;
MAug1_BMS.dDOXY_QC = MAug1_BMS.dDOXY; %DO gradient
MAug1_BMS.dTA_QC = MAug1_BMS.dTA;     %TA gradient
MAug1_BMS.NEP_WM_QC = MAug1_BMS.NEP_WM;
MAug1_BMS.NEC_WM_QC = MAug1_BMS.NEC_WM;

% set observations when DOXYstd>0.8 to NaN
MAug1_BMS.NEP_QC(MAug1_ihighdoxystd) = NaN;
MAug1_BMS.NEC_QC(:,MAug1_ihighdoxystd) = NaN;
MAug1_BMS.dDOXY_QC(MAug1_ihighdoxystd) = NaN; %DO gradient
MAug1_BMS.dTA_QC(:,MAug1_ihighdoxystd) = NaN;   %TA gradient
MAug1_BMS.NEP_WM_QC(MAug1_ihighdoxystd) = NaN;
MAug1_BMS.NEC_WM_QC(:,MAug1_ihighdoxystd) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(MAug1_BMS.SDN, MAug1_BMS.NEP, 'k');
NEPplotQC = plot(MAug1_BMS.SDN, MAug1_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(MAug1_BMS.SDN, zeros(size(MAug1_BMS.SDN)),'k');
set(gca, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug1_BMS.SDN(1) MAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');


%WM Plot
%close all
figure
hold on; box on;
NEPplot = plot(MAug1_BMS.SDN, MAug1_BMS.NEP_WM, 'k');
NEPplotQC = plot(MAug1_BMS.SDN, MAug1_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
plot(MAug1_BMS.SDN, zeros(size(MAug1_BMS.SDN)),'k');
set(gca, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug1_BMS.SDN(1) MAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');


%% Bin data to hourly intervals

X = floor(nanmin(MAug1_BMS.SDN)):1/24:ceil(nanmax(MAug1_BMS.SDN));

MAug1_BMSbin.SDN = X;
% variables from 2 differnet heights
MAug1_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.DOXY(1,:), X);
MAug1_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.DOXY(2,:), X);

MAug1_BMSbin.pH(1,:)       = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.pH(1,:), X);
MAug1_BMSbin.pH(2,:)       = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.pH(2,:), X);

MAug1_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.PSAL(1,:), X);
MAug1_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.PSAL(2,:), X);

MAug1_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.O2SATPER(1,:), X);
MAug1_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.O2SATPER(2,:), X);

MAug1_BMSbin.Pres(1,:)     = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.Pres(1,:), X);
MAug1_BMSbin.Pres(2,:)     = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.Pres(2,:), X);

MAug1_BMSbin.DENS(1,:)     = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.DENS(1,:), X);
MAug1_BMSbin.DENS(2,:)     = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.DENS(2,:), X);

MAug1_BMSbin.PAR(1,:)      = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.PAR(1,:), X);
MAug1_BMSbin.PAR(2,:)      = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.PAR(2,:), X);

MAug1_BMSbin.bin_depth     = MAug1_BMS.bin_depth;

for i = 1:108
    MAug1_BMSbin.uv(i,:)   = bin_data_to_X_GF(MAug1_BMS.SDN,MAug1_BMS.uv(i,:), X);
end
% bin data hourly. Vector variables 
MAug1_BMSbin.PRES  = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.Pres, X);
MAug1_BMSbin.U0       = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.U0, X);
MAug1_BMSbin.DIR      = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.direction, X);
MAug1_BMSbin.ustar    = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.ustar, X);
MAug1_BMSbin.ustar_rm = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.ustar_runmean, X);
MAug1_BMSbin.dTA      = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.dTA, X);
MAug1_BMSbin.dTA_QC   = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.dTA_QC, X);
MAug1_BMSbin.dpH      = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.dpH, X);
MAug1_BMSbin.dDOXY    = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.dDOXY, X);
MAug1_BMSbin.dDOXY_QC = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.dDOXY_QC, X);
MAug1_BMSbin.NEP      = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.NEP, X);
MAug1_BMSbin.NEP_QC   = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.NEP_QC, X);
MAug1_BMSbin.NEC      = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.NEC, X);
MAug1_BMSbin.NEC_QC   = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.NEC_QC, X);
MAug1_BMSbin.NEP_WM   = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.NEP_WM, X);
MAug1_BMSbin.NEP_WM_QC= bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.NEP_WM_QC, X);
MAug1_BMSbin.NEC_WM   = bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.NEC_WM, X);
MAug1_BMSbin.NEC_WM_QC= bin_data_to_X_GF(MAug1_BMS.SDN, MAug1_BMS.NEC_WM_QC, X);


%% Plot Binned Fluxes 
close all
clc
% figure
% hold on; box on;
% NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP, 'b'); 
% NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC, 'r-'); 
% plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 Hourly Binned Fluxes');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_QC, 'r-'); 
plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Hourly Binned Fluxes');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_WM, 'r-'); 
% plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 WM Original Hourly Binned Fluxes');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_WM_QC, 'r-'); 
plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 WM Full QC Hourly Binned Fluxes');


%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
DOplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.dDOXY, 'c'); 
TAplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.dTA, 'k-'); 
plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug1 Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 Aug1 2020 Binned Gradiets');



%% ***QC*** Remove sections when velocity is too slow
ibad = MAug1_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow

% MAug1_BMSbin.NEP(ibad) = NaN;
% MAug1_BMSbin.NEC(ibad) = NaN;
MAug1_BMSbin.NEP_QC(ibad) = NaN;
MAug1_BMSbin.NEC_QC(ibad) = NaN;

MAug1_BMSbin.dDOXY_QC(ibad) = NaN;
MAug1_BMSbin.dTA_QC(ibad) = NaN;

% MAug1_BMSbin.NEP_WM(ibad) = NaN;
% MAug1_BMSbin.NEC_WM(ibad) = NaN;
MAug1_BMSbin.NEP_WM_QC(ibad) = NaN;
MAug1_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
clc
% figure
% hold on; box on;
% NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP, 'b'); 
% NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC, 'r-'); 
% plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_QC, 'r-'); 
plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_WM, 'r-'); 
% plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 Hourly Binned WM Fluxes');

% clc
% figure
% hold on; box on;
% NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_WM_QC, 'r-'); 
% plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 Hourly Binned WM Fluxes Full QC');

%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
figure
hold on; box on;
boxplot(MAug1_BMSbin.NEP_QC)
ylabel('NEP')
title('Aug1 NEP Boxplots')

figure
hold on; box on;
boxplot(MAug1_BMSbin.NEC_QC)
ylabel('NEC')
title('Aug1 NEC Boxplots')

figure
hold on; box on;
boxplot(MAug1_BMSbin.dDOXY_QC)
ylabel('DO')
title('Aug1 dDO Boxplots')

figure
hold on; box on;
boxplot(MAug1_BMSbin.dTA_QC)
ylabel('TA')
title('Aug1 dTA Boxplots')

% Outliers 

% 64: NEP outlier (value: -30)
    %'07-Aug-2020 15:00:00' - 3pm, value unlikely
    % OUTLIER 
% 67: NEP outlier (value: -33)
    %'07-Aug-2020 18:00:00' - 6pm, value unlikely
    % OUTLIER 
% 15: NEC outlier (value: 12)

% 38: NEC outlier (value: 10)

% 18: NEC outlier (value:10) 

% 87: dTA outlier (value: 4.1)
    % BAD PROFILE 
    % OUTLIER

% 178: dTA outlier (value: 4)
    % BAD PROFILE 
    % OUTLIER

datestr(MAug1_BMSbin.SDN(18))


% Plot Profiles at outliers - 
close all 
for i = 178   %1:length(MAug2_BMSbin.SDN)
    figure (i)
    scatter(MAug1_BMSbin.uv(1:108,i), MAug1_BMSbin.bin_depth(1:108))
    title(['Cudjoe Aug1 Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end


%outliers to be removed: 64 67 87 178 
% MAug1_BMSbin.NEP_QC(64) = NaN;
% MAug1_BMSbin.NEC_QC(64) = NaN;
% MAug1_BMSbin.dDOXY_QC(64) = NaN;
% MAug1_BMSbin.dTA_QC(64) = NaN;
% 
% MAug1_BMSbin.NEP_QC(67) = NaN;
% MAug1_BMSbin.NEC_QC(67) = NaN;
% MAug1_BMSbin.dDOXY_QC(67) = NaN;
% MAug1_BMSbin.dTA_QC(67) = NaN;
% 
% MAug1_BMSbin.NEP_QC(87) = NaN;
% MAug1_BMSbin.NEC_QC(87) = NaN;
% MAug1_BMSbin.dDOXY_QC(87) = NaN;
% MAug1_BMSbin.dTA_QC(87) = NaN;
% 
% MAug1_BMSbin.NEP_QC(178) = NaN;
% MAug1_BMSbin.NEC_QC(178) = NaN;
% MAug1_BMSbin.dDOXY_QC(178) = NaN;
% MAug1_BMSbin.dTA_QC(178) = NaN;
% % 
close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('Oct Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_QC, 'r-'); 
plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned Fluxes Full QC');

%% Plot Diel Curves
close all
% nepdbin = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP, 24);
% necdbin = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin, 'bo', 'markersize', 3);
% plot(1:24, necdbin, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin,1), 'bo-');
% plot(1:24, nanmedian(necdbin,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('Aug Diel Plot 2020');
% xlabel('hour of day');

nepdbin_QC = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_QC, 24);
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
title('Aug Diel Plot 2020  QC');
xlabel('hour of day');

% WM data
% nepdbin_WM = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_WM, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('Aug 2020 WM');
% xlabel('hour of day');

nepdbin_WM_QC = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_WM_QC, 24);
figure
hold on; box on;
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, nepdbin_WM_QC, 'bo', 'markersize', 3);
plot(1:24, necdbin_WM_QC, 'ro', 'markersize', 3);
plot(1:24, nanmedian(nepdbin_WM_QC,1), 'bo-');
plot(1:24, nanmedian(necdbin_WM_QC,1), 'ro-')
xlim([1 24]);
ylabel(['NEP or \color{red}NEC']);
title('Aug 2020 WM Full QC');
xlabel('hour of day');



%% Extract Daytime data for Ratios
clc
% Extract daytime data using MAug1_BMSbin.PAR
MAug1_inight = MAug1_BMSbin.PAR(1,:) < 1; %find all nightime datapoints 

%create new arrays for daytime data
MAug1_BMSbin.SDN_day = MAug1_BMSbin.SDN;
MAug1_BMSbin.PAR_day = MAug1_BMSbin.PAR(1,:);

MAug1_BMSbin.NEP_day = MAug1_BMSbin.NEP;
MAug1_BMSbin.NEC_day = MAug1_BMSbin.NEC;
MAug1_BMSbin.NEP_day_QC = MAug1_BMSbin.NEP_QC;
MAug1_BMSbin.NEC_day_QC = MAug1_BMSbin.NEC_QC;

MAug1_BMSbin.dDOXY_day_QC = MAug1_BMSbin.dDOXY_QC;      %DO Gradient
MAug1_BMSbin.dTA_day_QC = MAug1_BMSbin.dTA_QC;          %TA Gradient

MAug1_BMSbin.NEP_WM_day = MAug1_BMSbin.NEP_WM;
MAug1_BMSbin.NEC_WM_day = MAug1_BMSbin.NEC_WM;
MAug1_BMSbin.NEP_WM_day_QC = MAug1_BMSbin.NEP_WM_QC;
MAug1_BMSbin.NEC_WM_day_QC = MAug1_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
MAug1_BMSbin.SDN_day(MAug1_inight) = NaN;
MAug1_BMSbin.PAR_day (MAug1_inight) = NaN;

MAug1_BMSbin.NEP_day(MAug1_inight) = NaN;
MAug1_BMSbin.NEC_day(MAug1_inight) = NaN;
MAug1_BMSbin.NEP_day_QC(MAug1_inight) = NaN;
MAug1_BMSbin.NEC_day_QC(MAug1_inight) = NaN;

MAug1_BMSbin.dDOXY_day_QC(MAug1_inight) = NaN;      %DO Gradient
MAug1_BMSbin.dTA_day_QC(MAug1_inight) = NaN;        %TA Gradient

MAug1_BMSbin.NEP_WM_day(MAug1_inight) = NaN;
MAug1_BMSbin.NEC_WM_day(MAug1_inight) = NaN;
MAug1_BMSbin.NEP_WM_day_QC(MAug1_inight) = NaN;
MAug1_BMSbin.NEC_WM_day_QC(MAug1_inight) = NaN;

%Plot to check only nighttime points removed
figure 
hold on
scatter(MAug1_BMSbin.SDN, MAug1_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MAug1_BMSbin.SDN_day, MAug1_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
MAug1_BMSbin.NEP_day(isnan(MAug1_BMSbin.NEP_day))=[];
MAug1_BMSbin.NEC_day(isnan(MAug1_BMSbin.NEC_day))=[];
MAug1_BMSbin.NEP_day_QC(isnan(MAug1_BMSbin.NEP_day_QC))=[];
MAug1_BMSbin.NEC_day_QC(isnan(MAug1_BMSbin.NEC_day_QC))=[];

MAug1_BMSbin.dDOXY_day_QC(isnan(MAug1_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
MAug1_BMSbin.dTA_day_QC(isnan(MAug1_BMSbin.dTA_day_QC))=[];       %TA Gradient

MAug1_BMSbin.NEP_WM_day(isnan(MAug1_BMSbin.NEP_WM_day))=[];
MAug1_BMSbin.NEC_WM_day(isnan(MAug1_BMSbin.NEC_WM_day))=[];
MAug1_BMSbin.NEP_WM_day_QC(isnan(MAug1_BMSbin.NEP_WM_day_QC))=[];
MAug1_BMSbin.NEC_WM_day_QC(isnan(MAug1_BMSbin.NEC_WM_day_QC))=[];

% create nighttime hours datasets
MAug1_BMSbin.SDN_night = MAug1_BMSbin.SDN;
MAug1_BMSbin.PAR_night = MAug1_BMSbin.PAR(1,:);

MAug1_BMSbin.NEP_night = MAug1_BMSbin.NEP;
MAug1_BMSbin.NEC_night = MAug1_BMSbin.NEC;
MAug1_BMSbin.NEP_night_QC = MAug1_BMSbin.NEP_QC;
MAug1_BMSbin.NEC_night_QC = MAug1_BMSbin.NEC_QC;

MAug1_BMSbin.dDOXY_night_QC = MAug1_BMSbin.dDOXY_QC;      %DO Gradient
MAug1_BMSbin.dTA_night_QC = MAug1_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
MAug1_BMSbin.SDN_night=MAug1_BMSbin.SDN_night(MAug1_inight);
MAug1_BMSbin.PAR_night=MAug1_BMSbin.PAR_night(MAug1_inight);

MAug1_BMSbin.NEP_night_QC=MAug1_BMSbin.NEP_night_QC(MAug1_inight);
MAug1_BMSbin.NEC_night_QC=MAug1_BMSbin.NEC_night_QC(MAug1_inight);

MAug1_BMSbin.dDOXY_night_QC=MAug1_BMSbin.dDOXY_night_QC(MAug1_inight);      %DO Gradient
MAug1_BMSbin.dTA_night_QC=MAug1_BMSbin.dTA_night_QC(MAug1_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(MAug1_BMSbin.SDN, MAug1_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MAug1_BMSbin.SDN_night, MAug1_BMSbin.PAR_night, 'r.'); % night plot

%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(MAug1_BMSbin.NEP_day,MAug1_BMSbin.NEC_day);
% MAug1_BMSbin.Reg_Line = m*MAug1_BMSbin.NEP_day + b;
% MAug1_BMSbin.Ratio = m;
% MAug1_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(MAug1_BMSbin.NEP_day,MAug1_BMSbin.NEC_day,'o')
% plot(MAug1_BMSbin.NEP_day,MAug1_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Aug 2020 Pre-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug1_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug1_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(MAug1_BMSbin.NEP_day_QC,MAug1_BMSbin.NEC_day_QC);
MAug1_BMSbin.Reg_Line_QC = m_QC*MAug1_BMSbin.NEP_day_QC + b_QC;
MAug1_BMSbin.Ratio_QC = m_QC;
MAug1_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(MAug1_BMSbin.NEP_day_QC,MAug1_BMSbin.NEC_day_QC,'o')
plot(MAug1_BMSbin.NEP_day_QC,MAug1_BMSbin.Reg_Line_QC,'r')
%ylim([-50 50])
%xlim([-25 25])
xlabel('NCP');
ylabel('NCC');
title('Marker 32 Aug 2020 Pre-Restoration NCC:NCP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug1_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug1_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(MAug1_BMSbin.NEP_WM_day,MAug1_BMSbin.NEC_WM_day);
% MAug1_BMSbin.Reg_Line_WM = m_WM*MAug1_BMSbin.NEP_WM_day + b_WM;
% MAug1_BMSbin.Ratio_WM = m_WM;
% MAug1_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(MAug1_BMSbin.NEP_WM_day,MAug1_BMSbin.NEC_WM_day,'o')
% plot(MAug1_BMSbin.NEP_WM_day,MAug1_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Aug 2020 Pre-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug1_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug1_BMSbin.R2_WM)
% 

[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(MAug1_BMSbin.NEP_WM_day_QC,MAug1_BMSbin.NEC_WM_day_QC);
MAug1_BMSbin.Reg_Line_WM_QC = m_WM_QC*MAug1_BMSbin.NEP_WM_day_QC + b_WM_QC;
MAug1_BMSbin.Ratio_WM_QC = m_WM_QC;
MAug1_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(MAug1_BMSbin.NEP_WM_day_QC,MAug1_BMSbin.NEC_WM_day_QC,'o')
% plot(MAug1_BMSbin.NEP_WM_day_QC,MAug1_BMSbin.Reg_Line_WM_QC,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Aug 2020 Pre-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug1_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug1_BMSbin.R2_WM_QC)

%% For NEC:NEP Regressions Using Gradients
close all 
clc

% multiply o2 gradient by -1 for O2 production
MAug1_BMSbin.dDOXY_Reg = -1.*MAug1_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
MAug1_BMSbin.dTA_Reg = 0.5.*MAug1_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(MAug1_BMSbin.dDOXY_day_QC);
figure 
hold on 
DOday = plot(1:Xlength, MAug1_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, MAug1_BMSbin.dDOXY_Reg);
xlabel('Aug1 Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Marker 32 Aug1 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, MAug1_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, MAug1_BMSbin.dTA_Reg);
xlabel('Aug1 Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Marker 32 Aug1 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(MAug1_BMSbin.dDOXY_Reg, MAug1_BMSbin.dTA_Reg);
MAug1_BMSbin.Reg_Line_G = m_G*MAug1_BMSbin.dDOXY_Reg + b_G;
MAug1_BMSbin.Ratio_G = m_G;
MAug1_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(MAug1_BMSbin.dDOXY_Reg,MAug1_BMSbin.dTA_Reg,'o')
plot(MAug1_BMSbin.dDOXY_Reg,MAug1_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Marker 32 Aug1 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug1_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug1_BMSbin.R2_G)

clc
disp('Finished with Aug 1 Metabolism Calculations');

save('MAug120_2.mat', 'MAug1_BMS', 'MAug1_BMSbin');

%% Plot Profiles 
% plot for profile within pump heights
% close all 
% for i =1:100 %length(MAug1_ADavg.SDN)
%     figure (i)
%     scatter(MAug1_ADavg.uv(1:108,i), MAug1_ADavg.bin_depth(1:108))
%     title(['Marker 32 Aug Velocity Profile Number ',num2str(i),])
%     xlabel('Velocity (m/s)');
%     ylabel('Height (m)');
% end


%% Subplots 
close all
clc

sgtitle('Marker 32 Aug1 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'08/05 12:00';'08/06 12:00';'08/07 12:00';'08/08 12:00';'08/09 12:00';...
    '08/10 12:00';'08/11 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
% legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(MAug1_BMSbin.SDN, zeros(size(MAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug1_good_Xrange, 'XTick', MAug1_tick, 'xticklabel', MAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug1 Days');
xticklabels({'08/05 12:00';'08/06 12:00';'08/07 12:00';'08/08 12:00';'08/09 12:00';...
    '08/10 12:00';'08/11 12:00'})
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MAug1_BMSbin.SDN, MAug1_BMSbin.NEC_QC, 24);
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
plot(MAug1_BMSbin.NEP_day_QC,MAug1_BMSbin.NEC_day_QC,'o')
plot(MAug1_BMSbin.NEP_day_QC,MAug1_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(MAug1_BMSbin.Ratio_QC,2);
str2 = num2str(MAug1_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(MAug1_BMSbin.dDOXY_Reg,MAug1_BMSbin.dTA_Reg,'o')
plot(MAug1_BMSbin.dDOXY_Reg,MAug1_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
title('dTA:dDO Ratio from Gradients');
str1 = num2str(MAug1_BMSbin.Ratio_G,2);
str2 = num2str(MAug1_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')



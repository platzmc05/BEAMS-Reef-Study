% July SeapHOx Data Analysis from Cudjoe Ledge Reef
% Michelle Platz - USF 
% 3/3/2021

% SeapHOx sensor deployed 6/30/2020 at 2pm EST
    % pump 1 height above benthos = 72 cm  
    % pump 2 height above benthos = 24 cm 
% ADP sensor deployed 6/30/2020 at 2pm EST
    % height from substrate to ADCP head = 21 cm

close all
clc
clear all
%height variables *** heights in meters! 
UJuly_z1 = 0.72;
UJuly_z2 = 0.24;
UJuly_ADheight = 0.21;
UJuly_ADbin_depth_1m = 1-UJuly_ADheight;
%% Initial look at data
% ***** create UJuly20_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
UJuly20_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('U_630BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
UJuly20_SPraw.DOXY = UJuly20_SPraw.O2SATPER.*calcO2sat(UJuly20_SPraw.MCAT_TC, UJuly20_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
UJuly20_SPraw.pHint_prelim = calc_dfet_pHint(UJuly20_SPraw.Vint, UJuly20_SPraw.DFET_TC, -0.4);

% ***** create UJuly20_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
UJuly20_SP = parse_to_pumpheights_ARM_2pump_Mar19(UJuly20_SPraw, [datenum('06-30-2020 15:00:00'), datenum('07-26-2020 8:59:30')]);

% Calculate Gradients 
UJuly20_SP = calc_TA_gradientV2(UJuly20_SP, 2369.19, [0.8:0.1:1.2], 1, 2); %old TA: 2304.87

% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

UJuly20_SP.dDOXY = UJuly20_SP.DOXY(1,:) - UJuly20_SP.DOXY(2,:); %Oxygen Gradient
UJuly20_SP.dpH = UJuly20_SP.pH(1,:) - UJuly20_SP.pH(2,:); %pH Gradient 
UJuly20_SP.dTA = UJuly20_SP.TAtop - UJuly20_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
UJuly20_DateString = {'07/01/2020 12:00:00';'07/02/2020 12:00:00';'07/03/2020 12:00:00';'07/04/2020 12:00:00';...
    '07/05/2020 12:00:00';'07/06/2020 12:00:00';'07/07/2020 12:00:00';'07/08/2020 12:00:00';'07/09/2020 12:00:00';...
    '07/10/2020 12:00:00';'07/11/2020 12:00:00';'07/12/2020 12:00:00';'07/13/2020 12:00:00';'07/14/2020 12:00:00';...
    '07/15/2020 12:00:00';'07/16/2020 12:00:00';'07/17/2020 12:00:00';'07/18/2020 12:00:00';'07/19/2020 12:00:00';...
    '07/20/2020 12:00:00';'07/21/2020 12:00:00';'07/22/2020 12:00:00';'07/23/2020 12:00:00';'07/24/2020 12:00:00';...
    '07/25/2020 12:00:00'};
formatIn = 'mm/dd/yyyy HH:MM:SS';
UJuly20_tick = datenum(UJuly20_DateString,formatIn);

UJuly20_Xrange = [datenum('06-30-2020 14:00:00'), datenum('07-26-2020 8:59:30')];

figure
hold on; box on;
plot(UJuly20_SP.SDN, UJuly20_SP.dDOXY); %oxygen gradient 
plot(UJuly20_SP.SDN, UJuly20_SP.dTA); %TA gradient
plot(UJuly20_SP.SDN, zeros(size(UJuly20_SP.SDN))); %zero line
set(gca, 'xlim', UJuly20_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe July 2020 Unbinned DO and TA Gradients');

%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
UJuly_SiteChar.SDN = UJuly20_SPraw.SDN;
UJuly_SiteChar.TC = UJuly20_SPraw.OPT_TC;
UJuly_SiteChar.PAR = UJuly20_SPraw.PAR;
UJuly_SiteChar.PSAL = UJuly20_SPraw.PSAL;
UJuly_SiteChar.Pres = UJuly20_SPraw.Pres;

close all
figure 
hold on; 
plot(UJuly_SiteChar.SDN, UJuly_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
UJuly_SiteChar.SDN = UJuly20_SPraw.SDN(113:74268);
UJuly_SiteChar.TC = UJuly20_SPraw.OPT_TC(113:74268);
UJuly_SiteChar.PAR = UJuly20_SPraw.PAR(113:74268);
UJuly_SiteChar.PSAL = UJuly20_SPraw.PSAL(113:74268);
UJuly_SiteChar.Pres = UJuly20_SPraw.Pres(113:74268);

%extract full length of ADCP datafile  
UJuly_ADfull=aquadoppraw2mat('U_6_30', 70, [datenum('06-30-2020 14:00:00'), datenum('07-09-2020 05:03:00')]);

% add AD variables to Site Char 
UJuly_SiteChar.AD_SDN = UJuly_ADfull.SDN;
UJuly_SiteChar.AD_Pres = UJuly_ADfull.Pres;
UJuly_SiteChar.AD_TC = UJuly_ADfull.TC;
UJuly_SiteChar.bin_depth = UJuly_ADfull.bin_depth;
UJuly_SiteChar.u = UJuly_ADfull.u;
UJuly_SiteChar.v = UJuly_ADfull.v;
UJuly_SiteChar.w = UJuly_ADfull.w;
UJuly_SiteChar.uv = UJuly_ADfull.uv;
UJuly_SiteChar.direction = UJuly_ADfull.direction;

%Plot to see when surface interval observations are
close all
figure 
hold on; 
plot(UJuly_SiteChar.AD_SDN, UJuly_SiteChar.AD_TC)

%clip ends of data to remove surfave interval observations 
UJuly_SiteChar.AD_SDN = UJuly_ADfull.SDN(108:end);
UJuly_SiteChar.AD_Pres = UJuly_ADfull.Pres(108:end);
UJuly_SiteChar.AD_TC = UJuly_ADfull.TC(108:end);
UJuly_SiteChar.bin_depth = UJuly_ADfull.bin_depth;
UJuly_SiteChar.u = UJuly_ADfull.u(:,108:end);
UJuly_SiteChar.v = UJuly_ADfull.v(:,108:end);
UJuly_SiteChar.w = UJuly_ADfull.w(:,108:end);
UJuly_SiteChar.uv = UJuly_ADfull.uv(:,108:end);
UJuly_SiteChar.direction = UJuly_ADfull.direction(:,108:end);

% find U0 
UJuly_z1 = 0.72;
UJuly_z2 = 0.24;
UJuly_ADheight = 0.21;
UJuly_ADbin_depth_1m = 1-UJuly_ADheight;
UJuly_i1m = find(UJuly_SiteChar.bin_depth==UJuly_ADbin_depth_1m);
UJuly_SiteChar.U0 = UJuly_SiteChar.uv(UJuly_i1m,:);

% save data in separate datastructure
%save('UJuly20_SiteChar_2.mat', 'UJuly_SiteChar')


%% Constrain Xrange from graph results and extract good gradient data - 
close all 

UJuly20_good_Xrange = [datenum('06-30-2020 16:15:00'), datenum('04-Jul-2020 02:45:00')];

UJuly20_DateString = {'07/01/2020 12:00:00';'07/02/2020 12:00:00';'07/03/2020 12:00:00';'07/04/2020 12:00:00';...
    '07/05/2020 12:00:00';'07/06/2020 12:00:00';'07/07/2020 12:00:00';'07/08/2020 12:00:00';'07/09/2020 12:00:00';...
    '07/10/2020 12:00:00';'07/11/2020 12:00:00';'07/12/2020 12:00:00';'07/13/2020 12:00:00';'07/14/2020 12:00:00';...
    '07/15/2020 12:00:00';'07/16/2020 12:00:00';'07/17/2020 12:00:00';'07/18/2020 12:00:00';'07/19/2020 12:00:00';...
    '07/20/2020 12:00:00';'07/21/2020 12:00:00';'07/22/2020 12:00:00';'07/23/2020 12:00:00';'07/24/2020 12:00:00';...
    '07/25/2020 12:00:00'};
formatIn = 'mm/dd/yyyy HH:MM:SS';
UJuly20_tick = datenum(UJuly20_DateString,formatIn);

% plot to check range is correct
figure
hold on; box on;
plot(UJuly20_SP.SDN, UJuly20_SP.dDOXY); %oxygen gradient 
plot(UJuly20_SP.SDN, UJuly20_SP.dTA); %TA gradient
plot(UJuly20_SP.SDN, zeros(size(UJuly20_SP.SDN))); %zero line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe July 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all

U_JulyBMS_idx_start = find(UJuly20_SP.SDN==datenum('06-30-2020 16:15:00'));
U_JulyBMS_idx_end = find(UJuly20_SP.SDN==datenum('07-04-2020 02:45:00'));

% Create new data vectors of just the good data
U_JulyBMS_good_data = U_JulyBMS_idx_start:U_JulyBMS_idx_end;
Initial_data_points = length(U_JulyBMS_good_data)

% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(UJuly20_SP);
for v = 1:length(vars)
    UJuly20_SP.(vars{v}) = (UJuly20_SP.(vars{v})(:,U_JulyBMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: UJuly_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
UJuly_AD=aquadoppraw2mat('U_6_30', 70, [datenum('06-30-2020 16:15:00'), datenum('04-Jul-2020 03:00:00')]);

%averages data to the middle of the minute interval spacified 
UJuly_ADavg = average_aquadopp(UJuly_AD, 15);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.72m (pump 1) and 0.24m (pump 2) 
% 0.21m from substrate to ACDP head - 
% adjusted height = 0.51m (bin 41) and 0.03m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period - 35.625
clc
% already removed data outside data frame so can take average of whole set
UJuly_Sal_est = mean(UJuly20_SP.PSAL(1,3:end));

[UJuly_ADavg] = ustar_from_aquadopp2(UJuly_ADavg,[0.51 0.11], UJuly_Sal_est); %bins adjusted 


%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[UJuly_ADavg] = ustar_McGillis_Method(UJuly_ADavg, 0.72, 0.24, 41, 1);
close all
figure
hold on 
ustar_plot = plot(UJuly_ADavg.SDN, UJuly_ADavg.ustar, 'r');
ustar_WM_plot = plot(UJuly_ADavg.SDN, UJuly_ADavg.ustar_WM);
plot(UJuly_ADavg.SDN, zeros(size(UJuly_ADavg.SDN)),'k');
set(gca, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlim([UJuly_ADavg.SDN(1) UJuly_ADavg.SDN(end)])
xlabel('July Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Cudjoe July 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: UJuly_BMS *****

ADavg_vars = fieldnames(UJuly_ADavg);
for v = 1:length(ADavg_vars)
    UJuly_BMS.(ADavg_vars{v}) = (UJuly_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(UJuly20_SP);
for v = 1:length(SP_vars)
    UJuly_BMS.(SP_vars{v}) = (UJuly20_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(UJuly_BMS.SDN)


%% %% *************** Calculate Fluxes ****************
% actual pump heights  = 0.72m (pump 1) and 0.24m (pump 2) 
% 0.21m from substrate to ACDP head in July at U 
% adjusted height = 0.51m (bin 41) and 0.03m (too shallow) (bin 1)  - bins with which to calc ustar 
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[UJuly_BMS] = calc_NCC_3(UJuly_BMS,[0.51 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(UJuly_BMS.DOXY(1,:))
[UJuly_BMS] = calc_NCP_3(UJuly_BMS, [0.51 0.11],C1guess); 

% Plot NCP and NCC
close all
figure
hold on; box on; 
NEPplot = plot(UJuly_BMS.SDN, UJuly_BMS.NEP);
NECplot = plot(UJuly_BMS.SDN, UJuly_BMS.NEC);
plot(UJuly_BMS.SDN, zeros(size(UJuly_BMS.SDN)),'k');
set(gca, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UJuly_BMS.SDN(1) UJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Cudjoe July 2020 Fluxes');

%McGillis method flux calculations 
[UJuly_BMS] = calc_NCP_McGillis_Method(UJuly_BMS, 0.72, 0.24, UJuly_Sal_est);
[UJuly_BMS] = calc_NCC_McGillis_Method(UJuly_BMS, 0.72, 0.24, UJuly_Sal_est);

% Plot NCP and NCC WM
% close all
% figure
% hold on; box on; 
% NEPplot = plot(UJuly_BMS.SDN, UJuly_BMS.NEP_WM);
% NECplot = plot(UJuly_BMS.SDN, UJuly_BMS.NEC_WM);
% plot(UJuly_BMS.SDN, zeros(size(UJuly_BMS.SDN)),'k');
% set(gca, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlim([UJuly_BMS.SDN(1) UJuly_BMS.SDN(end)])
% xlabel('July Days');
% ylabel('NCP or NCC [mmol/m2/hr]');
% legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
% title('Cudjoe July 2020 WM Fluxes');


%Compare Flux_fit vs WM Plots 

% NCP Plot  
% close all
figure
hold on; box on; 
NEPplot = plot(UJuly_BMS.SDN, UJuly_BMS.NEP);
NEPplotWM = plot(UJuly_BMS.SDN, UJuly_BMS.NEP_WM);
plot(UJuly_BMS.SDN, zeros(size(UJuly_BMS.SDN)),'k');
set(gca, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UJuly_BMS.SDN(1) UJuly_BMS.SDN(end)])
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Cudjoe July 2020 NEP');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(UJuly_BMS.SDN, UJuly_BMS.NEC);
NECplotWM = plot(UJuly_BMS.SDN, UJuly_BMS.NEC_WM);
plot(UJuly_BMS.SDN, zeros(size(UJuly_BMS.SDN)),'k');
set(gca, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UJuly_BMS.SDN(1) UJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Cudjoe July 2020 NEC');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
UJuly_BMS.DOXYstd = std(UJuly_BMS.DOXY);

% get DOXY std
UJuly_idoxystd = find(UJuly_BMS.DOXYstd > 2);% 58.4% of data is greater than 0.8 stdev
UJuly_ihighdoxystd = [];
for i = 1:length(UJuly_idoxystd)
    
    UJuly_ihighdoxystd = vertcat(UJuly_ihighdoxystd,[UJuly_idoxystd(i)-1:1:UJuly_idoxystd(i)+1]');
end
% get unique IDs
UJuly_ihighdoxystd = unique(UJuly_ihighdoxystd);
% remove 0's and out of index values
UJuly_ihighdoxystd(UJuly_ihighdoxystd==0) = [];
UJuly_ihighdoxystd(UJuly_ihighdoxystd> length(UJuly_BMS.SDN)) = [];

% make it into index
trex = false(size(UJuly_BMS.SDN));
trex(UJuly_ihighdoxystd) = true;
UJuly_ihighdoxystd = trex;
clear trex;

UJuly_BMS.NEP_QC = UJuly_BMS.NEP;
UJuly_BMS.NEC_QC = UJuly_BMS.NEC;
UJuly_BMS.dDOXY_QC = UJuly_BMS.dDOXY; %DO gradient
UJuly_BMS.dTA_QC = UJuly_BMS.dTA;     %TA gradient
UJuly_BMS.NEP_WM_QC = UJuly_BMS.NEP_WM;
UJuly_BMS.NEC_WM_QC = UJuly_BMS.NEC_WM;

% set observations when DOXYstd>0.8 to NaN
UJuly_BMS.NEP_QC(UJuly_ihighdoxystd) = NaN;
UJuly_BMS.NEC_QC(:,UJuly_ihighdoxystd) = NaN;
UJuly_BMS.dDOXY_QC(UJuly_ihighdoxystd) = NaN; %DO gradient
UJuly_BMS.dTA_QC(:,UJuly_ihighdoxystd) = NaN;   %TA gradient
UJuly_BMS.NEP_WM_QC(UJuly_ihighdoxystd) = NaN;
UJuly_BMS.NEC_WM_QC(:,UJuly_ihighdoxystd) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(UJuly_BMS.SDN, UJuly_BMS.NEP, 'k');
NEPplotQC = plot(UJuly_BMS.SDN, UJuly_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(UJuly_BMS.SDN, zeros(size(UJuly_BMS.SDN)),'k');
set(gca, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UJuly_BMS.SDN(1) UJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Cudjoe July 2020 Fluxes');


%WM Plot
% figure
% hold on; box on;
% NEPplot = plot(UJuly_BMS.SDN, UJuly_BMS.NEP_WM, 'k');
% NEPplotQC = plot(UJuly_BMS.SDN, UJuly_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
% plot(UJuly_BMS.SDN, zeros(size(UJuly_BMS.SDN)),'k');
% set(gca, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlim([UJuly_BMS.SDN(1) UJuly_BMS.SDN(end)])
% xlabel('July Days');
% ylabel('NCP [mmol/m2/hr]');
% legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
% title('Cudjoe July 2020 Fluxes');


%% Bin data to hourly intervals
clc

X = floor(nanmin(UJuly_BMS.SDN)):1/24:ceil(nanmax(UJuly_BMS.SDN));

UJuly_BMSbin.SDN = X;
% variables from 2 differnet heights
UJuly_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.DOXY(1,:), X);
UJuly_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.DOXY(2,:), X);

UJuly_BMSbin.pH(1,:)       = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.pH(1,:), X);
UJuly_BMSbin.pH(2,:)       = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.pH(2,:), X);

UJuly_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.PSAL(1,:), X);
UJuly_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.PSAL(2,:), X);

UJuly_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.O2SATPER(1,:), X);
UJuly_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.O2SATPER(2,:), X);

UJuly_BMSbin.Pres(1,:)     = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.Pres(1,:), X);
UJuly_BMSbin.Pres(2,:)     = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.Pres(2,:), X);

UJuly_BMSbin.DENS(1,:)     = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.DENS(1,:), X);
UJuly_BMSbin.DENS(2,:)     = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.DENS(2,:), X);

UJuly_BMSbin.PAR(1,:)      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.PAR(1,:), X);
UJuly_BMSbin.PAR(2,:)      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.PAR(2,:), X);

UJuly_BMSbin.bin_depth     = UJuly_BMS.bin_depth;

for i = 1:108
    UJuly_BMSbin.uv(i,:)   = bin_data_to_X_GF(UJuly_BMS.SDN,UJuly_BMS.uv(i,:), X);
end

% bin data hourly. Vector variables 
UJuly_BMSbin.PRES  = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.Pres, X);
UJuly_BMSbin.U0       = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.U0, X);
UJuly_BMSbin.DIR      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.direction, X);
UJuly_BMSbin.ustar    = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.ustar, X);
UJuly_BMSbin.ustar_rm = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.ustar_runmean, X);
UJuly_BMSbin.dTA      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.dTA, X);
UJuly_BMSbin.dTA_QC   = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.dTA_QC, X);
UJuly_BMSbin.dpH      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.dpH, X);
UJuly_BMSbin.dDOXY    = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.dDOXY, X);
UJuly_BMSbin.dDOXY_QC = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.dDOXY_QC, X);
UJuly_BMSbin.NEP      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.NEP, X);
UJuly_BMSbin.NEP_QC   = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.NEP_QC, X);
UJuly_BMSbin.NEC      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.NEC, X);
UJuly_BMSbin.NEC_QC   = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.NEC_QC, X);
UJuly_BMSbin.NEP_WM      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.NEP_WM, X);
UJuly_BMSbin.NEP_WM_QC   = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.NEP_WM_QC, X);
UJuly_BMSbin.NEC_WM      = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.NEC_WM, X);
UJuly_BMSbin.NEC_WM_QC   = bin_data_to_X_GF(UJuly_BMS.SDN, UJuly_BMS.NEC_WM_QC, X);

%% Plot Binned Fluxes 
close all
clc
% % figure
% % hold on; box on;
% % NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP, 'b'); 
% % NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC, 'r-'); 
% % plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
% % set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
% % datetick('x', 'dd', 'keeplimits', 'keepticks');
% % xlabel('July Days');
% % ylabel('NEP or NEC');
% % legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% % title('Cudjoe July 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_QC, 'r-'); 
plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe July 2020 Hourly Binned Fluxes QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_WM, 'r-'); 
% plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe July 2020 Hourly Binned Fluxes WM');

% clc
% figure
% hold on; box on;
% NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_WM_QC, 'r-'); 
% plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe July 2020 Hourly Binned Fluxes WM QC');


%% Plot Binned Gradients 
% close all
figure
hold on; box on;
DOplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
% DOplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dDOXY, 'c'); 
% TAplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dTA, 'k-'); 
plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe July 2020 Binned Gradiets');


%% ***QC*** Remove sections when velocity is too slow
ibad = UJuly_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow

% UJuly_BMSbin.NEP(ibad) = NaN;
% UJuly_BMSbin.NEC(ibad) = NaN;
UJuly_BMSbin.NEP_QC(ibad) = NaN;
UJuly_BMSbin.NEC_QC(ibad) = NaN;
UJuly_BMSbin.dDOXY_QC(ibad) = NaN;
UJuly_BMSbin.dTA_QC(ibad) = NaN;

% UJuly_BMSbin.NEP_WM(ibad) = NaN;
% UJuly_BMSbin.NEC_WM(ibad) = NaN;
UJuly_BMSbin.NEP_WM_QC(ibad) = NaN;
UJuly_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
clc
% figure
% hold on; box on;
% NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP, 'b'); 
% NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC, 'r-'); 
% plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe July 2020 Hourly Binned Fluxes');

close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe July 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_QC, 'r-'); 
plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe July 2020 Hourly Binned Fluxes Full QC');




% WM plots
% figure
% hold on; box on;
% NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_WM, 'r-'); 
% plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe July 2020 Hourly Binned WM Fluxes');

% close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_WM_QC, 'r-'); 
% plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe July 2020 Hourly Binned WM Fluxes Full QC');

%% Boxplots and Profile Plots 

% Fluxfit Calcs
figure
hold on; box on;
boxplot(UJuly_BMSbin.NEP_QC)
ylabel('DO')
title('July NEP Boxplots')

figure
hold on; box on;
boxplot(UJuly_BMSbin.NEC_QC)
ylabel('TA')
title('July NEC Boxplots')

figure
hold on; box on;
boxplot(UJuly_BMSbin.dDOXY_QC)
ylabel('DO')
title('July dDO Boxplots')

figure
hold on; box on;
boxplot(UJuly_BMSbin.dTA_QC)
ylabel('TA')
title('July dTA Boxplots')


% 78: +NEP outlier (value: 34), not a DO outlier - plot profile and SDN to make sure daytime
        % 03-Jul-2020 05:00:00 - 
        % point ruled an OUTLIER: dark observation 
        
% 36: +NEP outlier (value: 21.4), not a DO outlier - plot profile and SDN to make sure daytime 
        % '01-Jul-2020 11:00:00 - daytime observation
        % profile good - NOT AN OUTLIER 
        
% 82: +NEC outlier (value: 14.3) and dTA outlier (value: 4.78)
        % point ruled an OUTLIER: extreme outliers for both gradients and flux

% 58: +NEC outlier (value: 14.4) and dTA outlier (value: 2.61)
        % '02-Jul-2020 09:00:00' 
        % profile fast but good
        % point an OUTLIER - double outlier

% 60: +dTA outlier (value: 3.13)
        % QC violated very next point
        % '02-Jul-2020 11:00:00'
        % bad profile - OUTLIER 
        
% 21: +NEC outlier, not a dTA outlier - plot profile
        % profile good 
        % not an outlier

% plot to check:  
datestr(UJuly_BMSbin.SDN(60))

close all 
for i = 60  %1:length(UJuly_BMSbin.SDN)
    figure (i)
    scatter(UJuly_BMSbin.uv(1:108,i), UJuly_BMSbin.bin_depth(1:108))
    title(['Cudjoe July Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end


%remove outliers: 78 82 58 60 

% 
% UJuly_BMSbin.NEP_QC(58) = NaN;
% UJuly_BMSbin.NEC_QC(58) = NaN;
% UJuly_BMSbin.dDOXY_QC(58) = NaN;
% UJuly_BMSbin.dTA_QC(58) = NaN;
% 
% UJuly_BMSbin.NEP_QC(60) = NaN;
% UJuly_BMSbin.NEC_QC(60) = NaN;
% UJuly_BMSbin.dDOXY_QC(60) = NaN;
% UJuly_BMSbin.dTA_QC(60) = NaN;
% 
% UJuly_BMSbin.NEP_QC(78) = NaN;
% UJuly_BMSbin.NEC_QC(78) = NaN;
% UJuly_BMSbin.dDOXY_QC(78) = NaN;
% UJuly_BMSbin.dTA_QC(78) = NaN;
% 
% UJuly_BMSbin.NEP_QC(82) = NaN;
% UJuly_BMSbin.NEC_QC(82) = NaN;
% UJuly_BMSbin.dDOXY_QC(82) = NaN;
% UJuly_BMSbin.dTA_QC(82) = NaN;

% replot to see after outliers removed
close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe July 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_QC, 'r-'); 
plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe July 2020 Hourly Binned Fluxes Full QC');
%% Plot Diel Curves

% nepdbin = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP, 24);
% necdbin = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin, 'bo', 'markersize', 3);
% plot(1:24, necdbin, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin,1), 'bo-');
% plot(1:24, nanmedian(necdbin,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('July 2020');
% xlabel('hour of day');
close all 
nepdbin_QC = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_QC, 24);
figure
hold on; box on;
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, nepdbin_QC, 'bo', 'markersize', 3);
plot(1:24, necdbin_QC, 'ro', 'markersize', 3);
plot(1:24, nanmedian(nepdbin_QC,1), 'bo-');
plot(1:24, nanmedian(necdbin_QC,1), 'ro-')
ylabel(['NEP or \color{red}NEC']);
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
set(gca, 'XGrid', 'on');
xlim([0 24]);
title('July 2020 Fluxfit');
xlabel('hour of day');

% % WM data
% nepdbin_WM = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_WM, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('July 2020');
% xlabel('hour of day');

% nepdbin_WM_QC = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_WM_QC, 24);
% necdbin_WM_QC = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_WM_QC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM_QC, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM_QC, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM_QC,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM_QC,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('July 2020 WM');
% xlabel('hour of day');

%% Depth-NEP plots 

scatter(UJuly_BMSbin.PRES, UJuly_BMSbin.NEP_QC)


%% Extract Daytime data for Ratios
clc
% Extract daytime data using UJuly_BMSbin.PAR
UJuly_inight = UJuly_BMSbin.PAR(1,:) < 1; %find all nightime datapoints 

%create new arrays for daytime data
UJuly_BMSbin.SDN_day = UJuly_BMSbin.SDN;
UJuly_BMSbin.PAR_day = UJuly_BMSbin.PAR(1,:);

UJuly_BMSbin.NEP_day = UJuly_BMSbin.NEP;
UJuly_BMSbin.NEC_day = UJuly_BMSbin.NEC;
UJuly_BMSbin.NEP_day_QC = UJuly_BMSbin.NEP_QC;
UJuly_BMSbin.NEC_day_QC = UJuly_BMSbin.NEC_QC;

UJuly_BMSbin.dDOXY_day_QC = UJuly_BMSbin.dDOXY_QC;      %DO Gradient
UJuly_BMSbin.dTA_day_QC = UJuly_BMSbin.dTA_QC;          %TA Gradient

UJuly_BMSbin.NEP_WM_day = UJuly_BMSbin.NEP_WM;
UJuly_BMSbin.NEC_WM_day = UJuly_BMSbin.NEC_WM;
UJuly_BMSbin.NEP_WM_day_QC = UJuly_BMSbin.NEP_WM_QC;
UJuly_BMSbin.NEC_WM_day_QC = UJuly_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
UJuly_BMSbin.SDN_day(UJuly_inight) = NaN;
UJuly_BMSbin.PAR_day (UJuly_inight) = NaN;

UJuly_BMSbin.NEP_day(UJuly_inight) = NaN;
UJuly_BMSbin.NEC_day(UJuly_inight) = NaN;
UJuly_BMSbin.NEP_day_QC(UJuly_inight) = NaN;
UJuly_BMSbin.NEC_day_QC(UJuly_inight) = NaN;

UJuly_BMSbin.dDOXY_day_QC(UJuly_inight) = NaN;      %DO Gradient
UJuly_BMSbin.dTA_day_QC(UJuly_inight) = NaN;        %TA Gradient

UJuly_BMSbin.NEP_WM_day(UJuly_inight) = NaN;
UJuly_BMSbin.NEC_WM_day(UJuly_inight) = NaN;
UJuly_BMSbin.NEP_WM_day_QC(UJuly_inight) = NaN;
UJuly_BMSbin.NEC_WM_day_QC(UJuly_inight) = NaN;

%Plot to check only nighttime points removed
figure 
hold on
scatter(UJuly_BMSbin.SDN, UJuly_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(UJuly_BMSbin.SDN_day, UJuly_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from daytime fluxes
UJuly_BMSbin.NEP_day(isnan(UJuly_BMSbin.NEP_day))=[];
UJuly_BMSbin.NEC_day(isnan(UJuly_BMSbin.NEC_day))=[];
UJuly_BMSbin.NEP_day_QC(isnan(UJuly_BMSbin.NEP_day_QC))=[];
UJuly_BMSbin.NEC_day_QC(isnan(UJuly_BMSbin.NEC_day_QC))=[];

UJuly_BMSbin.dDOXY_day_QC(isnan(UJuly_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
UJuly_BMSbin.dTA_day_QC(isnan(UJuly_BMSbin.dTA_day_QC))=[];       %TA Gradient

UJuly_BMSbin.NEP_WM_day(isnan(UJuly_BMSbin.NEP_WM_day))=[];
UJuly_BMSbin.NEC_WM_day(isnan(UJuly_BMSbin.NEC_WM_day))=[];
UJuly_BMSbin.NEP_WM_day_QC(isnan(UJuly_BMSbin.NEP_WM_day_QC))=[];
UJuly_BMSbin.NEC_WM_day_QC(isnan(UJuly_BMSbin.NEC_WM_day_QC))=[];

% create nighttime hours datasets
UJuly_BMSbin.SDN_night = UJuly_BMSbin.SDN;
UJuly_BMSbin.PAR_night = UJuly_BMSbin.PAR(1,:);

UJuly_BMSbin.NEP_night = UJuly_BMSbin.NEP;
UJuly_BMSbin.NEC_night = UJuly_BMSbin.NEC;
UJuly_BMSbin.NEP_night_QC = UJuly_BMSbin.NEP_QC;
UJuly_BMSbin.NEC_night_QC = UJuly_BMSbin.NEC_QC;

UJuly_BMSbin.dDOXY_night_QC = UJuly_BMSbin.dDOXY_QC;      %DO Gradient
UJuly_BMSbin.dTA_night_QC = UJuly_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
UJuly_BMSbin.SDN_night=UJuly_BMSbin.SDN_night(UJuly_inight);
UJuly_BMSbin.PAR_night=UJuly_BMSbin.PAR_night(UJuly_inight);

UJuly_BMSbin.NEP_night_QC=UJuly_BMSbin.NEP_night_QC(UJuly_inight);
UJuly_BMSbin.NEC_night_QC=UJuly_BMSbin.NEC_night_QC(UJuly_inight);

UJuly_BMSbin.dDOXY_night_QC=UJuly_BMSbin.dDOXY_night_QC(UJuly_inight);      %DO Gradient
UJuly_BMSbin.dTA_night_QC=UJuly_BMSbin.dTA_night_QC(UJuly_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(UJuly_BMSbin.SDN, UJuly_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(UJuly_BMSbin.SDN_night, UJuly_BMSbin.PAR_night, 'r.'); % night plot

%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(UJuly_BMSbin.NEP_day,UJuly_BMSbin.NEC_day);
% UJuly_BMSbin.Reg_Line = m*UJuly_BMSbin.NEP_day + b;
% UJuly_BMSbin.Ratio = m;
% UJuly_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(UJuly_BMSbin.NEP_day,UJuly_BMSbin.NEC_day,'o')
% plot(UJuly_BMSbin.NEP_day,UJuly_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe July 2020 Pre-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UJuly_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UJuly_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(UJuly_BMSbin.NEP_day_QC,UJuly_BMSbin.NEC_day_QC);
UJuly_BMSbin.Reg_Line_QC = m_QC*UJuly_BMSbin.NEP_day_QC + b_QC;
UJuly_BMSbin.Ratio_QC = m_QC;
UJuly_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(UJuly_BMSbin.NEP_day_QC,UJuly_BMSbin.NEC_day_QC,'o')
plot(UJuly_BMSbin.NEP_day_QC,UJuly_BMSbin.Reg_Line_QC,'r')
%ylim([-50 50])
%xlim([-25 25])
xlabel('NCP');
ylabel('NCC');
title('Cudjoe July 2020 Pre-Restoration NCC:NCP Ratio Fluxfit Full QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UJuly_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UJuly_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(UJuly_BMSbin.NEP_WM_day,UJuly_BMSbin.NEC_WM_day);
% UJuly_BMSbin.Reg_Line_WM = m_WM*UJuly_BMSbin.NEP_WM_day + b_WM;
% UJuly_BMSbin.Ratio_WM = m_WM;
% UJuly_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(UJuly_BMSbin.NEP_WM_day,UJuly_BMSbin.NEC_WM_day,'o')
% plot(UJuly_BMSbin.NEP_WM_day,UJuly_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe July 2020 Pre-Restoration NCC:NCP Ratio WM');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UJuly_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UJuly_BMSbin.R2_WM)


[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(UJuly_BMSbin.NEP_WM_day_QC,UJuly_BMSbin.NEC_WM_day_QC);
UJuly_BMSbin.Reg_Line_WM_QC = m_WM_QC*UJuly_BMSbin.NEP_WM_day_QC + b_WM_QC;
UJuly_BMSbin.Ratio_WM_QC = m_WM_QC;
UJuly_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(UJuly_BMSbin.NEP_WM_day_QC,UJuly_BMSbin.NEC_WM_day_QC,'o')
% plot(UJuly_BMSbin.NEP_WM_day_QC,UJuly_BMSbin.Reg_Line_WM_QC,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe July 2020 Pre-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UJuly_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UJuly_BMSbin.R2_WM_QC)

%% For NEC:NEP Regressions Using Gradients
close all 
clc

% multiply o2 gradient by -1 for O2 production
UJuly_BMSbin.dDOXY_Reg = -1.*UJuly_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
UJuly_BMSbin.dTA_Reg = 0.5.*UJuly_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
figure 
hold on 
DOday = plot(1:36, UJuly_BMSbin.dDOXY_day_QC);
DOreg = plot(1:36, UJuly_BMSbin.dDOXY_Reg);
xlabel('July Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Cudjoe July 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:36, UJuly_BMSbin.dTA_day_QC);
DOreg = plot(1:36, UJuly_BMSbin.dTA_Reg);
xlabel('July Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Cudjoe July 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(UJuly_BMSbin.dDOXY_Reg, UJuly_BMSbin.dTA_Reg);
UJuly_BMSbin.Reg_Line_G = m_G*UJuly_BMSbin.dDOXY_Reg + b_G;
UJuly_BMSbin.Ratio_G = m_G;
UJuly_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(UJuly_BMSbin.dDOXY_Reg,UJuly_BMSbin.dTA_Reg,'o')
plot(UJuly_BMSbin.dDOXY_Reg,UJuly_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Cudjoe July 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UJuly_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UJuly_BMSbin.R2_G)

clc
disp('Finished with July Metabolism Calculations');

save('UJuly20_2.mat', 'UJuly_BMS', 'UJuly_BMSbin');



%% Subplots 
close all
clc

sgtitle('Cudjoe July 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
% DOplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dDOXY, 'c'); 
% TAplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.dTA, 'k-'); 
plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
% legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(UJuly_BMSbin.SDN, zeros(size(UJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UJuly20_good_Xrange, 'XTick', UJuly20_tick, 'xticklabel', UJuly20_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/01 12:00';'07/02 12:00';'07/03 12:00';'07/04 12:00'})
xlabel('July Days');
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(UJuly_BMSbin.SDN, UJuly_BMSbin.NEC_QC, 24);
plot(1:24, zeros(size(1:24)), 'k:');
plot(1:24, nepdbin_QC, 'bo', 'markersize', 3);
plot(1:24, necdbin_QC, 'ro', 'markersize', 3);
plot(1:24, nanmedian(nepdbin_QC,1), 'bo-');
plot(1:24, nanmedian(necdbin_QC,1), 'ro-')
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
set(gca, 'XGrid', 'on');
xlim([0 24]);
ylabel(['\color{blue}NEP \color{black} or \color{red}NEC']);
title('Diel Plot');
xlabel('hour of day');

subplot(3,3,8); % Ratio Plot using fluxes 
hold on; box on;
plot(UJuly_BMSbin.NEP_day_QC,UJuly_BMSbin.NEC_day_QC,'o')
plot(UJuly_BMSbin.NEP_day_QC,UJuly_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(UJuly_BMSbin.Ratio_QC,2);
str2 = num2str(UJuly_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(UJuly_BMSbin.dDOXY_Reg,UJuly_BMSbin.dTA_Reg,'o')
plot(UJuly_BMSbin.dDOXY_Reg,UJuly_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
title('dTA:dDO Ratio from Gradients');
str1 = num2str(UJuly_BMSbin.Ratio_G,2);
str2 = num2str(UJuly_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


% Sept SeapHOx Data Analysis from Cudjoe Ledge Reef
% Michelle Platz - USF 
% 3/9/2021

% SeapHOx sensor deployed 9/2/2020 
    % SP datafile: 'U_902BMS.txt'
    % pump 1 height above benthos = 70 cm  
    % pump 2 height above benthos = 20 cm 
% ADP sensor deployed 9/02/2020 
    % ADCP datafile: 'U_9_02'
    % height from substrate to ADCP head = 18 cm

close all
clc
clear all
%% Initial look at data
% ***** create USept_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
USept_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('U_902BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
USept_SPraw.DOXY = USept_SPraw.O2SATPER.*calcO2sat(USept_SPraw.MCAT_TC, USept_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
USept_SPraw.pHint_prelim = calc_dfet_pHint(USept_SPraw.Vint, USept_SPraw.DFET_TC, -0.4);

% ***** create USept_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
USept_SP = parse_to_pumpheights_ARM_2pump_Mar19(USept_SPraw, [datenum('09-02-2020 16:00:00'), datenum('09-28-2020 10:59:30')]);

% Calculate Gradients 
USept_SP = calc_TA_gradientV2(USept_SP, 2369.19, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

USept_SP.dDOXY = USept_SP.DOXY(1,:) - USept_SP.DOXY(2,:); %Oxygen Gradient
USept_SP.dpH = USept_SP.pH(1,:) - USept_SP.pH(2,:); %pH Gradient 
USept_SP.dTA = USept_SP.TAtop - USept_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
USept_DateString = {'09/03/2020 12:00:00';'09/04/2020 12:00:00';'09/05/2020 12:00:00';...
    '09/06/2020 12:00:00';'09/07/2020 12:00:00';'09/08/2020 12:00:00';'09/09/2020 12:00:00';'09/10/2020 12:00:00';'09/11/2020 12:00:00';'09/12/2020 12:00:00';...
    '09/13/2020 12:00:00';'09/14/2020 12:00:00';'09/15/2020 12:00:00';'09/16/2020 12:00:00';'09/17/2020 12:00:00';'09/18/2020 12:00:00';'09/19/2020 12:00:00';...
    '09/20/2020 12:00:00';'09/21/2020 12:00:00';'09/22/2020 12:00:00';'09/23/2020 12:00:00';'09/24/2020 12:00:00';'09/25/2020 12:00:00';'09/26/2020 12:00:00';...
    '09/27/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
USept_tick = datenum(USept_DateString,formatIn);

USept_Xrange = [datenum('09-02-2020 16:00:00'), datenum('09-28-2020 10:59:30')];

figure
hold on; box on;
plot(USept_SP.SDN, USept_SP.dDOXY); %oxygen gradient 
plot(USept_SP.SDN, USept_SP.dTA); %TA gradient
plot(USept_SP.SDN, zeros(size(USept_SP.SDN))); %zero line
set(gca, 'xlim', USept_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe Sept 2020 Unbinned DO and TA Gradients');


%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
USept_SiteChar.SDN = USept_SPraw.SDN;
USept_SiteChar.TC = USept_SPraw.OPT_TC;
USept_SiteChar.PAR = USept_SPraw.PAR;
USept_SiteChar.PSAL = USept_SPraw.PSAL;
USept_SiteChar.Pres = USept_SPraw.Pres;

%Plot pressure data to see when surface interval observations are
close all
figure 
hold on; 
plot(USept_SiteChar.SDN, USept_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
USept_SiteChar.SDN = USept_SPraw.SDN(423:74615);
USept_SiteChar.TC = USept_SPraw.OPT_TC(423:74615);
USept_SiteChar.PAR = USept_SPraw.PAR(423:74615);
USept_SiteChar.PSAL = USept_SPraw.PSAL(423:74615);
USept_SiteChar.Pres = USept_SPraw.Pres(423:74615);

%extract full length of ADCP datafile  
USept_ADfull=aquadoppraw2mat('U_9_02', 70, [datenum('09-02-2020 16:00:00'), datenum('10-30-2020 05:57:30')]);

% add AD variables to Site Char 
USept_SiteChar.AD_SDN = USept_ADfull.SDN;
USept_SiteChar.AD_Pres = USept_ADfull.Pres;
USept_SiteChar.AD_TC = USept_ADfull.TC;
USept_SiteChar.bin_depth = USept_ADfull.bin_depth;
USept_SiteChar.u = USept_ADfull.u;
USept_SiteChar.v = USept_ADfull.v;
USept_SiteChar.w = USept_ADfull.w;
USept_SiteChar.uv = USept_ADfull.uv;
USept_SiteChar.direction = USept_ADfull.direction;

%Plot pressure data to see when surface interval observations are
close all
figure 
hold on; 
plot(USept_SiteChar.AD_SDN, USept_SiteChar.AD_Pres)

%clip ends of data to remove surfave interval observations 
USept_SiteChar.AD_SDN = USept_ADfull.SDN(71:end);
USept_SiteChar.AD_Pres = USept_ADfull.Pres(71:end);
USept_SiteChar.AD_TC = USept_ADfull.TC(71:end);
USept_SiteChar.bin_depth = USept_ADfull.bin_depth;
USept_SiteChar.u = USept_ADfull.u(:,71:end);
USept_SiteChar.v = USept_ADfull.v(:,71:end);
USept_SiteChar.w = USept_ADfull.w(:,71:end);
USept_SiteChar.uv = USept_ADfull.uv(:,71:end);
USept_SiteChar.direction = USept_ADfull.direction(:,71:end);

% find U0 
USept_z1 = 0.70;
USept_z2 = 0.20;
USept_ADheight = 0.18;
USept_ADbin_depth_1m = 1-(USept_ADheight);% = 0.82
USept_i1m = find(USept_SiteChar.bin_depth==(0.82));
USept_SiteChar.U0 = USept_SiteChar.uv(USept_i1m,:);

% save data in separate datastructure
save('USept20_SiteChar_2.mat', 'USept_SiteChar' )


%% Constrain Xrange from graph results and extract good gradient data - 
close all 

% ****** end date was originally 12th, but need to assess data drift before
% including these data ***************
USept_good_Xrange = [datenum('09-02-2020 18:00:00'), datenum('09-6-2020 12:00:00')]; 
% big weather event 9/12 pm - lose ADCP measurements for 24 hrs 
% pH sensor drifts after this point -- likely sedimentation issue

% plot to check range is correct
figure
hold on; box on;
plot(USept_SP.SDN, USept_SP.dDOXY); %oxygen gradient 
plot(USept_SP.SDN, USept_SP.dTA); %TA gradient
plot(USept_SP.SDN, zeros(size(USept_SP.SDN))); %zero line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe Sept 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all

U_Sept_BMS_idx_start = find(USept_SP.SDN==datenum('09-02-2020 18:00:00'))
U_Sept_BMS_idx_end = find(USept_SP.SDN==datenum('09-06-2020 12:00:00'))

% Create new data vectors of just the good data
U_Sept_BMS_good_data = U_Sept_BMS_idx_start:U_Sept_BMS_idx_end;
Initial_data_points = length(U_Sept_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(USept_SP);
for v = 1:length(vars)
    USept_SP.(vars{v}) = (USept_SP.(vars{v})(:,U_Sept_BMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: USept_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
USept_AD=aquadoppraw2mat('U_9_02', 70, [datenum('09-02-2020 18:00:00'), datenum('09-06-2020 12:15:00')]);

%averages data to the middle of the minute interval spacified 
USept_ADavg = average_aquadopp(USept_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.7m (pump 1) and 0.2m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.52m (bin 42) and 0.02m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
USept_Sal_est = mean(USept_SP.PSAL(1,3:end));

[USept_ADavg] = ustar_from_aquadopp2(USept_ADavg,[0.52 0.11], USept_Sal_est); %bins adjusted 

clc
%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[USept_ADavg] = ustar_McGillis_Method(USept_ADavg, 0.70, 0.20, 42, 1);

%compare ustar calculation methods
close all
figure
hold on 
ustar_plot = plot(USept_ADavg.SDN, USept_ADavg.ustar, 'r');
ustar_WM_plot = plot(USept_ADavg.SDN, USept_ADavg.ustar_WM);
plot(USept_ADavg.SDN, zeros(size(USept_ADavg.SDN)),'k');
set(gca, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([USept_ADavg.SDN(1) USept_ADavg.SDN(end)])
xlabel('Sept Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: USept_BMS *****

ADavg_vars = fieldnames(USept_ADavg);
for v = 1:length(ADavg_vars)
    USept_BMS.(ADavg_vars{v}) = (USept_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(USept_SP);
for v = 1:length(SP_vars)
    USept_BMS.(SP_vars{v}) = (USept_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(USept_BMS.SDN)
% min: 13-28-43-58 becuase SP was restarted in the field rather than on the
% minute
%% %% *************** Calculate Fluxes ****************

% actual pump heights  = 0.70m (pump 1) and 0.20m (pump 2) 
% 0.18m from substrate to ACDP head in Sept at U 
% adjusted height = 0.52 m (bin 42) and 0.02 m (too shallow) (bin 1)  
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[USept_BMS] = calc_NCC_3(USept_BMS,[0.52 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(USept_BMS.DOXY(1,:));
[USept_BMS] = calc_NCP_3(USept_BMS, [0.52 0.11],C1guess); %estimate C1 guess using median DOXY(1,:) value  

% Plot NCP and NCC
%close all
figure
hold on; box on; 
NEPplot = plot(USept_BMS.SDN, USept_BMS.NEP);
NECplot = plot(USept_BMS.SDN, USept_BMS.NEC);
plot(USept_BMS.SDN, zeros(size(USept_BMS.SDN)),'k');
set(gca, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([USept_BMS.SDN(1) USept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Fluxes');

%McGillis method flux calculations 
[USept_BMS] = calc_NCP_McGillis_Method(USept_BMS, 0.70, 0.20, USept_Sal_est);
[USept_BMS] = calc_NCC_McGillis_Method(USept_BMS, 0.70, 0.20, USept_Sal_est);

% Plot NCP and NCC
% close all
figure
hold on; box on; 
NEPplot = plot(USept_BMS.SDN, USept_BMS.NEP_WM);
NECplot = plot(USept_BMS.SDN, USept_BMS.NEC_WM);
plot(USept_BMS.SDN, zeros(size(USept_BMS.SDN)),'k');
set(gca, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([USept_BMS.SDN(1) USept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Cudjoe Sept 2020 WM Fluxes');

%% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(USept_BMS.SDN, USept_BMS.NEP);
NEPplotWM = plot(USept_BMS.SDN, USept_BMS.NEP_WM);
plot(USept_BMS.SDN, zeros(size(USept_BMS.SDN)),'k');
set(gca, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([USept_BMS.SDN(1) USept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(USept_BMS.SDN, USept_BMS.NEC);
NECplotWM = plot(USept_BMS.SDN, USept_BMS.NEC_WM);
plot(USept_BMS.SDN, zeros(size(USept_BMS.SDN)),'k');
set(gca, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([USept_BMS.SDN(1) USept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
USept_BMS.DOXYstd = std(USept_BMS.DOXY);

% get DOXY std
USept_idoxystd = find(USept_BMS.DOXYstd > 2);% 58.4% of data is greater than 0.8 stdev
USept_ihighdoxystd = [];
for i = 1:length(USept_idoxystd)
    
    USept_ihighdoxystd = vertcat(USept_ihighdoxystd,[USept_idoxystd(i)-1:1:USept_idoxystd(i)+1]');
end
% get unique IDs
USept_ihighdoxystd = unique(USept_ihighdoxystd);
% remove 0's and out of index values
USept_ihighdoxystd(USept_ihighdoxystd==0) = [];
USept_ihighdoxystd(USept_ihighdoxystd> length(USept_BMS.SDN)) = [];

% make it into index
trex = false(size(USept_BMS.SDN));
trex(USept_ihighdoxystd) = true;
USept_ihighdoxystd = trex;
clear trex;

%create _QC datasets to preserve original data
USept_BMS.NEP_QC = USept_BMS.NEP;
USept_BMS.NEC_QC = USept_BMS.NEC;
USept_BMS.dDOXY_QC = USept_BMS.dDOXY; %DO gradient
USept_BMS.dTA_QC = USept_BMS.dTA;     %TA gradient
USept_BMS.NEP_WM_QC = USept_BMS.NEP_WM;
USept_BMS.NEC_WM_QC = USept_BMS.NEC_WM;

% set observations when DOXYstd>X to NaN - only remove NaNs from QC datasets
USept_BMS.NEP_QC(USept_ihighdoxystd) = NaN;
USept_BMS.NEC_QC(:,USept_ihighdoxystd) = NaN;
USept_BMS.dDOXY_QC(USept_ihighdoxystd) = NaN; %DO gradient
USept_BMS.dTA_QC(:,USept_ihighdoxystd) = NaN;   %TA gradient
USept_BMS.NEP_WM_QC(USept_ihighdoxystd) = NaN;
USept_BMS.NEC_WM_QC(:,USept_ihighdoxystd) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(USept_BMS.SDN, USept_BMS.NEP, 'k');
NEPplotQC = plot(USept_BMS.SDN, USept_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(USept_BMS.SDN, zeros(size(USept_BMS.SDN)),'k');
set(gca, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([USept_BMS.SDN(1) USept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Fluxes');


%WM Plot
%close all
% figure
% hold on; box on;
% NEPplot = plot(USept_BMS.SDN, USept_BMS.NEP_WM, 'k');
% NEPplotQC = plot(USept_BMS.SDN, USept_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
% plot(USept_BMS.SDN, zeros(size(USept_BMS.SDN)),'k');
% set(gca, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlim([USept_BMS.SDN(1) USept_BMS.SDN(end)])
% xlabel('Sept Days');
% ylabel('NCP [mmol/m2/hr]');
% legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
% title('Cudjoe Sept 2020 Fluxes');
% 

%% Bin data to hourly intervals

X = floor(nanmin(USept_BMS.SDN)):1/24:ceil(nanmax(USept_BMS.SDN));

USept_BMSbin.SDN = X;
% variables from 2 differnet heights
USept_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.DOXY(1,:), X);
USept_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.DOXY(2,:), X);

USept_BMSbin.pH(1,:)       = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.pH(1,:), X);
USept_BMSbin.pH(2,:)       = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.pH(2,:), X);

USept_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.PSAL(1,:), X);
USept_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.PSAL(2,:), X);

USept_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.O2SATPER(1,:), X);
USept_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.O2SATPER(2,:), X);

USept_BMSbin.Pres(1,:)     = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.Pres(1,:), X);
USept_BMSbin.Pres(2,:)     = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.Pres(2,:), X);

USept_BMSbin.DENS(1,:)     = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.DENS(1,:), X);
USept_BMSbin.DENS(2,:)     = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.DENS(2,:), X);

USept_BMSbin.PAR(1,:)      = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.PAR(1,:), X);
USept_BMSbin.PAR(2,:)      = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.PAR(2,:), X);

USept_BMSbin.bin_depth     = USept_BMS.bin_depth;

for i = 1:108
    USept_BMSbin.uv(i,:)   = bin_data_to_X_GF(USept_BMS.SDN,USept_BMS.uv(i,:), X);
end

% bin data hourly. Vector variables 
USept_BMSbin.PRES  = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.Pres, X);
USept_BMSbin.U0       = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.U0, X);
USept_BMSbin.DIR      = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.direction, X);
USept_BMSbin.ustar    = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.ustar, X);
USept_BMSbin.ustar_rm = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.ustar_runmean, X);
USept_BMSbin.dTA      = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.dTA, X);
USept_BMSbin.dTA_QC   = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.dTA_QC, X);
USept_BMSbin.dpH      = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.dpH, X);
USept_BMSbin.dDOXY    = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.dDOXY, X);
USept_BMSbin.dDOXY_QC = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.dDOXY_QC, X);
USept_BMSbin.NEP      = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.NEP, X);
USept_BMSbin.NEP_QC   = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.NEP_QC, X);
USept_BMSbin.NEC      = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.NEC, X);
USept_BMSbin.NEC_QC   = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.NEC_QC, X);
USept_BMSbin.NEP_WM   = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.NEP_WM, X);
USept_BMSbin.NEP_WM_QC= bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.NEP_WM_QC, X);
USept_BMSbin.NEC_WM   = bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.NEC_WM, X);
USept_BMSbin.NEC_WM_QC= bin_data_to_X_GF(USept_BMS.SDN, USept_BMS.NEC_WM_QC, X);


%% Plot Binned Fluxes 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP, 'b'); 
% NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC, 'r-'); 
% plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Sept 2020 Hourly Binned Fluxes FluxFit');

% close all
clc
figure
hold on; box on;
NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP_QC, 'b'); 
NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC_QC, 'r-'); 
plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Hourly Binned Fluxes FluxFit QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC_WM, 'r-'); 
% plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Sept 2020 WM Original Hourly Binned Fluxes');

% close all
clc
figure
hold on; box on;
NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC_WM_QC, 'r-'); 
plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Sept 2020 WM Full QC Hourly Binned Fluxes');

%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(USept_BMSbin.SDN, USept_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(USept_BMSbin.SDN, USept_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
% DOplot = plot(USept_BMSbin.SDN, USept_BMSbin.dDOXY, 'c'); 
% TAplot = plot(USept_BMSbin.SDN, USept_BMSbin.dTA, 'k-'); 
plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Binned Gradiets');



%% ***QC*** Remove sections when velocity is too slow
ibad = USept_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow

% USept_BMSbin.NEP(ibad) = NaN;
% USept_BMSbin.NEC(ibad) = NaN;
USept_BMSbin.NEP_QC(ibad) = NaN;
USept_BMSbin.NEC_QC(ibad) = NaN;
USept_BMSbin.dDOXY_QC(ibad) = NaN;
USept_BMSbin.dTA_QC(ibad) = NaN;

% USept_BMSbin.NEP_WM(ibad) = NaN;
% USept_BMSbin.NEC_WM(ibad) = NaN;
USept_BMSbin.NEP_WM_QC(ibad) = NaN;
USept_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP, 'b'); 
% NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC, 'r-'); 
% plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Sept 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP_QC, 'b'); 
NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC_QC, 'r-'); 
plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC_WM, 'r-'); 
% plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Sept 2020 Hourly Binned WM Fluxes');

% clc
% figure
% hold on; box on;
% NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC_WM_QC, 'r-'); 
% plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Sept 2020 Hourly Binned WM Fluxes Full QC');


%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
close all 

figure
hold on; box on;
boxplot(USept_BMSbin.NEP_QC)
ylabel('NEP')
title('Sept NEP Boxplots')

figure
hold on; box on;
boxplot(USept_BMSbin.NEC_QC)
ylabel('NEC')
title('Sept NEC Boxplots')

figure
hold on; box on;
boxplot(USept_BMSbin.dDOXY_QC)
ylabel('DO')
title('Sept dDO Boxplots')

figure
hold on; box on;
boxplot(USept_BMSbin.dTA_QC)
ylabel('TA')
title('Sept dTA Boxplots')


% Outliers 

% 61: NEC and dTA outlier 
    % bad profile
    % OUTLIER 
% 51: NEC outlier (value: -9.67)
    % bad profile
    % OUTLIER 
% 50: NEC outlier (value: -7.29)
    % bad profile
    % OUTLIER 
% 28: NEC outlier (value: -6.68)
    % bad profile
    % OUTLIER 
    
length(USept_BMSbin.SDN)
datestr(USept_BMSbin.SDN(62))
USept_BMSbin.NEC_QC(62)
% Plot Profiles at outliers - 
close all 
for i = 28  %1:length(USept_BMSbin.SDN)
    figure (i)
    scatter(USept_BMSbin.uv(1:108,i), USept_BMSbin.bin_depth(1:108))
    title(['Cudjoe Sept Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end
%outliers to be removed:  28 51 61 50
% USept_BMSbin.NEP_QC(28) = NaN;
% USept_BMSbin.NEC_QC(28) = NaN;
% USept_BMSbin.dDOXY_QC(28) = NaN;
% USept_BMSbin.dTA_QC(28) = NaN;
% 
% USept_BMSbin.NEP_QC(50) = NaN;
% USept_BMSbin.NEC_QC(50) = NaN;
% USept_BMSbin.dDOXY_QC(50) = NaN;
% USept_BMSbin.dTA_QC(50) = NaN;
% 
% USept_BMSbin.NEP_QC(51) = NaN;
% USept_BMSbin.NEC_QC(51) = NaN;
% USept_BMSbin.dDOXY_QC(51) = NaN;
% USept_BMSbin.dTA_QC(51) = NaN;
% 
% USept_BMSbin.NEP_QC(61) = NaN;
% USept_BMSbin.NEC_QC(61) = NaN;
% USept_BMSbin.dDOXY_QC(61) = NaN;
% USept_BMSbin.dTA_QC(61) = NaN;

clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(USept_BMSbin.SDN, USept_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(USept_BMSbin.SDN, USept_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP_QC, 'b'); 
NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC_QC, 'r-'); 
plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Hourly Binned Fluxes Full QC');


%% Plot Diel Curves

% nepdbin = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEP, 24);
% necdbin = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin, 'bo', 'markersize', 3);
% plot(1:24, necdbin, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin,1), 'bo-');
% plot(1:24, nanmedian(necdbin,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('Sept Diel Plot 2020');
% xlabel('hour of day');

nepdbin_QC = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEC_QC, 24);
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
title('Sept Diel Plot 2020  QC');
xlabel('hour of day');

% WM data
% nepdbin_WM = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEC_WM, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('Sept 2020 WM');
% xlabel('hour of day');

nepdbin_WM_QC = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEC_WM_QC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM_QC, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM_QC, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM_QC,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM_QC,1), 'ro-')
% xlim([1 24]);
% ylabel(['NEP or \color{red}NEC']);
% title('Sept 2020 WM Full QC');
% xlabel('hour of day');



%% Extract Daytime data for Ratios
clc
% Extract daytime data using USept_BMSbin.PAR
USept_inight = USept_BMSbin.PAR(1,:) < 5; %find all nightime datapoints 

%create new arrays for daytime data
USept_BMSbin.SDN_day = USept_BMSbin.SDN;
USept_BMSbin.PAR_day = USept_BMSbin.PAR(1,:);

USept_BMSbin.NEP_day = USept_BMSbin.NEP;
USept_BMSbin.NEC_day = USept_BMSbin.NEC;
USept_BMSbin.NEP_day_QC = USept_BMSbin.NEP_QC;
USept_BMSbin.NEC_day_QC = USept_BMSbin.NEC_QC;

USept_BMSbin.dDOXY_day_QC = USept_BMSbin.dDOXY_QC;      %DO Gradient
USept_BMSbin.dTA_day_QC = USept_BMSbin.dTA_QC;          %TA Gradient

USept_BMSbin.NEP_WM_day = USept_BMSbin.NEP_WM;
USept_BMSbin.NEC_WM_day = USept_BMSbin.NEC_WM;
USept_BMSbin.NEP_WM_day_QC = USept_BMSbin.NEP_WM_QC;
USept_BMSbin.NEC_WM_day_QC = USept_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
USept_BMSbin.SDN_day(USept_inight) = NaN;
USept_BMSbin.PAR_day (USept_inight) = NaN;

USept_BMSbin.NEP_day(USept_inight) = NaN;
USept_BMSbin.NEC_day(USept_inight) = NaN;
USept_BMSbin.NEP_day_QC(USept_inight) = NaN;
USept_BMSbin.NEC_day_QC(USept_inight) = NaN;

USept_BMSbin.dDOXY_day_QC(USept_inight) = NaN;      %DO Gradient
USept_BMSbin.dTA_day_QC(USept_inight) = NaN;        %TA Gradient

USept_BMSbin.NEP_WM_day(USept_inight) = NaN;
USept_BMSbin.NEC_WM_day(USept_inight) = NaN;
USept_BMSbin.NEP_WM_day_QC(USept_inight) = NaN;
USept_BMSbin.NEC_WM_day_QC(USept_inight) = NaN;

%Plot to check only nighttime points removed
close all
figure 
hold on
scatter(USept_BMSbin.SDN, USept_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(USept_BMSbin.SDN_day, USept_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
USept_BMSbin.NEP_day(isnan(USept_BMSbin.NEP_day))=[];
USept_BMSbin.NEC_day(isnan(USept_BMSbin.NEC_day))=[];
USept_BMSbin.NEP_day_QC(isnan(USept_BMSbin.NEP_day_QC))=[];
USept_BMSbin.NEC_day_QC(isnan(USept_BMSbin.NEC_day_QC))=[];

USept_BMSbin.dDOXY_day_QC(isnan(USept_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
USept_BMSbin.dTA_day_QC(isnan(USept_BMSbin.dTA_day_QC))=[];       %TA Gradient

USept_BMSbin.NEP_WM_day(isnan(USept_BMSbin.NEP_WM_day))=[];
USept_BMSbin.NEC_WM_day(isnan(USept_BMSbin.NEC_WM_day))=[];
USept_BMSbin.NEP_WM_day_QC(isnan(USept_BMSbin.NEP_WM_day_QC))=[];
USept_BMSbin.NEC_WM_day_QC(isnan(USept_BMSbin.NEC_WM_day_QC))=[];


USept_inight = USept_BMSbin.PAR(1,:) < 1 %find all nightime datapoints 
% create nighttime hours datasets
% post-restoration night: 21:00 - 6:00am --> 10 hours

USept_BMSbin.SDN_night = USept_BMSbin.SDN;
USept_BMSbin.PAR_night = USept_BMSbin.PAR(1,:);

USept_BMSbin.NEP_night = USept_BMSbin.NEP;
USept_BMSbin.NEC_night = USept_BMSbin.NEC;
USept_BMSbin.NEP_night_QC = USept_BMSbin.NEP_QC;
USept_BMSbin.NEC_night_QC = USept_BMSbin.NEC_QC;

USept_BMSbin.dDOXY_night_QC = USept_BMSbin.dDOXY_QC;      %DO Gradient
USept_BMSbin.dTA_night_QC = USept_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
USept_BMSbin.SDN_night=USept_BMSbin.SDN_night(USept_inight);
USept_BMSbin.PAR_night=USept_BMSbin.PAR_night(USept_inight);

USept_BMSbin.NEP_night_QC=USept_BMSbin.NEP_night_QC(USept_inight);
USept_BMSbin.NEC_night_QC=USept_BMSbin.NEC_night_QC(USept_inight);

USept_BMSbin.dDOXY_night_QC=USept_BMSbin.dDOXY_night_QC(USept_inight);      %DO Gradient
USept_BMSbin.dTA_night_QC=USept_BMSbin.dTA_night_QC(USept_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(USept_BMSbin.SDN, USept_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(USept_BMSbin.SDN_night, USept_BMSbin.PAR_night, 'r.'); % day plot

datestr(USept_BMSbin.SDN_night)




%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(USept_BMSbin.NEP_day,USept_BMSbin.NEC_day);
% USept_BMSbin.Reg_Line = m*USept_BMSbin.NEP_day + b;
% USept_BMSbin.Ratio = m;
% USept_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(USept_BMSbin.NEP_day,USept_BMSbin.NEC_day,'o')
% plot(USept_BMSbin.NEP_day,USept_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Sept 2020 Pre-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + USept_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + USept_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(USept_BMSbin.NEP_day_QC,USept_BMSbin.NEC_day_QC);
USept_BMSbin.Reg_Line_QC = m_QC*USept_BMSbin.NEP_day_QC + b_QC;
USept_BMSbin.Ratio_QC = m_QC;
USept_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(USept_BMSbin.NEP_day_QC,USept_BMSbin.NEC_day_QC,'o')
plot(USept_BMSbin.NEP_day_QC,USept_BMSbin.Reg_Line_QC,'r')
%ylim([-50 50])
%xlim([-25 25])
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Sept 2020 Pre-Restoration NCC:NCP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + USept_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + USept_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(USept_BMSbin.NEP_WM_day,USept_BMSbin.NEC_WM_day);
% USept_BMSbin.Reg_Line_WM = m_WM*USept_BMSbin.NEP_WM_day + b_WM;
% USept_BMSbin.Ratio_WM = m_WM;
% USept_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(USept_BMSbin.NEP_WM_day,USept_BMSbin.NEC_WM_day,'o')
% plot(USept_BMSbin.NEP_WM_day,USept_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Sept 2020 Pre-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + USept_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + USept_BMSbin.R2_WM)


[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(USept_BMSbin.NEP_WM_day_QC,USept_BMSbin.NEC_WM_day_QC);
USept_BMSbin.Reg_Line_WM_QC = m_WM_QC*USept_BMSbin.NEP_WM_day_QC + b_WM_QC;
USept_BMSbin.Ratio_WM_QC = m_WM_QC;
USept_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(USept_BMSbin.NEP_WM_day_QC,USept_BMSbin.NEC_WM_day_QC,'o')
% plot(USept_BMSbin.NEP_WM_day_QC,USept_BMSbin.Reg_Line_WM_QC,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Sept 2020 Pre-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + USept_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + USept_BMSbin.R2_WM_QC)


%% For NEC:NEP Regressions Using Gradients
% close all 
clc

% multiply o2 gradient by -1 for O2 production
USept_BMSbin.dDOXY_Reg = -1.*USept_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
USept_BMSbin.dTA_Reg = 0.5.*USept_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(USept_BMSbin.dDOXY_day_QC);
figure 
hold on 
DOday = plot(1:Xlength, USept_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, USept_BMSbin.dDOXY_Reg);
xlabel('Sept Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, USept_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, USept_BMSbin.dTA_Reg);
xlabel('Sept Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Cudjoe Sept 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(USept_BMSbin.dDOXY_Reg, USept_BMSbin.dTA_Reg);
USept_BMSbin.Reg_Line_G = m_G*USept_BMSbin.dDOXY_Reg + b_G;
USept_BMSbin.Ratio_G = m_G;
USept_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(USept_BMSbin.dDOXY_Reg,USept_BMSbin.dTA_Reg,'o')
plot(USept_BMSbin.dDOXY_Reg,USept_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Sept 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + USept_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + USept_BMSbin.R2_G)

clc
disp('Finished with Sept Metabolism Calculations');

save('USept20_2.mat', 'USept_BMS', 'USept_BMSbin');

%% Plot Profiles 
% plot for profile within pump heights
% close all 
% for i =1:100 %length(USept_ADavg.SDN)
%     figure (i)
%     scatter(USept_ADavg.uv(1:108,i), USept_ADavg.bin_depth(1:108))
%     title(['Cudjoe Sept Velocity Profile Number ',num2str(i),])
%     xlabel('Velocity (m/s)');
%     ylabel('Height (m)');
% end


%% Subplots 
close all
clc


sgtitle('Cudjoe September 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(USept_BMSbin.SDN, USept_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(USept_BMSbin.SDN, USept_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
% DOplot = plot(USept_BMSbin.SDN, USept_BMSbin.dDOXY, 'c'); 
% TAplot = plot(USept_BMSbin.SDN, USept_BMSbin.dTA, 'k-'); 
plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'09/03 12:00';'09/04 12:00';'09/05 12:00';...
    '09/06 12:00';'09/07 12:00';'09/08 12:00';...
    '09/09 12:00';'09/10 12:00';'09/11 12:00';'09/12 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
% legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(USept_BMSbin.SDN, USept_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(USept_BMSbin.SDN, zeros(size(USept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', USept_good_Xrange, 'XTick', USept_tick, 'xticklabel', USept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
xticklabels({'09/03 12:00';'09/04 12:00';'09/05 12:00';...
    '09/06 12:00';'09/07 12:00';'09/08 12:00';...
    '09/09 12:00';'09/10 12:00';'09/11 12:00';'09/12 12:00'})
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(USept_BMSbin.SDN, USept_BMSbin.NEC_QC, 24);
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
plot(USept_BMSbin.NEP_day_QC,USept_BMSbin.NEC_day_QC,'o')
plot(USept_BMSbin.NEP_day_QC,USept_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(USept_BMSbin.Ratio_QC,2);
str2 = num2str(USept_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(USept_BMSbin.dDOXY_Reg,USept_BMSbin.dTA_Reg,'o')
plot(USept_BMSbin.dDOXY_Reg,USept_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
title('NEC:NEP Ratio from Gradients');
str1 = num2str(USept_BMSbin.Ratio_G,2);
str2 = num2str(USept_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')

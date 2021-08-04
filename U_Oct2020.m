% Oct SeapHOx Data Analysis from Cudjoe Ledge Reef
% Michelle Platz - USF 
% 3/9/2021

% SeapHOx sensor deployed 9/2/2020 
    % SP datafile: 'U_902BMS.txt'
    % pump 1 height above benthos = 70 cm  
    % pump 2 height above benthos = 20 cm 
% ADP sensor deployed 9/02/2020 
    % ADCP datafile: 'U_9_02'
    % height from substrate to ADCP head = 18 cm
clear all
close all
clc

%% Initial look at data
% ***** create UOct_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
UOct_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('U_928BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
UOct_SPraw.DOXY = UOct_SPraw.O2SATPER.*calcO2sat(UOct_SPraw.MCAT_TC, UOct_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
UOct_SPraw.pHint_prelim = calc_dfet_pHint(UOct_SPraw.Vint, UOct_SPraw.DFET_TC, -0.4);

% ***** create UOct_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
UOct_SP = parse_to_pumpheights_ARM_2pump_Mar19(UOct_SPraw, [datenum('09-28-2020 12:23:00'), datenum('10-18-2020 19:22:30')]);

% Calculate Gradients 
UOct_SP = calc_TA_gradientV2(UOct_SP, 2369.19, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

UOct_SP.dDOXY = UOct_SP.DOXY(1,:) - UOct_SP.DOXY(2,:); %Oxygen Gradient
UOct_SP.dpH = UOct_SP.pH(1,:) - UOct_SP.pH(2,:); %pH Gradient 
UOct_SP.dTA = UOct_SP.TAtop - UOct_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
UOct_DateString = {'09/28/2020 12:00:00';'09/29/2020 12:00:00';'09/30/2020 12:00:00';...
    '10/01/2020 12:00:00';'10/02/2020 12:00:00';'10/03/2020 12:00:00';'10/04/2020 12:00:00';'10/05/2020 12:00:00';...
    '10/06/2020 12:00:00';'10/07/2020 12:00:00';'10/08/2020 12:00:00';'10/09/2020 12:00:00';'10/10/2020 12:00:00';...
    '10/11/2020 12:00:00';'10/12/2020 12:00:00';'10/13/2020 12:00:00';'10/14/2020 12:00:00';'10/15/2020 12:00:00';...
    '10/16/2020 12:00:00';'10/17/2020 12:00:00';'10/18/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
UOct_tick = datenum(UOct_DateString,formatIn);

UOct_Xrange = [datenum('09-28-2020 12:23:00'), datenum('10-18-2020 19:22:30')];

figure
hold on; box on;
plot(UOct_SP.SDN, UOct_SP.dDOXY); %oxygen gradient 
plot(UOct_SP.SDN, UOct_SP.dTA); %TA gradient
plot(UOct_SP.SDN, zeros(size(UOct_SP.SDN))); %zero line
set(gca, 'xlim', UOct_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe Oct 2020 Unbinned DO and TA Gradients');


%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
UOct_SiteChar.SDN = UOct_SPraw.SDN;
UOct_SiteChar.TC = UOct_SPraw.OPT_TC;
UOct_SiteChar.PAR = UOct_SPraw.PAR;
UOct_SiteChar.PSAL = UOct_SPraw.PSAL;
UOct_SiteChar.Pres = UOct_SPraw.Pres;

%Plot pressure data to see when surface interval observations are
close all
figure 
hold on; 
plot(UOct_SiteChar.SDN, UOct_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
UOct_SiteChar.SDN = UOct_SPraw.SDN(15:end);
UOct_SiteChar.TC = UOct_SPraw.OPT_TC(15:end);
UOct_SiteChar.PAR = UOct_SPraw.PAR(15:end);
UOct_SiteChar.PSAL = UOct_SPraw.PSAL(15:end);
UOct_SiteChar.Pres = UOct_SPraw.Pres(15:end);

%extract full length of ADCP datafile  --> already done in Oct 
% % % UOct_ADfull=aquadoppraw2mat('U_9_02', 70, [datenum('09-02-2020 16:00:00'), datenum('10-30-2020 05:57:30')]);
% % % 
% % % % add AD variables to Site Char 
% % % UOct_SiteChar.AD_SDN = UOct_ADfull.SDN;
% % % UOct_SiteChar.AD_Pres = UOct_ADfull.Pres;
% % % UOct_SiteChar.AD_TC = UOct_ADfull.TC;
% % % UOct_SiteChar.bin_depth = UOct_ADfull.bin_depth;
% % % UOct_SiteChar.u = UOct_ADfull.u;
% % % UOct_SiteChar.v = UOct_ADfull.v;
% % % UOct_SiteChar.w = UOct_ADfull.w;
% % % UOct_SiteChar.uv = UOct_ADfull.uv;
% % % UOct_SiteChar.direction = UOct_ADfull.direction;
% % % 
% % % %Plot pressure data to see when surface interval observations are
% % % close all
% % % figure 
% % % hold on; 
% % % plot(UOct_SiteChar.AD_SDN, UOct_SiteChar.AD_Pres)
% % % 
% % % %clip ends of data to remove surfave interval observations 
% % % UOct_SiteChar.AD_SDN = UOct_ADfull.SDN(71:end);
% % % UOct_SiteChar.AD_Pres = UOct_ADfull.Pres(71:end);
% % % UOct_SiteChar.AD_TC = UOct_ADfull.TC(71:end);
% % % UOct_SiteChar.bin_depth = UOct_ADfull.bin_depth;
% % % UOct_SiteChar.u = UOct_ADfull.u(:,71:end);
% % % UOct_SiteChar.v = UOct_ADfull.v(:,71:end);
% % % UOct_SiteChar.w = UOct_ADfull.w(:,71:end);
% % % UOct_SiteChar.uv = UOct_ADfull.uv(:,71:end);
% % % UOct_SiteChar.direction = UOct_ADfull.direction(:,71:end);
% % % 
% % % % find U0 
% % % UOct_z1 = 0.70;
% % % UOct_z2 = 0.20;
% % % UOct_ADheight = 0.18;
% % % UOct_ADbin_depth_1m = 1-(UOct_ADheight);% = 0.82
% % % UOct_i1m = find(UOct_SiteChar.bin_depth==(0.82));
% % % UOct_SiteChar.U0 = UOct_SiteChar.uv(UOct_i1m,:);

% save data in separate datastructure
save('UOct20_SiteChar_2.mat', 'UOct_SiteChar' )

%% Constrain Xrange from graph results and extract good gradient data - 
close all 

UOct_good_Xrange = [datenum('28-Sep-2020 13:23:00'), datenum('13-Oct-2020 08:08:00')];

% plot to check range is correct
figure
hold on; box on;
plot(UOct_SP.SDN, UOct_SP.dDOXY); %oxygen gradient 
plot(UOct_SP.SDN, UOct_SP.dTA); %TA gradient
plot(UOct_SP.SDN, zeros(size(UOct_SP.SDN))); %zero line
set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe Oct 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all

% two segments of data 

U_Oct_BMS_idx_start = find(UOct_SP.SDN==datenum('28-Sep-2020 13:38:00'))
% U_Oct_BMS_idx_end1 = find(UOct_SP.SDN==datenum('10-01-2020 00:23:00'))
% U_Oct_BMS_idx_start2 = find(UOct_SP.SDN==datenum('10-03-2020 12:23:00'))
U_Oct_BMS_idx_end = find(UOct_SP.SDN==datenum('13-Oct-2020 08:23:00'))

% Create new data vectors of just the good data
U_Oct_BMS_good_data = U_Oct_BMS_idx_start:U_Oct_BMS_idx_end;
Initial_data_points = length(U_Oct_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(UOct_SP);
for v = 1:length(vars)
    UOct_SP.(vars{v}) = (UOct_SP.(vars{v})(:,U_Oct_BMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: UOct_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
UOct_AD=aquadoppraw2mat('U_9_02', 70, [datenum('28-Sep-2020 13:38:00'), datenum('13-Oct-2020 08:38:00')]);%***********add 15 mins to SP end time

%averages data to the middle of the minute interval spacified 
UOct_ADavg = average_aquadopp(UOct_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.7m (pump 1) and 0.2m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.52m (bin 42) and 0.02m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
UOct_Sal_est = mean(UOct_SP.PSAL(1,3:end));

[UOct_ADavg] = ustar_from_aquadopp2(UOct_ADavg,[0.52 0.11], UOct_Sal_est); %bins adjusted 

clc
%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[UOct_ADavg] = ustar_McGillis_Method(UOct_ADavg, 0.70, 0.20, 42, 1);

%compare
close all
figure
hold on 
ustar_plot = plot(UOct_ADavg.SDN, UOct_ADavg.ustar, 'r');
ustar_WM_plot = plot(UOct_ADavg.SDN, UOct_ADavg.ustar_WM);
plot(UOct_ADavg.SDN, zeros(size(UOct_ADavg.SDN)),'k');
set(gca, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UOct_ADavg.SDN(1) UOct_ADavg.SDN(end)])
xlabel('Oct Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: UOct_BMS *****
% BMS = AD + SP datasets (full datasets)

ADavg_vars = fieldnames(UOct_ADavg);
for v = 1:length(ADavg_vars)
    UOct_BMS.(ADavg_vars{v}) = (UOct_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(UOct_SP);
for v = 1:length(SP_vars)
    UOct_BMS.(SP_vars{v}) = (UOct_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(UOct_BMS.SDN)
% min: 13-28-43-58 becuase SP was restarted in the field rather than on the
% minute
%% %% *************** Calculate Fluxes ****************

% actual pump heights  = 0.70m (pump 1) and 0.20m (pump 2) 
% 0.18m from substrate to ACDP head in Oct at U 
% adjusted height = 0.52 m (bin 42) and 0.02 m (too shallow) (bin 1)  
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[UOct_BMS] = calc_NCC_3(UOct_BMS,[0.52 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(UOct_BMS.DOXY(1,:));
[UOct_BMS] = calc_NCP_3(UOct_BMS, [0.52 0.11],C1guess); %estimate C1 guess using median DOXY(1,:) value  

% Plot NCP and NCC
close all
figure
hold on; box on; 
NEPplot = plot(UOct_BMS.SDN, UOct_BMS.NEP);
NECplot = plot(UOct_BMS.SDN, UOct_BMS.NEC);
plot(UOct_BMS.SDN, zeros(size(UOct_BMS.SDN)),'k');
set(gca, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UOct_BMS.SDN(1) UOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Fluxes');

%McGillis method flux calculations 
[UOct_BMS] = calc_NCP_McGillis_Method(UOct_BMS, 0.70, 0.20, UOct_Sal_est);
[UOct_BMS] = calc_NCC_McGillis_Method(UOct_BMS, 0.70, 0.20, UOct_Sal_est);

% Plot NCP and NCC
% close all
figure
hold on; box on; 
NEPplot = plot(UOct_BMS.SDN, UOct_BMS.NEP_WM);
NECplot = plot(UOct_BMS.SDN, UOct_BMS.NEC_WM);
plot(UOct_BMS.SDN, zeros(size(UOct_BMS.SDN)),'k');
set(gca, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UOct_BMS.SDN(1) UOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Cudjoe Oct 2020 WM Fluxes');

%% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(UOct_BMS.SDN, UOct_BMS.NEP);
NEPplotWM = plot(UOct_BMS.SDN, UOct_BMS.NEP_WM);
plot(UOct_BMS.SDN, zeros(size(UOct_BMS.SDN)),'k');
set(gca, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UOct_BMS.SDN(1) UOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Cudjoe Oct 2020 NEP Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(UOct_BMS.SDN, UOct_BMS.NEC);
NECplotWM = plot(UOct_BMS.SDN, UOct_BMS.NEC_WM);
plot(UOct_BMS.SDN, zeros(size(UOct_BMS.SDN)),'k');
set(gca, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UOct_BMS.SDN(1) UOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Cudjoe Oct 2020 NEC Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
UOct_BMS.DOXYstd = std(UOct_BMS.DOXY);

% get DOXY std
UOct_idoxystd = find(UOct_BMS.DOXYstd > 2);
UOct_ihighdoxystd = [];
for i = 1:length(UOct_idoxystd)
    
    UOct_ihighdoxystd = vertcat(UOct_ihighdoxystd,[UOct_idoxystd(i)-1:1:UOct_idoxystd(i)+1]');
end
% get unique IDs
UOct_ihighdoxystd = unique(UOct_ihighdoxystd);
% remove 0's and out of index values
UOct_ihighdoxystd(UOct_ihighdoxystd==0) = [];
UOct_ihighdoxystd(UOct_ihighdoxystd> length(UOct_BMS.SDN)) = [];

% make it into index
trex = false(size(UOct_BMS.SDN));
trex(UOct_ihighdoxystd) = true;
UOct_ihighdoxystd = trex;
clear trex;

UOct_BMS.NEP_QC = UOct_BMS.NEP;
UOct_BMS.NEC_QC = UOct_BMS.NEC;
UOct_BMS.dDOXY_QC = UOct_BMS.dDOXY; %DO gradient
UOct_BMS.dTA_QC = UOct_BMS.dTA;     %TA gradient
UOct_BMS.NEP_WM_QC = UOct_BMS.NEP_WM;
UOct_BMS.NEC_WM_QC = UOct_BMS.NEC_WM;

% set observations when DOXYstd>0.8 to NaN
UOct_BMS.NEP_QC(UOct_ihighdoxystd) = NaN;  %NEP flux
UOct_BMS.NEC_QC(UOct_ihighdoxystd) = NaN;%NEC flux
UOct_BMS.dDOXY_QC(UOct_ihighdoxystd) = NaN; %DO gradient
UOct_BMS.dTA_QC(:,UOct_ihighdoxystd) = NaN;   %TA gradient
UOct_BMS.NEP_WM_QC(UOct_ihighdoxystd) = NaN;
UOct_BMS.NEC_WM_QC(:,UOct_ihighdoxystd) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(UOct_BMS.SDN, UOct_BMS.NEP, 'k');
NEPplotQC = plot(UOct_BMS.SDN, UOct_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(UOct_BMS.SDN, zeros(size(UOct_BMS.SDN)),'k');
set(gca, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UOct_BMS.SDN(1) UOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Fluxes');


%WM Plot
%close all
figure
hold on; box on;
NEPplot = plot(UOct_BMS.SDN, UOct_BMS.NEP_WM, 'k');
NEPplotQC = plot(UOct_BMS.SDN, UOct_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
plot(UOct_BMS.SDN, zeros(size(UOct_BMS.SDN)),'k');
set(gca, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UOct_BMS.SDN(1) UOct_BMS.SDN(end)])
xlabel('Oct Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Fluxes WM');


%% Bin data to hourly intervals

X = floor(nanmin(UOct_BMS.SDN)):1/24:ceil(nanmax(UOct_BMS.SDN));

UOct_BMSbin.SDN = X;
% variables from 2 differnet heights
UOct_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.DOXY(1,:), X);
UOct_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.DOXY(2,:), X);

UOct_BMSbin.pH(1,:)       = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.pH(1,:), X);
UOct_BMSbin.pH(2,:)       = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.pH(2,:), X);

UOct_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.PSAL(1,:), X);
UOct_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.PSAL(2,:), X);

UOct_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.O2SATPER(1,:), X);
UOct_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.O2SATPER(2,:), X);

UOct_BMSbin.Pres(1,:)     = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.Pres(1,:), X);
UOct_BMSbin.Pres(2,:)     = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.Pres(2,:), X);

UOct_BMSbin.DENS(1,:)     = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.DENS(1,:), X);
UOct_BMSbin.DENS(2,:)     = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.DENS(2,:), X);

UOct_BMSbin.PAR(1,:)      = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.PAR(1,:), X);
UOct_BMSbin.PAR(2,:)      = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.PAR(2,:), X);

UOct_BMSbin.bin_depth     = UOct_BMS.bin_depth;

for i = 1:108
    UOct_BMSbin.uv(i,:)   = bin_data_to_X_GF(UOct_BMS.SDN,UOct_BMS.uv(i,:), X);
end

% bin data hourly. Vector variables 
UOct_BMSbin.PRES  = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.Pres, X);
UOct_BMSbin.U0       = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.U0, X);
UOct_BMSbin.DIR      = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.direction, X);
UOct_BMSbin.ustar    = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.ustar, X);
UOct_BMSbin.ustar_rm = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.ustar_runmean, X);
UOct_BMSbin.dTA      = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.dTA, X);
UOct_BMSbin.dTA_QC   = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.dTA_QC, X);
UOct_BMSbin.dpH      = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.dpH, X);
UOct_BMSbin.dDOXY    = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.dDOXY, X);
UOct_BMSbin.dDOXY_QC = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.dDOXY_QC, X);
UOct_BMSbin.NEP      = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.NEP, X);
UOct_BMSbin.NEP_QC   = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.NEP_QC, X);
UOct_BMSbin.NEC      = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.NEC, X);
UOct_BMSbin.NEC_QC   = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.NEC_QC, X);
UOct_BMSbin.NEP_WM   = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.NEP_WM, X);
UOct_BMSbin.NEP_WM_QC= bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.NEP_WM_QC, X);
UOct_BMSbin.NEC_WM   = bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.NEC_WM, X);
UOct_BMSbin.NEC_WM_QC= bin_data_to_X_GF(UOct_BMS.SDN, UOct_BMS.NEC_WM_QC, X);


%% Plot Binned Fluxes 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP, 'b'); 
% NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC, 'r-'); 
% plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Oct 2020 Hourly Binned Fluxes FluxFit');

% close all
clc
figure
hold on; box on;
NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC_QC, 'r-'); 
plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Hourly Binned Fluxes FluxFit QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC_WM, 'r-'); 
% plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Oct 2020 WM Original Hourly Binned Fluxes');

% close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC_WM_QC, 'r-'); 
% plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Oct 2020 WM Full QC Hourly Binned Fluxes');
% 

%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
DOplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.dDOXY, 'c'); 
TAplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.dTA, 'k-'); 
plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Binned Gradiets');
%% ***QC*** Remove sections when velocity is too slow
ibad = UOct_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow
sum(ibad)
% UOct_BMSbin.NEP(ibad) = NaN;
% UOct_BMSbin.NEC(ibad) = NaN;
UOct_BMSbin.NEP_QC(ibad) = NaN;
UOct_BMSbin.NEC_QC(ibad) = NaN;
UOct_BMSbin.dDOXY_QC(ibad) = NaN;
UOct_BMSbin.dTA_QC(ibad) = NaN;

% UOct_BMSbin.NEP_WM(ibad) = NaN;
% UOct_BMSbin.NEC_WM(ibad) = NaN;
UOct_BMSbin.NEP_WM_QC(ibad) = NaN;
UOct_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP, 'b'); 
% NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC, 'r-'); 
% plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Oct 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC_QC, 'r-'); 
plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC_WM, 'r-'); 
% plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Oct 2020 Hourly Binned WM Fluxes');

% clc
% figure
% hold on; box on;
% NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC_WM_QC, 'r-'); 
% plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Oct Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Oct 2020 Hourly Binned WM Fluxes Full QC');

%% remove data anomaly between 10/1 and 10/3 and spotty ustar section between 8/4 and 8/8

U_Oct_BMS_idx_start2 = find(UOct_BMSbin.SDN==datenum('10-01-2020 00:00:00'));
U_Oct_BMS_idx_end2 = find(UOct_BMSbin.SDN==datenum('10-08-2020 01:00:00'));

ibad2 = U_Oct_BMS_idx_start2:U_Oct_BMS_idx_end2;
length(ibad2)

UOct_BMSbin.NEP(ibad2) = NaN;
UOct_BMSbin.NEC(ibad2) = NaN;
UOct_BMSbin.NEP_QC(ibad2) = NaN;
UOct_BMSbin.NEC_QC(ibad2) = NaN;

UOct_BMSbin.dDOXY_QC(ibad2) = NaN;
UOct_BMSbin.dTA_QC(ibad2) = NaN;

UOct_BMSbin.NEP_WM(ibad2) = NaN;
UOct_BMSbin.NEC_WM(ibad2) = NaN;
UOct_BMSbin.NEP_WM_QC(ibad2) = NaN;
UOct_BMSbin.NEC_WM_QC(ibad2) = NaN;


%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
figure
hold on; box on;
boxplot(UOct_BMSbin.NEP_QC)
ylabel('NEP')
title('Oct NEP Boxplots')

figure
hold on; box on;
boxplot(UOct_BMSbin.NEC_QC)
ylabel('NEC')
title('Oct NEC Boxplots')

figure
hold on; box on;
boxplot(UOct_BMSbin.dDOXY_QC)
ylabel('DO')
title('Oct dDO Boxplots')

figure
hold on; box on;
boxplot(UOct_BMSbin.dTA_QC)
ylabel('TA')
title('Oct dTA Boxplots')

%Outliers: 

% 355: NEP outlier (value: -32)
    % '12-Oct-2020 18:00:00'
    
% 301: dDO outlier (value-1.7) and dTA outlier 
    % bad profile 
    % OUTLIER 

% 330: NEC and dTA outlier 
    % OUTLIER 

% 253: NEC outlier
    % not an outlier 
    
% 326: NEC outlier
    % bad profile - OUTLIER 
    
% 331: NEC outlier
    % '11-Oct-2020 18:00:00'
    
% 254: NEC outlier
    % bad profile - OUTLIER 

% 334: NEC outlier and : dTA outlier 
    % '11-Oct-2020 21:00:00'
    % bad profile - OUTLIER 

% 266: NEC outlier
    % bad profile - OUTLIER 

% 369: dTA outlier 
    % GOOD PROFILE 
    
% 325: dTA outlier 
    % bad profile - OUTLIER 

datestr(UOct_BMSbin.SDN(334))

% Plot Profiles at outliers - 
close all 
for i = 325   %1:length(UAug2_BMSbin.SDN)
    figure (i)
    scatter(UOct_BMSbin.uv(1:108,i), UOct_BMSbin.bin_depth(1:108))
    title(['Cudjoe Oct Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end
%outliers to be removed: 355, 301, 330, 326, 254, 334, 266, 325

% UOct_BMSbin.NEP_QC(355) = NaN;
% UOct_BMSbin.NEC_QC(355) = NaN;
% UOct_BMSbin.dDOXY_QC(355) = NaN;
% UOct_BMSbin.dTA_QC(355) = NaN;
% 
% UOct_BMSbin.NEP_QC(301) = NaN;
% UOct_BMSbin.NEC_QC(301) = NaN;
% UOct_BMSbin.dDOXY_QC(301) = NaN;
% UOct_BMSbin.dTA_QC(301) = NaN;
% 
% UOct_BMSbin.NEP_QC(330) = NaN;
% UOct_BMSbin.NEC_QC(330) = NaN;
% UOct_BMSbin.dDOXY_QC(330) = NaN;
% UOct_BMSbin.dTA_QC(330) = NaN;
% 
% UOct_BMSbin.NEP_QC(326) = NaN;
% UOct_BMSbin.NEC_QC(326) = NaN;
% UOct_BMSbin.dDOXY_QC(326) = NaN;
% UOct_BMSbin.dTA_QC(326) = NaN;
% 
% UOct_BMSbin.NEP_QC(254) = NaN;
% UOct_BMSbin.NEC_QC(254) = NaN;
% UOct_BMSbin.dDOXY_QC(254) = NaN;
% UOct_BMSbin.dTA_QC(254) = NaN;
% 
% UOct_BMSbin.NEP_QC(334) = NaN;
% UOct_BMSbin.NEC_QC(334) = NaN;
% UOct_BMSbin.dDOXY_QC(334) = NaN;
% UOct_BMSbin.dTA_QC(334) = NaN;
% 
% UOct_BMSbin.NEP_QC(266) = NaN;
% UOct_BMSbin.NEC_QC(266) = NaN;
% UOct_BMSbin.dDOXY_QC(266) = NaN;
% UOct_BMSbin.dTA_QC(266) = NaN;
% 
% UOct_BMSbin.NEP_QC(325) = NaN;
% UOct_BMSbin.NEC_QC(325) = NaN;
% UOct_BMSbin.dDOXY_QC(325) = NaN;
% UOct_BMSbin.dTA_QC(325) = NaN;


close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC_QC, 'r-'); 
plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Hourly Binned Fluxes Full QC');





%% Plot Diel Curves

% nepdbin = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEP, 24);
% necdbin = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEC, 24);
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

nepdbin_QC = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEC_QC, 24);
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
% nepdbin_WM = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEC_WM, 24);
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

nepdbin_WM_QC = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEC_WM_QC, 24);
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
% Extract daytime data using UOct_BMSbin.PAR
UOct_inight = UOct_BMSbin.PAR(1,:) < 5; %find all nightime datapoints 

%create new arrays for daytime data
UOct_BMSbin.SDN_day = UOct_BMSbin.SDN;
UOct_BMSbin.PAR_day = UOct_BMSbin.PAR(1,:);

UOct_BMSbin.NEP_day = UOct_BMSbin.NEP;
UOct_BMSbin.NEC_day = UOct_BMSbin.NEC;
UOct_BMSbin.NEP_day_QC = UOct_BMSbin.NEP_QC;
UOct_BMSbin.NEC_day_QC = UOct_BMSbin.NEC_QC;

UOct_BMSbin.dDOXY_day_QC = UOct_BMSbin.dDOXY_QC;      %DO Gradient
UOct_BMSbin.dTA_day_QC = UOct_BMSbin.dTA_QC;          %TA Gradient

UOct_BMSbin.NEP_WM_day = UOct_BMSbin.NEP_WM;
UOct_BMSbin.NEC_WM_day = UOct_BMSbin.NEC_WM;
UOct_BMSbin.NEP_WM_day_QC = UOct_BMSbin.NEP_WM_QC;
UOct_BMSbin.NEC_WM_day_QC = UOct_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
UOct_BMSbin.SDN_day(UOct_inight) = NaN;
UOct_BMSbin.PAR_day (UOct_inight) = NaN;

UOct_BMSbin.NEP_day(UOct_inight) = NaN;
UOct_BMSbin.NEC_day(UOct_inight) = NaN;
UOct_BMSbin.NEP_day_QC(UOct_inight) = NaN;
UOct_BMSbin.NEC_day_QC(UOct_inight) = NaN;

UOct_BMSbin.dDOXY_day_QC(UOct_inight) = NaN;      %DO Gradient
UOct_BMSbin.dTA_day_QC(UOct_inight) = NaN;        %TA Gradient

UOct_BMSbin.NEP_WM_day(UOct_inight) = NaN;
UOct_BMSbin.NEC_WM_day(UOct_inight) = NaN;
UOct_BMSbin.NEP_WM_day_QC(UOct_inight) = NaN;
UOct_BMSbin.NEC_WM_day_QC(UOct_inight) = NaN;

%Plot to check only nighttime points removed
figure 
hold on
scatter(UOct_BMSbin.SDN, UOct_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(UOct_BMSbin.SDN_day, UOct_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
UOct_BMSbin.NEP_day(isnan(UOct_BMSbin.NEP_day))=[];
UOct_BMSbin.NEC_day(isnan(UOct_BMSbin.NEC_day))=[];
UOct_BMSbin.NEP_day_QC(isnan(UOct_BMSbin.NEP_day_QC))=[];
UOct_BMSbin.NEC_day_QC(isnan(UOct_BMSbin.NEC_day_QC))=[];

UOct_BMSbin.dDOXY_day_QC(isnan(UOct_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
UOct_BMSbin.dTA_day_QC(isnan(UOct_BMSbin.dTA_day_QC))=[];       %TA Gradient

UOct_BMSbin.NEP_WM_day(isnan(UOct_BMSbin.NEP_WM_day))=[];
UOct_BMSbin.NEC_WM_day(isnan(UOct_BMSbin.NEC_WM_day))=[];
UOct_BMSbin.NEP_WM_day_QC(isnan(UOct_BMSbin.NEP_WM_day_QC))=[];
UOct_BMSbin.NEC_WM_day_QC(isnan(UOct_BMSbin.NEC_WM_day_QC))=[];

UOct_inight = UOct_BMSbin.PAR(1,:) < 1; %find all nightime datapoints 
% post-restoration night: 20:00 - 6:00am --> 11 hours

% create nighttime hours datasets
UOct_BMSbin.SDN_night = UOct_BMSbin.SDN;
UOct_BMSbin.PAR_night = UOct_BMSbin.PAR(1,:);

UOct_BMSbin.NEP_night = UOct_BMSbin.NEP;
UOct_BMSbin.NEC_night = UOct_BMSbin.NEC;
UOct_BMSbin.NEP_night_QC = UOct_BMSbin.NEP_QC;
UOct_BMSbin.NEC_night_QC = UOct_BMSbin.NEC_QC;

UOct_BMSbin.dDOXY_night_QC = UOct_BMSbin.dDOXY_QC;      %DO Gradient
UOct_BMSbin.dTA_night_QC = UOct_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
UOct_BMSbin.SDN_night=UOct_BMSbin.SDN_night(UOct_inight);
UOct_BMSbin.PAR_night=UOct_BMSbin.PAR_night(UOct_inight);

UOct_BMSbin.NEP_night_QC=UOct_BMSbin.NEP_night_QC(UOct_inight);
UOct_BMSbin.NEC_night_QC=UOct_BMSbin.NEC_night_QC(UOct_inight);

UOct_BMSbin.dDOXY_night_QC=UOct_BMSbin.dDOXY_night_QC(UOct_inight);      %DO Gradient
UOct_BMSbin.dTA_night_QC=UOct_BMSbin.dTA_night_QC(UOct_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(UOct_BMSbin.SDN, UOct_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(UOct_BMSbin.SDN_night, UOct_BMSbin.PAR_night, 'r.'); % day plot

%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(UOct_BMSbin.NEP_day,UOct_BMSbin.NEC_day);
% UOct_BMSbin.Reg_Line = m*UOct_BMSbin.NEP_day + b;
% UOct_BMSbin.Ratio = m;
% UOct_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(UOct_BMSbin.NEP_day,UOct_BMSbin.NEC_day,'o')
% plot(UOct_BMSbin.NEP_day,UOct_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Oct 2020 Pre-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UOct_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UOct_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(UOct_BMSbin.NEP_day_QC,UOct_BMSbin.NEC_day_QC);
UOct_BMSbin.Reg_Line_QC = m_QC*UOct_BMSbin.NEP_day_QC + b_QC;
UOct_BMSbin.Ratio_QC = m_QC;
UOct_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(UOct_BMSbin.NEP_day_QC,UOct_BMSbin.NEC_day_QC,'o')
plot(UOct_BMSbin.NEP_day_QC,UOct_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('Cudjoe Oct 2020 Post-Restoration NEC:NEP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + UOct_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UOct_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(UOct_BMSbin.NEP_WM_day,UOct_BMSbin.NEC_WM_day);
% UOct_BMSbin.Reg_Line_WM = m_WM*UOct_BMSbin.NEP_WM_day + b_WM;
% UOct_BMSbin.Ratio_WM = m_WM;
% UOct_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(UOct_BMSbin.NEP_WM_day,UOct_BMSbin.NEC_WM_day,'o')
% plot(UOct_BMSbin.NEP_WM_day,UOct_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Oct 2020 Pre-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UOct_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UOct_BMSbin.R2_WM)


[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(UOct_BMSbin.NEP_WM_day_QC,UOct_BMSbin.NEC_WM_day_QC);
UOct_BMSbin.Reg_Line_WM_QC = m_WM_QC*UOct_BMSbin.NEP_WM_day_QC + b_WM_QC;
UOct_BMSbin.Ratio_WM_QC = m_WM_QC;
UOct_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(UOct_BMSbin.NEP_WM_day_QC,UOct_BMSbin.NEC_WM_day_QC,'o')
% plot(UOct_BMSbin.NEP_WM_day_QC,UOct_BMSbin.Reg_Line_WM_QC,'r')
% xlabel('NEP');
% ylabel('NEC');
% title('Cudjoe Oct 2020 Post-Restoration NEC:NEP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NEC:NEP =" + UOct_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UOct_BMSbin.R2_WM_QC)
% 


%% For NEC:NEP Regressions Using Gradients
close all 
clc

% multiply o2 gradient by -1 for O2 production
UOct_BMSbin.dDOXY_Reg = -1.*UOct_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
UOct_BMSbin.dTA_Reg = 0.5.*UOct_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(UOct_BMSbin.dDOXY_day_QC);
figure 
hold on 
DOday = plot(1:Xlength, UOct_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, UOct_BMSbin.dDOXY_Reg);
xlabel('Oct Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, UOct_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, UOct_BMSbin.dTA_Reg);
xlabel('Oct Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Cudjoe Oct 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(UOct_BMSbin.dDOXY_Reg, UOct_BMSbin.dTA_Reg);
UOct_BMSbin.Reg_Line_G = m_G*UOct_BMSbin.dDOXY_Reg + b_G;
UOct_BMSbin.Ratio_G = m_G;
UOct_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(UOct_BMSbin.dDOXY_Reg,UOct_BMSbin.dTA_Reg,'o')
plot(UOct_BMSbin.dDOXY_Reg,UOct_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Oct 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UOct_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UOct_BMSbin.R2_G)


clc
disp('Finished with Oct Metabolism Calculations');
save('UOct20_2.mat', 'UOct_BMS', 'UOct_BMSbin');

%% Plot Profiles 
% plot for profile within pump heights
% close all 
% for i =1:100 %length(UOct_ADavg.SDN)
%     figure (i)
%     scatter(UOct_ADavg.uv(1:108,i), UOct_ADavg.bin_depth(1:108))
%     title(['Cudjoe Oct Velocity Profile Number ',num2str(i),])
%     xlabel('Velocity (m/s)');
%     ylabel('Height (m)');
% end


%% Subplots 
close all
clc

sgtitle('Cudjoe October 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'09/28 12:00';'09/29 12:00';'09/30 12:00';...
    '10/01 12:00';'10/02 12:00';'10/03 12:00';'10/04 12:00';'10/05 12:00';...
    '10/06 12:00';'10/07 12:00';'10/08 12:00';'10/09 12:00';'10/10 12:00';...
    '10/11 12:00';'10/12 12:00';'10/13 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
% legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(UOct_BMSbin.SDN, UOct_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(UOct_BMSbin.SDN, zeros(size(UOct_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UOct_good_Xrange, 'XTick', UOct_tick, 'xticklabel', UOct_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
xticklabels({'09/28 12:00';'09/29 12:00';'09/30 12:00';...
    '10/01 12:00';'10/02 12:00';'10/03 12:00';'10/04 12:00';'10/05 12:00';...
    '10/06 12:00';'10/07 12:00';'10/08 12:00';'10/09 12:00';'10/10 12:00';...
    '10/11 12:00';'10/12 12:00';'10/13 12:00'})
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(UOct_BMSbin.SDN, UOct_BMSbin.NEC_QC, 24);
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
plot(UOct_BMSbin.NEP_day_QC,UOct_BMSbin.NEC_day_QC,'o')
plot(UOct_BMSbin.NEP_day_QC,UOct_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
% ylim([-35 35]);
% xlim([-35 35]);
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(UOct_BMSbin.Ratio_QC,2);
str2 = num2str(UOct_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(UOct_BMSbin.dDOXY_Reg,UOct_BMSbin.dTA_Reg,'o')
plot(UOct_BMSbin.dDOXY_Reg,UOct_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
% ylim([-3 3]);
% xlim([-3 3]);
title('NEC:NEP Ratio from Gradients');
str1 = num2str(UOct_BMSbin.Ratio_G,2);
str2 = num2str(UOct_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')



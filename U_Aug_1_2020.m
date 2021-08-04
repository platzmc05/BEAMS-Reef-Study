% Aug_1 SeapHOx Data Analysis from Cudjoe Ledge Reef
% Michelle Platz - USF 
% 3/5/2021

% SeapHOx sensor deployed 8/05/2020 at 9:15 am EST
    % SP datafile: 'U_805BMS.txt'
    % pump 1 height above benthos = 70 cm  
    % pump 2 height above benthos = 20 cm 
% ADP sensor deployed 8/05/2020 at 9:15 am EST
    % ADCP datafile: 'U_8_05'
    % height from substrate to ADCP head = 18 cm

close all
clc
clear all

%% Initial look at data
% ***** create UAug1_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
UAug1_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('U_805BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
UAug1_SPraw.DOXY = UAug1_SPraw.O2SATPER.*calcO2sat(UAug1_SPraw.MCAT_TC, UAug1_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
UAug1_SPraw.pHint_prelim = calc_dfet_pHint(UAug1_SPraw.Vint, UAug1_SPraw.DFET_TC, -0.4);

% ***** create UAug1_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
UAug1_SP = parse_to_pumpheights_ARM_2pump_Mar19(UAug1_SPraw, [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 11:59:30')]);

% Calculate Gradients 
UAug1_SP = calc_TA_gradientV2(UAug1_SP, 2369.19, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

UAug1_SP.dDOXY = UAug1_SP.DOXY(1,:) - UAug1_SP.DOXY(2,:); %Oxygen Gradient
UAug1_SP.dpH = UAug1_SP.pH(1,:) - UAug1_SP.pH(2,:); %pH Gradient 
UAug1_SP.dTA = UAug1_SP.TAtop - UAug1_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
UAug1_DateString = {'08/05/2020 12:00:00';'08/06/2020 12:00:00';'08/07/2020 12:00:00';'08/08/2020 12:00:00';'08/09/2020 12:00:00';...
    '08/10/2020 12:00:00';'08/11/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
UAug1_tick = datenum(UAug1_DateString,formatIn);

UAug1_Xrange = [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 11:59:30')];

figure
hold on; box on;
plot(UAug1_SP.SDN, UAug1_SP.dDOXY); %oxygen gradient 
plot(UAug1_SP.SDN, UAug1_SP.dTA); %TA gradient
plot(UAug1_SP.SDN, zeros(size(UAug1_SP.SDN))); %zero line
set(gca, 'xlim', UAug1_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe Aug 2020 Unbinned DO and TA Gradients');

%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
UAug1_SiteChar.SDN = UAug1_SPraw.SDN;
UAug1_SiteChar.TC = UAug1_SPraw.OPT_TC;
UAug1_SiteChar.PAR = UAug1_SPraw.PAR;
UAug1_SiteChar.PSAL = UAug1_SPraw.PSAL;
UAug1_SiteChar.Pres = UAug1_SPraw.Pres;

close all
figure 
hold on; 
plot(UAug1_SiteChar.SDN, UAug1_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
UAug1_SiteChar.SDN = UAug1_SPraw.SDN(90:20641);
UAug1_SiteChar.TC = UAug1_SPraw.OPT_TC(90:20641);
UAug1_SiteChar.PAR = UAug1_SPraw.PAR(90:20641);
UAug1_SiteChar.PSAL = UAug1_SPraw.PSAL(90:20641);
UAug1_SiteChar.Pres = UAug1_SPraw.Pres(90:20641);

%extract full length of ADCP datafile  
UAug1_ADfull=aquadoppraw2mat('U_8_05', 70, [datenum('08-05-2020 08:00:00'), datenum('09-01-2020 08:00:00')]);

% add AD variables to Site Char 
UAug1_SiteChar.AD_SDN = UAug1_ADfull.SDN;
UAug1_SiteChar.AD_Pres = UAug1_ADfull.Pres;
UAug1_SiteChar.AD_TC = UAug1_ADfull.TC;
UAug1_SiteChar.bin_depth = UAug1_ADfull.bin_depth;
UAug1_SiteChar.u = UAug1_ADfull.u;
UAug1_SiteChar.v = UAug1_ADfull.v;
UAug1_SiteChar.w = UAug1_ADfull.w;
UAug1_SiteChar.uv = UAug1_ADfull.uv;
UAug1_SiteChar.direction = UAug1_ADfull.direction;

%Plot to see when surface interval observations are
close all
figure 
hold on; 
plot(UAug1_SiteChar.AD_SDN, UAug1_SiteChar.AD_Pres)

%clip ends of data to remove surfave interval observations 
UAug1_SiteChar.AD_SDN = UAug1_ADfull.SDN(95:end);
UAug1_SiteChar.AD_Pres = UAug1_ADfull.Pres(95:end);
UAug1_SiteChar.AD_TC = UAug1_ADfull.TC(95:end);
UAug1_SiteChar.bin_depth = UAug1_ADfull.bin_depth;
UAug1_SiteChar.u = UAug1_ADfull.u(:,95:end);
UAug1_SiteChar.v = UAug1_ADfull.v(:,95:end);
UAug1_SiteChar.w = UAug1_ADfull.w(:,95:end);
UAug1_SiteChar.uv = UAug1_ADfull.uv(:,95:end);
UAug1_SiteChar.direction = UAug1_ADfull.direction(:,95:end);

% find U0 
UAug1_z1 = 0.70;
UAug1_z2 = 0.20;
UAug1_ADheight = 0.18;
UAug1_ADbin_depth_1m = 1-(UAug1_ADheight);% = 0.82
UAug1_i1m = find(UAug1_SiteChar.bin_depth==(0.82));
UAug1_SiteChar.U0 = UAug1_SiteChar.uv(UAug1_i1m,:);

% save data in separate datastructure
save('UAug120_SiteChar_2.mat', 'UAug1_SiteChar' )


%% Constrain Xrange from graph results and extract good gradient data - 
close all 

UAug1_good_Xrange = [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 10:00:00')];

% plot to check range is correct
figure
hold on; box on;
plot(UAug1_SP.SDN, UAug1_SP.dDOXY); %oxygen gradient 
plot(UAug1_SP.SDN, UAug1_SP.dTA); %TA gradient
plot(UAug1_SP.SDN, zeros(size(UAug1_SP.SDN))); %zero line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe Aug 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all

U_Aug_1_BMS_idx_start = find(UAug1_SP.SDN==datenum('08-05-2020 10:00:00'))
U_Aug_1_BMS_idx_end = find(UAug1_SP.SDN==datenum('08-12-2020 10:00:00'))

% Create new data vectors of just the good data
U_Aug_1_BMS_good_data = U_Aug_1_BMS_idx_start:U_Aug_1_BMS_idx_end;
Initial_data_points = length(U_Aug_1_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(UAug1_SP);
for v = 1:length(vars)
    UAug1_SP.(vars{v}) = (UAug1_SP.(vars{v})(:,U_Aug_1_BMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: UAug1_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
UAug1_AD=aquadoppraw2mat('U_8_05', 70, [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 10:15:00')]);

%averages data to the middle of the minute interval spacified 
UAug1_ADavg = average_aquadopp(UAug1_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.7m (pump 1) and 0.2m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.52m (bin 42) and 0.02m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
UAug1_Sal_est = mean(UAug1_SP.PSAL(1,3:end));

[UAug1_ADavg] = ustar_from_aquadopp2(UAug1_ADavg,[0.52 0.11], UAug1_Sal_est); %bins adjusted 

clc
%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[UAug1_ADavg] = ustar_McGillis_Method(UAug1_ADavg, 0.70, 0.20, 42, 1);

%compare
close all
figure
hold on 
ustar_plot = plot(UAug1_ADavg.SDN, UAug1_ADavg.ustar, 'r');
ustar_WM_plot = plot(UAug1_ADavg.SDN, UAug1_ADavg.ustar_WM);
plot(UAug1_ADavg.SDN, zeros(size(UAug1_ADavg.SDN)),'k');
set(gca, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug1_ADavg.SDN(1) UAug1_ADavg.SDN(end)])
xlabel('Aug Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: UAug1_BMS *****

ADavg_vars = fieldnames(UAug1_ADavg);
for v = 1:length(ADavg_vars)
    UAug1_BMS.(ADavg_vars{v}) = (UAug1_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(UAug1_SP);
for v = 1:length(SP_vars)
    UAug1_BMS.(SP_vars{v}) = (UAug1_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(UAug1_BMS.SDN)


%% %% *************** Calculate Fluxes ****************

% actual pump heights  = 0.70m (pump 1) and 0.20m (pump 2) 
% 0.18m from substrate to ACDP head in Aug at U 
% adjusted height = 0.52 m (bin 42) and 0.02 m (too shallow) (bin 1)  
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[UAug1_BMS] = calc_NCC_3(UAug1_BMS,[0.52 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(UAug1_BMS.DOXY(1,2:end))
[UAug1_BMS] = calc_NCP_3(UAug1_BMS, [0.52 0.11],C1guess); %C1 guess - DOXY(1,x) guess is current saved, 

% Plot NCP and NCC
%close all
figure
hold on; box on; 
NEPplot = plot(UAug1_BMS.SDN, UAug1_BMS.NEP);
NECplot = plot(UAug1_BMS.SDN, UAug1_BMS.NEC);
plot(UAug1_BMS.SDN, zeros(size(UAug1_BMS.SDN)),'k');
set(gca, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug1_BMS.SDN(1) UAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');

%McGillis method flux calculations 
[UAug1_BMS] = calc_NCP_McGillis_Method(UAug1_BMS, 0.70, 0.20, UAug1_Sal_est);
[UAug1_BMS] = calc_NCC_McGillis_Method(UAug1_BMS, 0.70, 0.20, UAug1_Sal_est);

% Plot NCP and NCC
% close all
% figure
% hold on; box on; 
% NEPplot = plot(UAug1_BMS.SDN, UAug1_BMS.NEP_WM);
% NECplot = plot(UAug1_BMS.SDN, UAug1_BMS.NEC_WM);
% plot(UAug1_BMS.SDN, zeros(size(UAug1_BMS.SDN)),'k');
% set(gca, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlim([UAug1_BMS.SDN(1) UAug1_BMS.SDN(end)])
% xlabel('Aug Days');
% ylabel('NCP or NCC [mmol/m2/hr]');
% legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 WM Fluxes');

% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(UAug1_BMS.SDN, UAug1_BMS.NEP);
NEPplotWM = plot(UAug1_BMS.SDN, UAug1_BMS.NEP_WM);
plot(UAug1_BMS.SDN, zeros(size(UAug1_BMS.SDN)),'k');
set(gca, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug1_BMS.SDN(1) UAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(UAug1_BMS.SDN, UAug1_BMS.NEC);
NECplotWM = plot(UAug1_BMS.SDN, UAug1_BMS.NEC_WM);
plot(UAug1_BMS.SDN, zeros(size(UAug1_BMS.SDN)),'k');
set(gca, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug1_BMS.SDN(1) UAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
UAug1_BMS.DOXYstd = std(UAug1_BMS.DOXY);

% get DOXY std
UAug1_idoxystd = find(UAug1_BMS.DOXYstd > 2);
UAug1_ihighdoxystd = [];
for i = 1:length(UAug1_idoxystd)
    
    UAug1_ihighdoxystd = vertcat(UAug1_ihighdoxystd,[UAug1_idoxystd(i)-1:1:UAug1_idoxystd(i)+1]');
end
% get unique IDs
UAug1_ihighdoxystd = unique(UAug1_ihighdoxystd);
% remove 0's and out of index values
UAug1_ihighdoxystd(UAug1_ihighdoxystd==0) = [];
UAug1_ihighdoxystd(UAug1_ihighdoxystd> length(UAug1_BMS.SDN)) = [];

% make it into index
trex = false(size(UAug1_BMS.SDN));
trex(UAug1_ihighdoxystd) = true;
UAug1_ihighdoxystd = trex;
clear trex;

UAug1_BMS.NEP_QC = UAug1_BMS.NEP;
UAug1_BMS.NEC_QC = UAug1_BMS.NEC;
UAug1_BMS.dDOXY_QC = UAug1_BMS.dDOXY; %DO gradient
UAug1_BMS.dTA_QC = UAug1_BMS.dTA;     %TA gradient
UAug1_BMS.NEP_WM_QC = UAug1_BMS.NEP_WM;
UAug1_BMS.NEC_WM_QC = UAug1_BMS.NEC_WM;

% set observations when DOXYstd>0.8 to NaN
UAug1_BMS.NEP_QC(UAug1_ihighdoxystd) = NaN;
UAug1_BMS.NEC_QC(:,UAug1_ihighdoxystd) = NaN;
UAug1_BMS.dDOXY_QC(UAug1_ihighdoxystd) = NaN; %DO gradient
UAug1_BMS.dTA_QC(:,UAug1_ihighdoxystd) = NaN;   %TA gradient
UAug1_BMS.NEP_WM_QC(UAug1_ihighdoxystd) = NaN;
UAug1_BMS.NEC_WM_QC(:,UAug1_ihighdoxystd) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(UAug1_BMS.SDN, UAug1_BMS.NEP, 'k');
NEPplotQC = plot(UAug1_BMS.SDN, UAug1_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(UAug1_BMS.SDN, zeros(size(UAug1_BMS.SDN)),'k');
set(gca, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug1_BMS.SDN(1) UAug1_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');


%WM Plot
%close all
% figure
% hold on; box on;
% NEPplot = plot(UAug1_BMS.SDN, UAug1_BMS.NEP_WM, 'k');
% NEPplotQC = plot(UAug1_BMS.SDN, UAug1_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
% plot(UAug1_BMS.SDN, zeros(size(UAug1_BMS.SDN)),'k');
% set(gca, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlim([UAug1_BMS.SDN(1) UAug1_BMS.SDN(end)])
% xlabel('Aug Days');
% ylabel('NCP [mmol/m2/hr]');
% legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Fluxes');


%% Bin data to hourly intervals

X = floor(nanmin(UAug1_BMS.SDN)):1/24:ceil(nanmax(UAug1_BMS.SDN));

UAug1_BMSbin.SDN = X;
% variables from 2 differnet heights
UAug1_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.DOXY(1,:), X);
UAug1_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.DOXY(2,:), X);

UAug1_BMSbin.pH(1,:)       = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.pH(1,:), X);
UAug1_BMSbin.pH(2,:)       = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.pH(2,:), X);

UAug1_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.PSAL(1,:), X);
UAug1_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.PSAL(2,:), X);

UAug1_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.O2SATPER(1,:), X);
UAug1_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.O2SATPER(2,:), X);

UAug1_BMSbin.Pres(1,:)     = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.Pres(1,:), X);
UAug1_BMSbin.Pres(2,:)     = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.Pres(2,:), X);

UAug1_BMSbin.DENS(1,:)     = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.DENS(1,:), X);
UAug1_BMSbin.DENS(2,:)     = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.DENS(2,:), X);

UAug1_BMSbin.PAR(1,:)      = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.PAR(1,:), X);
UAug1_BMSbin.PAR(2,:)      = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.PAR(2,:), X);

UAug1_BMSbin.bin_depth     = UAug1_BMS.bin_depth;

for i = 1:108
    UAug1_BMSbin.uv(i,:)   = bin_data_to_X_GF(UAug1_BMS.SDN,UAug1_BMS.uv(i,:), X);
end

% bin data hourly. Vector variables 
UAug1_BMSbin.PRES  = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.Pres, X);
UAug1_BMSbin.U0       = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.U0, X);
UAug1_BMSbin.DIR      = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.direction, X);
UAug1_BMSbin.ustar    = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.ustar, X);
UAug1_BMSbin.ustar_rm = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.ustar_runmean, X);
UAug1_BMSbin.dTA      = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.dTA, X);
UAug1_BMSbin.dTA_QC   = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.dTA_QC, X);
UAug1_BMSbin.dpH      = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.dpH, X);
UAug1_BMSbin.dDOXY    = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.dDOXY, X);
UAug1_BMSbin.dDOXY_QC = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.dDOXY_QC, X);
UAug1_BMSbin.NEP      = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.NEP, X);
UAug1_BMSbin.NEP_QC   = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.NEP_QC, X);
UAug1_BMSbin.NEC      = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.NEC, X);
UAug1_BMSbin.NEC_QC   = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.NEC_QC, X);
UAug1_BMSbin.NEP_WM   = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.NEP_WM, X);
UAug1_BMSbin.NEP_WM_QC= bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.NEP_WM_QC, X);
UAug1_BMSbin.NEC_WM   = bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.NEC_WM, X);
UAug1_BMSbin.NEC_WM_QC= bin_data_to_X_GF(UAug1_BMS.SDN, UAug1_BMS.NEC_WM_QC, X);


%% Plot Binned Fluxes 
close all
clc
% figure
% hold on; box on;
% NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP, 'b'); 
% NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC, 'r-'); 
% plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Hourly Binned Fluxes');

% close all
clc
figure
hold on; box on;
NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_QC, 'r-'); 
plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Hourly Binned Fluxes');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_WM, 'r-'); 
% plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 WM Original Hourly Binned Fluxes');

% close all
clc
figure
hold on; box on;
NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_WM_QC, 'r-'); 
plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 WM Full QC Hourly Binned Fluxes');


%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
DOplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.dDOXY, 'c'); 
TAplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.dTA, 'k-'); 
plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug1 Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe Aug1 2020 Binned Gradiets');



%% ***QC*** Remove sections when velocity is too slow
ibad = UAug1_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow

% UAug1_BMSbin.NEP(ibad) = NaN;
% UAug1_BMSbin.NEC(ibad) = NaN;
UAug1_BMSbin.NEP_QC(ibad) = NaN;
UAug1_BMSbin.NEC_QC(ibad) = NaN;

UAug1_BMSbin.dDOXY_QC(ibad) = NaN;
UAug1_BMSbin.dTA_QC(ibad) = NaN;

% UAug1_BMSbin.NEP_WM(ibad) = NaN;
% UAug1_BMSbin.NEC_WM(ibad) = NaN;
UAug1_BMSbin.NEP_WM_QC(ibad) = NaN;
UAug1_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
clc
% figure
% hold on; box on;
% NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP, 'b'); 
% NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC, 'r-'); 
% plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_QC, 'r-'); 
plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_WM, 'r-'); 
% plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Hourly Binned WM Fluxes');

% clc
% figure
% hold on; box on;
% NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_WM_QC, 'r-'); 
% plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Hourly Binned WM Fluxes Full QC');

%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
close all

figure
hold on; box on;
boxplot(UAug1_BMSbin.NEP_QC)
ylabel('NEP')
title('Aug1 NEP Boxplots')

figure
hold on; box on;
boxplot(UAug1_BMSbin.NEC_QC)
ylabel('NEC')
title('Aug1 NEC Boxplots')

figure
hold on; box on;
boxplot(UAug1_BMSbin.dDOXY_QC)
ylabel('DO')
title('Aug1 dDO Boxplots')

figure
hold on; box on;
boxplot(UAug1_BMSbin.dTA_QC)
ylabel('TA')
title('Aug1 dTA Boxplots')

% Outliers:   

% 120: NEC outlier (value: -46) and dTA outlier (value: -12) 
    %  OUTLIER 
    
% 119: NEC outlier (value: -14.2) and dTA outlier (value: -4.59)
    %  OUTLIER 
    
% 117: NEC outlier (value: -15.7)
    %'09-Aug-2020 20:00:00'
    % profile good - not an outlier 

% 116: NEC outlier (value: -8.3)
    % profile good - not an outlier 

% 111: NEP outlier (vaue: 17.7) and NEC outlier (value: 6.59) - check SDN 
    % profile good - not an outlier 
    % '09-Aug-2020 14:00:00' - 2pm - values make sense

% 108: NEC outlier (value: 9.23)
    % profile good - not an outlier 

% 92: NEC outlier (value: 9.36)
    % profile good - not an outlier 

% 54: NEC outlier (value: -8.8)
    % profile good - not an outlier 


% plot to check:  
datestr(UAug1_BMSbin.SDN(111))

close all 
for i =54 %1:length(UAug1_BMSbin.SDN)
    figure (i)
    scatter(UAug1_BMSbin.uv(1:108,i), UAug1_BMSbin.bin_depth(1:108))
    title(['Cudjoe Aug1 Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end


%remove outliers: 120 119

% UAug1_BMSbin.NEP_QC(120) = NaN;
% UAug1_BMSbin.NEC_QC(120) = NaN;
% UAug1_BMSbin.dDOXY_QC(120) = NaN;
% UAug1_BMSbin.dTA_QC(120) = NaN;
% % 
% UAug1_BMSbin.NEP_QC(119) = NaN;
% UAug1_BMSbin.NEC_QC(119) = NaN;
% UAug1_BMSbin.dDOXY_QC(119) = NaN;
% UAug1_BMSbin.dTA_QC(119) = NaN;

% 
close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug1 Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe Aug1 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_QC, 'r-'); 
plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug1 Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Aug1 2020 Hourly Binned Fluxes Full QC');
% 
%% Plot Diel Curves
close all
% nepdbin = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP, 24);
% necdbin = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC, 24);
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

nepdbin_QC = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_QC, 24);
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
% nepdbin_WM = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_WM, 24);
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

nepdbin_WM_QC = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_WM_QC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM_QC, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM_QC, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM_QC,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM_QC,1), 'ro-')
% xlim([1 24]);
% ylabel(['NEP or \color{red}NEC']);
% title('Aug 2020 WM Full QC');
% xlabel('hour of day');
% 


%% Extract Daytime data for Ratios
clc
% Extract daytime data using UAug1_BMSbin.PAR
UAug1_inight = UAug1_BMSbin.PAR(1,:) < 1; %find all nightime datapoints 

%create new arrays for daytime data
UAug1_BMSbin.SDN_day = UAug1_BMSbin.SDN;
UAug1_BMSbin.PAR_day = UAug1_BMSbin.PAR(1,:);

UAug1_BMSbin.NEP_day = UAug1_BMSbin.NEP;
UAug1_BMSbin.NEC_day = UAug1_BMSbin.NEC;
UAug1_BMSbin.NEP_day_QC = UAug1_BMSbin.NEP_QC;
UAug1_BMSbin.NEC_day_QC = UAug1_BMSbin.NEC_QC;

UAug1_BMSbin.dDOXY_day_QC = UAug1_BMSbin.dDOXY_QC;      %DO Gradient
UAug1_BMSbin.dTA_day_QC = UAug1_BMSbin.dTA_QC;          %TA Gradient

UAug1_BMSbin.NEP_WM_day = UAug1_BMSbin.NEP_WM;
UAug1_BMSbin.NEC_WM_day = UAug1_BMSbin.NEC_WM;
UAug1_BMSbin.NEP_WM_day_QC = UAug1_BMSbin.NEP_WM_QC;
UAug1_BMSbin.NEC_WM_day_QC = UAug1_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
UAug1_BMSbin.SDN_day(UAug1_inight) = NaN;
UAug1_BMSbin.PAR_day (UAug1_inight) = NaN;

UAug1_BMSbin.NEP_day(UAug1_inight) = NaN;
UAug1_BMSbin.NEC_day(UAug1_inight) = NaN;
UAug1_BMSbin.NEP_day_QC(UAug1_inight) = NaN;
UAug1_BMSbin.NEC_day_QC(UAug1_inight) = NaN;

UAug1_BMSbin.dDOXY_day_QC(UAug1_inight) = NaN;      %DO Gradient
UAug1_BMSbin.dTA_day_QC(UAug1_inight) = NaN;        %TA Gradient

UAug1_BMSbin.NEP_WM_day(UAug1_inight) = NaN;
UAug1_BMSbin.NEC_WM_day(UAug1_inight) = NaN;
UAug1_BMSbin.NEP_WM_day_QC(UAug1_inight) = NaN;
UAug1_BMSbin.NEC_WM_day_QC(UAug1_inight) = NaN;

%Plot to check only nighttime points removed
figure 
hold on
scatter(UAug1_BMSbin.SDN, UAug1_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(UAug1_BMSbin.SDN_day, UAug1_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
UAug1_BMSbin.NEP_day(isnan(UAug1_BMSbin.NEP_day))=[];
UAug1_BMSbin.NEC_day(isnan(UAug1_BMSbin.NEC_day))=[];
UAug1_BMSbin.NEP_day_QC(isnan(UAug1_BMSbin.NEP_day_QC))=[];
UAug1_BMSbin.NEC_day_QC(isnan(UAug1_BMSbin.NEC_day_QC))=[];

UAug1_BMSbin.dDOXY_day_QC(isnan(UAug1_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
UAug1_BMSbin.dTA_day_QC(isnan(UAug1_BMSbin.dTA_day_QC))=[];       %TA Gradient

UAug1_BMSbin.NEP_WM_day(isnan(UAug1_BMSbin.NEP_WM_day))=[];
UAug1_BMSbin.NEC_WM_day(isnan(UAug1_BMSbin.NEC_WM_day))=[];
UAug1_BMSbin.NEP_WM_day_QC(isnan(UAug1_BMSbin.NEP_WM_day_QC))=[];
UAug1_BMSbin.NEC_WM_day_QC(isnan(UAug1_BMSbin.NEC_WM_day_QC))=[];

% create nighttime hours datasets
UAug1_BMSbin.SDN_night = UAug1_BMSbin.SDN;
UAug1_BMSbin.PAR_night = UAug1_BMSbin.PAR(1,:);

UAug1_BMSbin.NEP_night = UAug1_BMSbin.NEP;
UAug1_BMSbin.NEC_night = UAug1_BMSbin.NEC;
UAug1_BMSbin.NEP_night_QC = UAug1_BMSbin.NEP_QC;
UAug1_BMSbin.NEC_night_QC = UAug1_BMSbin.NEC_QC;

UAug1_BMSbin.dDOXY_night_QC = UAug1_BMSbin.dDOXY_QC;      %DO Gradient
UAug1_BMSbin.dTA_night_QC = UAug1_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
UAug1_BMSbin.SDN_night=UAug1_BMSbin.SDN_night(UAug1_inight);
UAug1_BMSbin.PAR_night=UAug1_BMSbin.PAR_night(UAug1_inight);

UAug1_BMSbin.NEP_night_QC=UAug1_BMSbin.NEP_night_QC(UAug1_inight);
UAug1_BMSbin.NEC_night_QC=UAug1_BMSbin.NEC_night_QC(UAug1_inight);

UAug1_BMSbin.dDOXY_night_QC=UAug1_BMSbin.dDOXY_night_QC(UAug1_inight);      %DO Gradient
UAug1_BMSbin.dTA_night_QC=UAug1_BMSbin.dTA_night_QC(UAug1_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(UAug1_BMSbin.SDN, UAug1_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(UAug1_BMSbin.SDN_night, UAug1_BMSbin.PAR_night, 'r.'); % night plot

%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(UAug1_BMSbin.NEP_day,UAug1_BMSbin.NEC_day);
% UAug1_BMSbin.Reg_Line = m*UAug1_BMSbin.NEP_day + b;
% UAug1_BMSbin.Ratio = m;
% UAug1_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(UAug1_BMSbin.NEP_day,UAug1_BMSbin.NEC_day,'o')
% plot(UAug1_BMSbin.NEP_day,UAug1_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Aug 2020 Pre-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug1_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug1_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(UAug1_BMSbin.NEP_day_QC,UAug1_BMSbin.NEC_day_QC);
UAug1_BMSbin.Reg_Line_QC = m_QC*UAug1_BMSbin.NEP_day_QC + b_QC;
UAug1_BMSbin.Ratio_QC = m_QC;
UAug1_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(UAug1_BMSbin.NEP_day_QC,UAug1_BMSbin.NEC_day_QC,'o')
plot(UAug1_BMSbin.NEP_day_QC,UAug1_BMSbin.Reg_Line_QC,'r')
%ylim([-50 50])
%xlim([-25 25])
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Aug 2020 Pre-Restoration NCC:NCP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug1_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug1_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(UAug1_BMSbin.NEP_WM_day,UAug1_BMSbin.NEC_WM_day);
% UAug1_BMSbin.Reg_Line_WM = m_WM*UAug1_BMSbin.NEP_WM_day + b_WM;
% UAug1_BMSbin.Ratio_WM = m_WM;
% UAug1_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(UAug1_BMSbin.NEP_WM_day,UAug1_BMSbin.NEC_WM_day,'o')
% plot(UAug1_BMSbin.NEP_WM_day,UAug1_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Aug 2020 Pre-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug1_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug1_BMSbin.R2_WM)
% 

[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(UAug1_BMSbin.NEP_WM_day_QC,UAug1_BMSbin.NEC_WM_day_QC);
UAug1_BMSbin.Reg_Line_WM_QC = m_WM_QC*UAug1_BMSbin.NEP_WM_day_QC + b_WM_QC;
UAug1_BMSbin.Ratio_WM_QC = m_WM_QC;
UAug1_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(UAug1_BMSbin.NEP_WM_day_QC,UAug1_BMSbin.NEC_WM_day_QC,'o')
% plot(UAug1_BMSbin.NEP_WM_day_QC,UAug1_BMSbin.Reg_Line_WM_QC,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Aug 2020 Pre-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug1_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug1_BMSbin.R2_WM_QC)

%% For NEC:NEP Regressions Using Gradients
close all 
clc

% multiply o2 gradient by -1 for O2 production
UAug1_BMSbin.dDOXY_Reg = -1.*UAug1_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
UAug1_BMSbin.dTA_Reg = 0.5.*UAug1_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(UAug1_BMSbin.dDOXY_day_QC);
figure 
hold on 
DOday = plot(1:Xlength, UAug1_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, UAug1_BMSbin.dDOXY_Reg);
xlabel('Aug1 Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Cudjoe Aug1 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, UAug1_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, UAug1_BMSbin.dTA_Reg);
xlabel('Aug1 Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Cudjoe Aug1 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(UAug1_BMSbin.dDOXY_Reg, UAug1_BMSbin.dTA_Reg);
UAug1_BMSbin.Reg_Line_G = m_G*UAug1_BMSbin.dDOXY_Reg + b_G;
UAug1_BMSbin.Ratio_G = m_G;
UAug1_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(UAug1_BMSbin.dDOXY_Reg,UAug1_BMSbin.dTA_Reg,'o')
plot(UAug1_BMSbin.dDOXY_Reg,UAug1_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Aug1 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug1_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug1_BMSbin.R2_G)

clc
disp('Finished with Aug 1 Metabolism Calculations');
save('UAug120_2.mat', 'UAug1_BMS', 'UAug1_BMSbin');

%% Plot Profiles 
% plot for profile within pump heights
% close all 
% for i =1:100 %length(UAug1_ADavg.SDN)
%     figure (i)
%     scatter(UAug1_ADavg.uv(1:108,i), UAug1_ADavg.bin_depth(1:108))
%     title(['Cudjoe Aug Velocity Profile Number ',num2str(i),])
%     xlabel('Velocity (m/s)');
%     ylabel('Height (m)');
% end


%% Subplots 
close all
clc


sgtitle('Cudjoe Aug1 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'08/05 12:00';'08/06 12:00';'08/07 12:00';'08/08 12:00';'08/09 12:00';...
    '08/10 12:00';'08/11 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
% legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'southwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(UAug1_BMSbin.SDN, zeros(size(UAug1_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug1_good_Xrange, 'XTick', UAug1_tick, 'xticklabel', UAug1_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'08/05 12:00';'08/06 12:00';'08/07 12:00';'08/08 12:00';'08/09 12:00';...
    '08/10 12:00';'08/11 12:00'})
xlabel('Aug1 Days');
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(UAug1_BMSbin.SDN, UAug1_BMSbin.NEC_QC, 24);
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
plot(UAug1_BMSbin.NEP_day_QC,UAug1_BMSbin.NEC_day_QC,'o')
plot(UAug1_BMSbin.NEP_day_QC,UAug1_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
ylim([-30 30]);
xlim([-30 30]);
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(UAug1_BMSbin.Ratio_QC,2);
str2 = num2str(UAug1_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(UAug1_BMSbin.dDOXY_Reg,UAug1_BMSbin.dTA_Reg,'o')
plot(UAug1_BMSbin.dDOXY_Reg,UAug1_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
ylim([-4 4]);
xlim([-4 4]);
title('NEC:NEP Ratio from Gradients');
str1 = num2str(UAug1_BMSbin.Ratio_G,2);
str2 = num2str(UAug1_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')



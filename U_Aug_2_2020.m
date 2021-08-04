% Aug_2 SeapHOx Data Analysis from Cudjoe Ledge Reef
% Michelle Platz - USF 
% 3/2/2021

% SeapHOx sensor deployed 8/12/2020 at 9:15 am EST
    % SP datafile: 'U_812BMS.txt'
    % pump 1 height above benthos = 70 cm  
    % pump 2 height above benthos = 20 cm 
% ADP sensor deployed 8/05/2020 at 9:15 am EST
    % ADCP datafile: 'U_8_05'
    % height from substrate to ADCP head = 18 cm

close all
clc
clear all
%% Initial look at data
% ***** create UAug2_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
UAug2_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('U_812BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
UAug2_SPraw.DOXY = UAug2_SPraw.O2SATPER.*calcO2sat(UAug2_SPraw.MCAT_TC, UAug2_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
UAug2_SPraw.pHint_prelim = calc_dfet_pHint(UAug2_SPraw.Vint, UAug2_SPraw.DFET_TC, -0.4);

% ***** create UAug2_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
UAug2_SP = parse_to_pumpheights_ARM_2pump_Mar19(UAug2_SPraw, [datenum('08-12-2020 13:58:00'), datenum('09-1-2020 17:27:30')]);

% Calculate Gradients 
UAug2_SP = calc_TA_gradientV2(UAug2_SP, 2369.19, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

UAug2_SP.dDOXY = UAug2_SP.DOXY(1,:) - UAug2_SP.DOXY(2,:); %Oxygen Gradient
UAug2_SP.dpH = UAug2_SP.pH(1,:) - UAug2_SP.pH(2,:); %pH Gradient 
UAug2_SP.dTA = UAug2_SP.TAtop - UAug2_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
UAug2_DateString = {'08/12/2020 12:00:00';'08/13/2020 12:00:00';'08/14/2020 12:00:00';'08/15/2020 12:00:00';'08/16/2020 12:00:00';...
    '08/17/2020 12:00:00';'08/18/2020 12:00:00';'08/19/2020 12:00:00';'08/20/2020 12:00:00';'08/21/2020 12:00:00';'08/22/2020 12:00:00';...
    '08/23/2020 12:00:00';'08/24/2020 12:00:00';'08/25/2020 12:00:00';'08/26/2020 12:00:00';'08/27/2020 12:00:00';'08/28/2020 12:00:00';...
    '08/29/2020 12:00:00';'08/30/2020 12:00:00';'08/31/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
UAug2_tick = datenum(UAug2_DateString,formatIn);

UAug2_Xrange = [datenum('08-12-2020 13:58:00'), datenum('09-1-2020 17:27:30')];

figure
hold on; box on;
plot(UAug2_SP.SDN, UAug2_SP.dDOXY); %oxygen gradient 
plot(UAug2_SP.SDN, UAug2_SP.dTA); %TA gradient
plot(UAug2_SP.SDN, zeros(size(UAug2_SP.SDN))); %zero line
set(gca, 'xlim', UAug2_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe Aug 2020 Unbinned DO and TA Gradients');


%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
UAug2_SiteChar.SDN = UAug2_SPraw.SDN;
UAug2_SiteChar.TC = UAug2_SPraw.OPT_TC;
UAug2_SiteChar.PAR = UAug2_SPraw.PAR;
UAug2_SiteChar.PSAL = UAug2_SPraw.PSAL;
UAug2_SiteChar.Pres = UAug2_SPraw.Pres;

close all
figure 
hold on; 
plot(UAug2_SiteChar.SDN, UAug2_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
UAug2_SiteChar.SDN = UAug2_SPraw.SDN(16:57138);
UAug2_SiteChar.TC = UAug2_SPraw.OPT_TC(16:57138);
UAug2_SiteChar.PAR = UAug2_SPraw.PAR(16:57138);
UAug2_SiteChar.PSAL = UAug2_SPraw.PSAL(16:57138);
UAug2_SiteChar.Pres = UAug2_SPraw.Pres(16:57138);

%extract full length of ADCP datafile  --> already done in Aug2 script
% % % UAug2_ADfull=aquadoppraw2mat('U_8_05', 70, [datenum('06-05-2020 08:00:00'), datenum('09-01-2020 08:00:00')]);
% % % 
% % % % add AD variables to Site Char 
% % % UAug2_SiteChar.AD_SDN = UAug2_ADfull.SDN;
% % % UAug2_SiteChar.AD_Pres = UAug2_ADfull.Pres;
% % % UAug2_SiteChar.AD_TC = UAug2_ADfull.TC;
% % % UAug2_SiteChar.bin_depth = UAug2_ADfull.bin_depth;
% % % UAug2_SiteChar.u = UAug2_ADfull.u;
% % % UAug2_SiteChar.v = UAug2_ADfull.v;
% % % UAug2_SiteChar.w = UAug2_ADfull.w;
% % % UAug2_SiteChar.uv = UAug2_ADfull.uv;
% % % UAug2_SiteChar.direction = UAug2_ADfull.direction;
% % % 
% % % %Plot to see when surface interval observations are
% % % close all
% % % figure 
% % % hold on; 
% % % plot(UAug2_SiteChar.AD_SDN, UAug2_SiteChar.AD_Pres)
% % % 
% % % %clip ends of data to remove surfave interval observations 
% % % UAug2_SiteChar.AD_SDN = UAug2_ADfull.SDN(95:end);
% % % UAug2_SiteChar.AD_Pres = UAug2_ADfull.Pres(95:end);
% % % UAug2_SiteChar.AD_TC = UAug2_ADfull.TC(95:end);
% % % UAug2_SiteChar.bin_depth = UAug2_ADfull.bin_depth;
% % % UAug2_SiteChar.u = UAug2_ADfull.u(:,95:end);
% % % UAug2_SiteChar.v = UAug2_ADfull.v(:,95:end);
% % % UAug2_SiteChar.w = UAug2_ADfull.w(:,95:end);
% % % UAug2_SiteChar.uv = UAug2_ADfull.uv(:,95:end);
% % % UAug2_SiteChar.direction = UAug2_ADfull.direction(:,95:end);
% % % 
% % % % find U0 
% % % UAug2_z1 = 0.70;
% % % UAug2_z2 = 0.20;
% % % UAug2_ADheight = 0.18;
% % % UAug2_ADbin_depth_1m = 1-(UAug2_ADheight);% = 0.82
% % % UAug2_i1m = find(UAug2_SiteChar.bin_depth==(0.82));
% % % UAug2_SiteChar.U0 = UAug2_SiteChar.uv(UAug2_i1m,:);

% save data in separate datastructure
save('UAug220_SiteChar_2.mat', 'UAug2_SiteChar' )


%% Constrain Xrange from graph results and extract good gradient data - 
close all 

UAug2_good_Xrange = [datenum('08-12-2020 14:58:00'), datenum('08-21-2020 13:00:00')];

% plot to check range is correct
figure
hold on; box on;
plot(UAug2_SP.SDN, UAug2_SP.dDOXY); %oxygen gradient 
plot(UAug2_SP.SDN, UAug2_SP.dTA); %TA gradient
plot(UAug2_SP.SDN, zeros(size(UAug2_SP.SDN))); %zero line
set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Cudjoe Aug 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all

U_Aug_2_BMS_idx_start = find(UAug2_SP.SDN==datenum('08-12-2020 14:58:00'))
U_Aug_2_BMS_idx_end = find(UAug2_SP.SDN==datenum('08-21-2020 13:58:00'))

% Create new data vectors of just the good data
U_Aug_2_BMS_good_data = U_Aug_2_BMS_idx_start:U_Aug_2_BMS_idx_end;
Initial_data_points = length(U_Aug_2_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(UAug2_SP);
for v = 1:length(vars)
    UAug2_SP.(vars{v}) = (UAug2_SP.(vars{v})(:,U_Aug_2_BMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: UAug2_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
UAug2_AD=aquadoppraw2mat('U_8_05', 70, [datenum('08-12-2020 14:58:00'), datenum('08-21-2020 14:13:00')]);

%averages data to the middle of the minute interval spacified 
UAug2_ADavg = average_aquadopp(UAug2_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.7m (pump 1) and 0.2m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.52m (bin 42) and 0.02m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
UAug2_Sal_est = mean(UAug2_SP.PSAL(1,3:end));

[UAug2_ADavg] = ustar_from_aquadopp2(UAug2_ADavg,[0.52 0.11], UAug2_Sal_est); %bins adjusted 

clc
%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[UAug2_ADavg] = ustar_McGillis_Method(UAug2_ADavg, 0.70, 0.20, 42, 1);

%compare
close all
figure
hold on 
ustar_plot = plot(UAug2_ADavg.SDN, UAug2_ADavg.ustar, 'r');
ustar_WM_plot = plot(UAug2_ADavg.SDN, UAug2_ADavg.ustar_WM);
plot(UAug2_ADavg.SDN, zeros(size(UAug2_ADavg.SDN)),'k');
set(gca, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug2_ADavg.SDN(1) UAug2_ADavg.SDN(end)])
xlabel('Aug Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: UAug2_BMS *****

ADavg_vars = fieldnames(UAug2_ADavg);
for v = 1:length(ADavg_vars)
    UAug2_BMS.(ADavg_vars{v}) = (UAug2_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(UAug2_SP);
for v = 1:length(SP_vars)
    UAug2_BMS.(SP_vars{v}) = (UAug2_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(UAug2_BMS.SDN)
% min: 13-28-43-58 becuase SP was restarted in the field rather than on the
% minute
%% %% *************** Calculate Fluxes ****************

% actual pump heights  = 0.70m (pump 1) and 0.20m (pump 2) 
% 0.18m from substrate to ACDP head in Aug at U 
% adjusted height = 0.52 m (bin 42) and 0.02 m (too shallow) (bin 1)  
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[UAug2_BMS] = calc_NCC_3(UAug2_BMS,[0.52 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(UAug2_BMS.DOXY(1,:))
[UAug2_BMS] = calc_NCP_3(UAug2_BMS, [0.52 0.11],C1guess); %estimate C1 guess using median DOXY(1,:) value  

% Plot NCP and NCC
%close all
figure
hold on; box on; 
NEPplot = plot(UAug2_BMS.SDN, UAug2_BMS.NEP);
NECplot = plot(UAug2_BMS.SDN, UAug2_BMS.NEC);
plot(UAug2_BMS.SDN, zeros(size(UAug2_BMS.SDN)),'k');
set(gca, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug2_BMS.SDN(1) UAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');

%McGillis method flux calculations 
[UAug2_BMS] = calc_NCP_McGillis_Method(UAug2_BMS, 0.70, 0.20, UAug2_Sal_est);
[UAug2_BMS] = calc_NCC_McGillis_Method(UAug2_BMS, 0.70, 0.20, UAug2_Sal_est);

% Plot NCP and NCC
% close all
% figure
% hold on; box on; 
% NEPplot = plot(UAug2_BMS.SDN, UAug2_BMS.NEP_WM);
% NECplot = plot(UAug2_BMS.SDN, UAug2_BMS.NEC_WM);
% plot(UAug2_BMS.SDN, zeros(size(UAug2_BMS.SDN)),'k');
% set(gca, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlim([UAug2_BMS.SDN(1) UAug2_BMS.SDN(end)])
% xlabel('Aug Days');
% ylabel('NCP or NCC [mmol/m2/hr]');
% legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 WM Fluxes');

%% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(UAug2_BMS.SDN, UAug2_BMS.NEP);
NEPplotWM = plot(UAug2_BMS.SDN, UAug2_BMS.NEP_WM);
plot(UAug2_BMS.SDN, zeros(size(UAug2_BMS.SDN)),'k');
set(gca, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug2_BMS.SDN(1) UAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(UAug2_BMS.SDN, UAug2_BMS.NEC);
NECplotWM = plot(UAug2_BMS.SDN, UAug2_BMS.NEC_WM);
plot(UAug2_BMS.SDN, zeros(size(UAug2_BMS.SDN)),'k');
set(gca, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug2_BMS.SDN(1) UAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
UAug2_BMS.DOXYstd = std(UAug2_BMS.DOXY);

% get DOXY std
UAug2_idoxystd = find(UAug2_BMS.DOXYstd > 2);
UAug2_ihighdoxystd2 = [];
for i = 1:length(UAug2_idoxystd)
    
    UAug2_ihighdoxystd2 = vertcat(UAug2_ihighdoxystd2,[UAug2_idoxystd(i)-1:1:UAug2_idoxystd(i)+1]');
end
% get unique IDs
UAug2_ihighdoxystd2 = unique(UAug2_ihighdoxystd2);
% remove 0's and out of index values
UAug2_ihighdoxystd2(UAug2_ihighdoxystd2==0) = [];
UAug2_ihighdoxystd2(UAug2_ihighdoxystd2> length(UAug2_BMS.SDN)) = [];

% make it into index
trex = false(size(UAug2_BMS.SDN));
trex(UAug2_ihighdoxystd2) = true;
UAug2_ihighdoxystd2 = trex;
clear trex;

UAug2_BMS.NEP_QC = UAug2_BMS.NEP;
UAug2_BMS.NEC_QC = UAug2_BMS.NEC;
UAug2_BMS.dDOXY_QC = UAug2_BMS.dDOXY; %DO gradient
UAug2_BMS.dTA_QC = UAug2_BMS.dTA;     %TA gradient
UAug2_BMS.NEP_WM_QC = UAug2_BMS.NEP_WM;
UAug2_BMS.NEC_WM_QC = UAug2_BMS.NEC_WM;

% set observations when DOXYstd>0.8 to NaN
UAug2_BMS.NEP_QC(UAug2_ihighdoxystd2) = NaN;
UAug2_BMS.NEC_QC(:,UAug2_ihighdoxystd2) = NaN;
UAug2_BMS.dDOXY_QC(UAug2_ihighdoxystd2) = NaN; %DO gradient
UAug2_BMS.dTA_QC(:,UAug2_ihighdoxystd2) = NaN;   %TA gradient
UAug2_BMS.NEP_WM_QC(UAug2_ihighdoxystd2) = NaN;
UAug2_BMS.NEC_WM_QC(:,UAug2_ihighdoxystd2) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(UAug2_BMS.SDN, UAug2_BMS.NEP, 'k');
NEPplotQC = plot(UAug2_BMS.SDN, UAug2_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(UAug2_BMS.SDN, zeros(size(UAug2_BMS.SDN)),'k');
set(gca, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug2_BMS.SDN(1) UAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');


%WM Plot
%close all
figure
hold on; box on;
NEPplot = plot(UAug2_BMS.SDN, UAug2_BMS.NEP_WM, 'k');
NEPplotQC = plot(UAug2_BMS.SDN, UAug2_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
plot(UAug2_BMS.SDN, zeros(size(UAug2_BMS.SDN)),'k');
set(gca, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([UAug2_BMS.SDN(1) UAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Fluxes');


%% Bin data to hourly intervals

X = floor(nanmin(UAug2_BMS.SDN)):1/24:ceil(nanmax(UAug2_BMS.SDN));

UAug2_BMSbin.SDN = X;
% variables from 2 differnet heights
UAug2_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.DOXY(1,:), X);
UAug2_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.DOXY(2,:), X);

UAug2_BMSbin.pH(1,:)       = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.pH(1,:), X);
UAug2_BMSbin.pH(2,:)       = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.pH(2,:), X);

UAug2_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.PSAL(1,:), X);
UAug2_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.PSAL(2,:), X);

UAug2_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.O2SATPER(1,:), X);
UAug2_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.O2SATPER(2,:), X);

UAug2_BMSbin.Pres(1,:)     = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.Pres(1,:), X);
UAug2_BMSbin.Pres(2,:)     = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.Pres(2,:), X);

UAug2_BMSbin.DENS(1,:)     = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.DENS(1,:), X);
UAug2_BMSbin.DENS(2,:)     = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.DENS(2,:), X);

UAug2_BMSbin.PAR(1,:)      = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.PAR(1,:), X);
UAug2_BMSbin.PAR(2,:)      = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.PAR(2,:), X);


UAug2_BMSbin.bin_depth     = UAug2_BMS.bin_depth;
for i = 1:108
    UAug2_BMSbin.uv(i,:)   = bin_data_to_X_GF(UAug2_BMS.SDN,UAug2_BMS.uv(i,:), X);
end


% bin data hourly. Vector variables 
UAug2_BMSbin.PRES  = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.Pres, X);
UAug2_BMSbin.U0       = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.U0, X);
UAug2_BMSbin.DIR      = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.direction, X);
UAug2_BMSbin.ustar    = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.ustar, X);
UAug2_BMSbin.ustar_rm = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.ustar_runmean, X);
UAug2_BMSbin.dTA      = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.dTA, X);
UAug2_BMSbin.dTA_QC   = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.dTA_QC, X);
UAug2_BMSbin.dpH      = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.dpH, X);
UAug2_BMSbin.dDOXY    = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.dDOXY, X);
UAug2_BMSbin.dDOXY_QC = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.dDOXY_QC, X);
UAug2_BMSbin.NEP      = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.NEP, X);
UAug2_BMSbin.NEP_QC   = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.NEP_QC, X);
UAug2_BMSbin.NEC      = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.NEC, X);
UAug2_BMSbin.NEC_QC   = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.NEC_QC, X);
UAug2_BMSbin.NEP_WM   = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.NEP_WM, X);
UAug2_BMSbin.NEP_WM_QC= bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.NEP_WM_QC, X);
UAug2_BMSbin.NEC_WM   = bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.NEC_WM, X);
UAug2_BMSbin.NEC_WM_QC= bin_data_to_X_GF(UAug2_BMS.SDN, UAug2_BMS.NEC_WM_QC, X);

%% Plot Binned Fluxes 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP, 'b'); 
% NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC, 'r-'); 
% plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Hourly Binned Fluxes FluxFit');

% close all
clc
figure
hold on; box on;
NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_QC, 'r-'); 
plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Cudjoe Aug 2020 Hourly Binned Fluxes FluxFit QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_WM, 'r-'); 
% plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 WM Original Hourly Binned Fluxes');

% close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_WM_QC, 'r-'); 
% plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 WM Full QC Hourly Binned Fluxes');

%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
DOplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.dDOXY, 'c'); 
TAplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.dTA, 'k-'); 
plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug2 Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe Aug2 2020 Binned Gradiets');



%% ***QC*** Remove sections when velocity is too slow
ibad = UAug2_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow

% UAug2_BMSbin.NEP(ibad) = NaN;
% UAug2_BMSbin.NEC(ibad) = NaN;
UAug2_BMSbin.NEP_QC(ibad) = NaN;
UAug2_BMSbin.NEC_QC(ibad) = NaN;
UAug2_BMSbin.dDOXY_QC(ibad) = NaN;
UAug2_BMSbin.dTA_QC(ibad) = NaN;

% UAug2_BMSbin.NEP_WM(ibad) = NaN;
% UAug2_BMSbin.NEC_WM(ibad) = NaN;
UAug2_BMSbin.NEP_WM_QC(ibad) = NaN;
UAug2_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP, 'b'); 
% NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC, 'r-'); 
% plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_QC, 'r-'); 
plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Aug 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_WM, 'r-'); 
% plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Hourly Binned WM Fluxes');

% clc
% figure
% hold on; box on;
% NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_WM_QC, 'r-'); 
% plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Cudjoe Aug 2020 Hourly Binned WM Fluxes Full QC');

%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
figure
hold on; box on;
boxplot(UAug2_BMSbin.NEP_QC)
ylabel('NEP')
title('Aug2 NEP Boxplots')

figure
hold on; box on;
boxplot(UAug2_BMSbin.NEC_QC)
ylabel('NEC')
title('Aug2 NEC Boxplots')

figure
hold on; box on;
boxplot(UAug2_BMSbin.dDOXY_QC)
ylabel('DO')
title('Aug2 dDO Boxplots')

figure
hold on; box on;
boxplot(UAug2_BMSbin.dTA_QC)
ylabel('TA')
title('Aug2 dTA Boxplots')

%Outliers 

% 162: NEP outlier (value: -55.1)
    % profile good 
    %'18-Aug-2020 17:00:00' - daytime point, value unlikely  
    % -  point is OUTLIER

% 163: NEP outlier (value: -37.7)
    % profile good 
    % daytime point, value unlikely 
    % -  point is OUTLIER

% 186: NEP outlier (value: -43.8)
    % '19-Aug-2020 17:00:00' - daytime point, value unlikely 
    % bad profile -  point is OUTLIER

% 187: NEC and dTA outlier 
    % point is OUTLIER
    
% 125: dTA outlier (value: -2.14) 
    %not an outlier

% 93: NEC outlier (value: -5.38) 
    %not an outlier

% 211:  NEP outlier (value: -36.2)
    % Bad profile - OUTLIER 

datestr(UAug2_BMSbin.SDN(162))

% Plot Profiles at outliers - 187 115 162 163 164
close all 
for i = 211   %1:length(UAug2_BMSbin.SDN)
    figure (i)
    scatter(UAug2_BMSbin.uv(1:108,i), UAug2_BMSbin.bin_depth(1:108))
    title(['Cudjoe Aug2 Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end
%outliers to be removed:  
% 162, 163, 186 187, 211


% UAug2_BMSbin.NEP_QC(162) = NaN;
% UAug2_BMSbin.NEC_QC(162) = NaN;
% UAug2_BMSbin.dDOXY_QC(162) = NaN;
% UAug2_BMSbin.dTA_QC(162) = NaN;
% 
% UAug2_BMSbin.NEP_QC(163) = NaN;
% UAug2_BMSbin.NEC_QC(163) = NaN;
% UAug2_BMSbin.dDOXY_QC(163) = NaN;
% UAug2_BMSbin.dTA_QC(163) = NaN;
% 
% UAug2_BMSbin.NEP_QC(186) = NaN;
% UAug2_BMSbin.NEC_QC(186) = NaN;
% UAug2_BMSbin.dDOXY_QC(186) = NaN;
% UAug2_BMSbin.dTA_QC(186) = NaN;
% 
% UAug2_BMSbin.NEP_QC(187) = NaN;
% UAug2_BMSbin.NEC_QC(187) = NaN;
% UAug2_BMSbin.dDOXY_QC(187) = NaN;
% UAug2_BMSbin.dTA_QC(187) = NaN;
% 
% UAug2_BMSbin.NEP_QC(211) = NaN;
% UAug2_BMSbin.NEC_QC(211) = NaN;
% UAug2_BMSbin.dDOXY_QC(211) = NaN;
% UAug2_BMSbin.dTA_QC(211) = NaN;
% 


close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug2 Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Cudjoe Aug2 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_QC, 'b'); 
NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_QC, 'r-'); 
plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug2 Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Cudjoe Aug2 2020 Hourly Binned Fluxes Full QC');


%% Plot Diel Curves

% nepdbin = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP, 24);
% necdbin = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC, 24);
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

nepdbin_QC = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_QC, 24);
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
title('Aug Diel Plot 2020 Fluxfit QC');
xlabel('hour of day');

% WM data
% nepdbin_WM = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_WM, 24);
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

nepdbin_WM_QC = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_WM_QC, 24);
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



%% Extract Daytime data for Ratios
clc
% Extract daytime data using UAug2_BMSbin.PAR
UAug2_inight = UAug2_BMSbin.PAR(1,:) < 10; %find all nightime datapoints 

%create new arrays for daytime data
UAug2_BMSbin.SDN_day = UAug2_BMSbin.SDN;
UAug2_BMSbin.PAR_day = UAug2_BMSbin.PAR(1,:);

UAug2_BMSbin.NEP_day = UAug2_BMSbin.NEP;
UAug2_BMSbin.NEC_day = UAug2_BMSbin.NEC;
UAug2_BMSbin.NEP_day_QC = UAug2_BMSbin.NEP_QC;
UAug2_BMSbin.NEC_day_QC = UAug2_BMSbin.NEC_QC;

UAug2_BMSbin.dDOXY_day_QC = UAug2_BMSbin.dDOXY_QC;      %DO Gradient
UAug2_BMSbin.dTA_day_QC = UAug2_BMSbin.dTA_QC;          %TA Gradient

UAug2_BMSbin.NEP_WM_day = UAug2_BMSbin.NEP_WM;
UAug2_BMSbin.NEC_WM_day = UAug2_BMSbin.NEC_WM;
UAug2_BMSbin.NEP_WM_day_QC = UAug2_BMSbin.NEP_WM_QC;
UAug2_BMSbin.NEC_WM_day_QC = UAug2_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
UAug2_BMSbin.SDN_day(UAug2_inight) = NaN;
UAug2_BMSbin.PAR_day (UAug2_inight) = NaN;

UAug2_BMSbin.NEP_day(UAug2_inight) = NaN;
UAug2_BMSbin.NEC_day(UAug2_inight) = NaN;
UAug2_BMSbin.NEP_day_QC(UAug2_inight) = NaN;
UAug2_BMSbin.NEC_day_QC(UAug2_inight) = NaN;

UAug2_BMSbin.dDOXY_day_QC(UAug2_inight) = NaN;      %DO Gradient
UAug2_BMSbin.dTA_day_QC(UAug2_inight) = NaN;        %TA Gradient

UAug2_BMSbin.NEP_WM_day(UAug2_inight) = NaN;
UAug2_BMSbin.NEC_WM_day(UAug2_inight) = NaN;
UAug2_BMSbin.NEP_WM_day_QC(UAug2_inight) = NaN;
UAug2_BMSbin.NEC_WM_day_QC(UAug2_inight) = NaN;
%Plot to check only nighttime points removed
figure 
hold on
scatter(UAug2_BMSbin.SDN, UAug2_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(UAug2_BMSbin.SDN_day, UAug2_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
UAug2_BMSbin.NEP_day(isnan(UAug2_BMSbin.NEP_day))=[];
UAug2_BMSbin.NEC_day(isnan(UAug2_BMSbin.NEC_day))=[];
UAug2_BMSbin.NEP_day_QC(isnan(UAug2_BMSbin.NEP_day_QC))=[];
UAug2_BMSbin.NEC_day_QC(isnan(UAug2_BMSbin.NEC_day_QC))=[];

UAug2_BMSbin.dDOXY_day_QC(isnan(UAug2_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
UAug2_BMSbin.dTA_day_QC(isnan(UAug2_BMSbin.dTA_day_QC))=[];       %TA Gradient

UAug2_BMSbin.NEP_WM_day(isnan(UAug2_BMSbin.NEP_WM_day))=[];
UAug2_BMSbin.NEC_WM_day(isnan(UAug2_BMSbin.NEC_WM_day))=[];
UAug2_BMSbin.NEP_WM_day_QC(isnan(UAug2_BMSbin.NEP_WM_day_QC))=[];
UAug2_BMSbin.NEC_WM_day_QC(isnan(UAug2_BMSbin.NEC_WM_day_QC))=[];

% create nighttime hours datasets
UAug2_BMSbin.SDN_night = UAug2_BMSbin.SDN;
UAug2_BMSbin.PAR_night = UAug2_BMSbin.PAR(1,:);

UAug2_BMSbin.NEP_night = UAug2_BMSbin.NEP;
UAug2_BMSbin.NEC_night = UAug2_BMSbin.NEC;
UAug2_BMSbin.NEP_night_QC = UAug2_BMSbin.NEP_QC;
UAug2_BMSbin.NEC_night_QC = UAug2_BMSbin.NEC_QC;

UAug2_BMSbin.dDOXY_night_QC = UAug2_BMSbin.dDOXY_QC;      %DO Gradient
UAug2_BMSbin.dTA_night_QC = UAug2_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
UAug2_BMSbin.SDN_night=UAug2_BMSbin.SDN_night(UAug2_inight);
UAug2_BMSbin.PAR_night=UAug2_BMSbin.PAR_night(UAug2_inight);

UAug2_BMSbin.NEP_night_QC=UAug2_BMSbin.NEP_night_QC(UAug2_inight);
UAug2_BMSbin.NEC_night_QC=UAug2_BMSbin.NEC_night_QC(UAug2_inight);

UAug2_BMSbin.dDOXY_night_QC=UAug2_BMSbin.dDOXY_night_QC(UAug2_inight);      %DO Gradient
UAug2_BMSbin.dTA_night_QC=UAug2_BMSbin.dTA_night_QC(UAug2_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(UAug2_BMSbin.SDN, UAug2_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(UAug2_BMSbin.SDN_night, UAug2_BMSbin.PAR_night, 'r.'); % night plot


%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

%Shapiro Wilk normality test: null hypothesis of composite normality
%       H = 0 => Do not reject the null hypothesis at significance level ALPHA.
%       H = 1 => Reject the null hypothesis at significance level ALPHA.
% clc
[H, pValue, SWstatistic] = swtest(UAug2_BMSbin.NEP_day_QC, 0.05)
[H, pValue, SWstatistic] = swtest(UAug2_BMSbin.NEC_day_QC, 0.05)

%Lilliefors Normal Distribution Test - returns h = 0, indicating that the null hypothesis 
%-- that the data is a sample from a normal distribution -- is not rejected.
[h,p] = lillietest(UAug2_BMSbin.NEP_day_QC)
[h,p] = lillietest(UAug2_BMSbin.NEC_day_QC)

% *************** both test results confirm non-normal datasets ************


% [m,b,r,sm,sb]=lsqfitgm(UAug2_BMSbin.NEP_day,UAug2_BMSbin.NEC_day);
% UAug2_BMSbin.Reg_Line = m*UAug2_BMSbin.NEP_day + b;
% UAug2_BMSbin.Ratio = m;
% UAug2_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(UAug2_BMSbin.NEP_day,UAug2_BMSbin.NEC_day,'o')
% plot(UAug2_BMSbin.NEP_day,UAug2_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Aug 2020 Post-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug2_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug2_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(UAug2_BMSbin.NEP_day_QC,UAug2_BMSbin.NEC_day_QC);
UAug2_BMSbin.Reg_Line_QC = m_QC*UAug2_BMSbin.NEP_day_QC + b_QC;
UAug2_BMSbin.Ratio_QC = m_QC;
UAug2_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(UAug2_BMSbin.NEP_day_QC,UAug2_BMSbin.NEC_day_QC,'o')
plot(UAug2_BMSbin.NEP_day_QC,UAug2_BMSbin.Reg_Line_QC,'r')
%ylim([-50 50])
%xlim([-25 25])
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Aug 2020 Post-Restoration NCC:NCP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug2_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug2_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(UAug2_BMSbin.NEP_WM_day,UAug2_BMSbin.NEC_WM_day);
% UAug2_BMSbin.Reg_Line_WM = m_WM*UAug2_BMSbin.NEP_WM_day + b_WM;
% UAug2_BMSbin.Ratio_WM = m_WM;
% UAug2_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(UAug2_BMSbin.NEP_WM_day,UAug2_BMSbin.NEC_WM_day,'o')
% plot(UAug2_BMSbin.NEP_WM_day,UAug2_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Aug 2020 Post-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug2_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug2_BMSbin.R2_WM)


[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(UAug2_BMSbin.NEP_WM_day_QC,UAug2_BMSbin.NEC_WM_day_QC);
UAug2_BMSbin.Reg_Line_WM_QC = m_WM_QC*UAug2_BMSbin.NEP_WM_day_QC + b_WM_QC;
UAug2_BMSbin.Ratio_WM_QC = m_WM_QC;
UAug2_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(UAug2_BMSbin.NEP_WM_day_QC,UAug2_BMSbin.NEC_WM_day_QC,'o')
% plot(UAug2_BMSbin.NEP_WM_day_QC,UAug2_BMSbin.Reg_Line_WM_QC,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Cudjoe Aug 2020 Post-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug2_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug2_BMSbin.R2_WM_QC)
% 

%% For NEC:NEP Regressions Using Gradients
close all 
clc

% multiply o2 gradient by -1 for O2 production
UAug2_BMSbin.dDOXY_Reg = -1.*UAug2_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
UAug2_BMSbin.dTA_Reg = 0.5.*UAug2_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(UAug2_BMSbin.dDOXY_day_QC);

figure 
hold on 
DOday = plot(1:Xlength, UAug2_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, UAug2_BMSbin.dDOXY_Reg);
xlabel('Aug2 Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Cudjoe Aug2 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, UAug2_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, UAug2_BMSbin.dTA_Reg);
xlabel('Aug2 Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Cudjoe Aug2 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(UAug2_BMSbin.dDOXY_Reg, UAug2_BMSbin.dTA_Reg);
UAug2_BMSbin.Reg_Line_G = m_G*UAug2_BMSbin.dDOXY_Reg + b_G;
UAug2_BMSbin.Ratio_G = m_G;
UAug2_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(UAug2_BMSbin.dDOXY_Reg,UAug2_BMSbin.dTA_Reg,'o')
plot(UAug2_BMSbin.dDOXY_Reg,UAug2_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Aug2 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + UAug2_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + UAug2_BMSbin.R2_G)

clc
disp('Finished with Aug 2 Metabolism Calculations');

save('UAug220_2.mat', 'UAug2_BMS', 'UAug2_BMSbin');



%% Subplots 
close all
clc

sgtitle('Cudjoe Aug2 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'08/12 12:00';'08/13 12:00';'08/14 12:00';'08/15 12:00';'08/16 12:00';...
    '08/17 12:00';'08/18 12:00';'08/19 12:00';'08/20 12:00';'08/21 12:00'})
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'southwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(UAug2_BMSbin.SDN, zeros(size(UAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', UAug2_good_Xrange, 'XTick', UAug2_tick, 'xticklabel', UAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug2 Days');
xticklabels({'08/12 12:00';'08/13 12:00';'08/14 12:00';'08/15 12:00';'08/16 12:00';...
    '08/17 12:00';'08/18 12:00';'08/19 12:00';'08/20 12:00';'08/21 12:00'})
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(UAug2_BMSbin.SDN, UAug2_BMSbin.NEC_QC, 24);
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
plot(UAug2_BMSbin.NEP_day_QC,UAug2_BMSbin.NEC_day_QC,'o')
plot(UAug2_BMSbin.NEP_day_QC,UAug2_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(UAug2_BMSbin.Ratio_QC,2);
str2 = num2str(UAug2_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(UAug2_BMSbin.dDOXY_Reg,UAug2_BMSbin.dTA_Reg,'o')
plot(UAug2_BMSbin.dDOXY_Reg,UAug2_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
title('NEC:NEP Ratio from Gradients');
str1 = num2str(UAug2_BMSbin.Ratio_G,2);
str2 = num2str(UAug2_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')

%% drag plots 

U_U02 = (UAug2_BMSbin.U0).^2
U_ustar2 = (UAug2_BMSbin.ustar).^2
% find the slope or regression coeficient:
U_mdl = fitlm(U_U02,U_ustar2)
U_reg =  (0.0075524*U_U02)-2.7e-07

close all
figure 
hold on 
scatter(U_U02,U_ustar2) 
plot (U_U02, U_reg)
xlabel('U0^2');
ylabel('ustar^2');
title('Cudjoe Aug 2 Drag') 
annotation('textbox', [0.132552083333333 0.83739406779661 0.150260416666667 0.0832627118644117],...
    'String',{'Slope:   0.0075524','Root Mean Squared Error: 2.15e-05','R^2: 0.643'},...
    'FitBoxToText','off')



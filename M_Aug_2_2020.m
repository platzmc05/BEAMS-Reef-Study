% Aug_2 Data Analysis from Marker 32 Reef
% Michelle Platz - USF 
% 3/10/2021

% SeapHOx sensor deployed 8/12/2020 at 9:15 am EST
    % SP datafile: 'M_812BMS.txt'
    % pump 1 height above benthos = 70 cm  
    % pump 2 height above benthos = 20 cm 
% ADP sensor deployed 8/05/2020 at 9:15 am EST
    % ADCP datafile: 'M_8_05'
    % height from substrate to ADCP head = 18 cm

close all
clc
clear all
%% Initial look at data
% ***** create MAug2_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
MAug2_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('M_812BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
MAug2_SPraw.DOXY = MAug2_SPraw.O2SATPER.*calcO2sat(MAug2_SPraw.MCAT_TC, MAug2_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
MAug2_SPraw.pHint_prelim = calc_dfet_pHint(MAug2_SPraw.Vint, MAug2_SPraw.DFET_TC, -0.4);

% ***** create MAug2_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
MAug2_SP = parse_to_pumpheights_ARM_2pump_Mar19(MAug2_SPraw, [datenum('08-12-2020 14:01:30'), datenum('08-15-2020 02:31:00')]);

% Calculate Gradients 
MAug2_SP = calc_TA_gradientV2(MAug2_SP, 2368.31, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

MAug2_SP.dDOXY = MAug2_SP.DOXY(1,:) - MAug2_SP.DOXY(2,:); %Oxygen Gradient
MAug2_SP.dpH = MAug2_SP.pH(1,:) - MAug2_SP.pH(2,:); %pH Gradient 
MAug2_SP.dTA = MAug2_SP.TAtop - MAug2_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
MAug2_DateString = {'08/12/2020 12:00:00';'08/13/2020 12:00:00';'08/14/2020 12:00:00';'08/15/2020 12:00:00';'08/16/2020 12:00:00';...
    '08/17/2020 12:00:00';'08/18/2020 12:00:00';'08/19/2020 12:00:00';'08/20/2020 12:00:00';'08/21/2020 12:00:00';'08/22/2020 12:00:00';...
    '08/23/2020 12:00:00';'08/24/2020 12:00:00';'08/25/2020 12:00:00';'08/26/2020 12:00:00';'08/27/2020 12:00:00';'08/28/2020 12:00:00';...
    '08/29/2020 12:00:00';'08/30/2020 12:00:00';'08/31/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
MAug2_tick = datenum(MAug2_DateString,formatIn);

MAug2_Xrange = [datenum('08-12-2020 14:00:00'), datenum('08-15-2020 02:31:00')];

figure
hold on; box on;
plot(MAug2_SP.SDN, MAug2_SP.dDOXY); %oxygen gradient 
plot(MAug2_SP.SDN, MAug2_SP.dTA); %TA gradient
plot(MAug2_SP.SDN, zeros(size(MAug2_SP.SDN))); %zero line
set(gca, 'xlim', MAug2_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 Aug 2020 Unbinned DO and TA Gradients');


%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
MAug2_SiteChar.SDN = MAug2_SPraw.SDN;
MAug2_SiteChar.TC = MAug2_SPraw.OPT_TC;
MAug2_SiteChar.PAR = MAug2_SPraw.PAR;
MAug2_SiteChar.PSAL = MAug2_SPraw.PSAL;
MAug2_SiteChar.Pres = MAug2_SPraw.Pres;

close all
figure 
hold on; 
plot(MAug2_SiteChar.SDN, MAug2_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
MAug2_SiteChar.SDN = MAug2_SPraw.SDN(25:end);
MAug2_SiteChar.TC = MAug2_SPraw.OPT_TC(25:end);
MAug2_SiteChar.PAR = MAug2_SPraw.PAR(25:end);
MAug2_SiteChar.PSAL = MAug2_SPraw.PSAL(25:end);
MAug2_SiteChar.Pres = MAug2_SPraw.Pres(25:end);

%extract full length of ADCP datafile  --> already done in Aug1 script
% % % MAug2_ADfull=aquadoppraw2mat('U_8_05', 70, [datenum('06-05-2020 08:00:00'), datenum('09-01-2020 08:00:00')]);
% % % 
% % % % add AD variables to Site Char 
% % % MAug2_SiteChar.AD_SDN = MAug2_ADfull.SDN;
% % % MAug2_SiteChar.AD_Pres = MAug2_ADfull.Pres;
% % % MAug2_SiteChar.AD_TC = MAug2_ADfull.TC;
% % % MAug2_SiteChar.bin_depth = MAug2_ADfull.bin_depth;
% % % MAug2_SiteChar.u = MAug2_ADfull.u;
% % % MAug2_SiteChar.v = MAug2_ADfull.v;
% % % MAug2_SiteChar.w = MAug2_ADfull.w;
% % % MAug2_SiteChar.uv = MAug2_ADfull.uv;
% % % MAug2_SiteChar.direction = MAug2_ADfull.direction;
% % % 
% % % %Plot to see when surface interval observations are
% % % close all
% % % figure 
% % % hold on; 
% % % plot(MAug2_SiteChar.AD_SDN, MAug2_SiteChar.AD_Pres)
% % % 
% % % %clip ends of data to remove surfave interval observations 
% % % MAug2_SiteChar.AD_SDN = MAug2_ADfull.SDN(95:end);
% % % MAug2_SiteChar.AD_Pres = MAug2_ADfull.Pres(95:end);
% % % MAug2_SiteChar.AD_TC = MAug2_ADfull.TC(95:end);
% % % MAug2_SiteChar.bin_depth = MAug2_ADfull.bin_depth;
% % % MAug2_SiteChar.u = MAug2_ADfull.u(:,95:end);
% % % MAug2_SiteChar.v = MAug2_ADfull.v(:,95:end);
% % % MAug2_SiteChar.w = MAug2_ADfull.w(:,95:end);
% % % MAug2_SiteChar.uv = MAug2_ADfull.uv(:,95:end);
% % % MAug2_SiteChar.direction = MAug2_ADfull.direction(:,95:end);
% % % 
% % % % find U0 
% % % MAug2_z1 = 0.70;
% % % MAug2_z2 = 0.20;
% % % MAug2_ADheight = 0.18;
% % % MAug2_ADbin_depth_1m = 1-(MAug2_ADheight);% = 0.82
% % % MAug2_i1m = find(MAug2_SiteChar.bin_depth==(0.82));
% % % MAug2_SiteChar.U0 = MAug2_SiteChar.uv(MAug2_i1m,:);

% save data in separate datastructure
save('MAug220_SiteChar_2.mat', 'MAug2_SiteChar' )
%% Constrain Xrange from graph results and extract good gradient data - 
close all 

MAug2_good_Xrange = [datenum('08-12-2020 14:00:00'), datenum('08-15-2020 02:31:00')];

% plot to check range is correct
figure
hold on; box on;
plot(MAug2_SP.SDN, MAug2_SP.dDOXY); %oxygen gradient 
plot(MAug2_SP.SDN, MAug2_SP.dTA); %TA gradient
plot(MAug2_SP.SDN, zeros(size(MAug2_SP.SDN))); %zero line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 Aug 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all

M_Aug_2_BMS_idx_start = find(MAug2_SP.SDN==datenum('08-12-2020 14:01:30'))
M_Aug_2_BMS_idx_end = find(MAug2_SP.SDN==datenum('08-15-2020 02:31:00'))

% Create new data vectors of just the good data
M_Aug_2_BMS_good_data = 3:241;
Initial_data_points = length(M_Aug_2_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(MAug2_SP);
for v = 1:length(vars)
    MAug2_SP.(vars{v}) = (MAug2_SP.(vars{v})(:,M_Aug_2_BMS_good_data));
end

datestr(MAug2_SP.SDN(1))
datestr(MAug2_SP.SDN(end))

%% *************** ADCP DATA ****************
% ***** create new data structure: MAug2_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
MAug2_AD=aquadoppraw2mat('M_8_05', 70, [datenum('12-Aug-2020 14:31:30'), datenum('15-Aug-2020 02:16:30')]);

%averages data to the middle of the minute interval spacified 
MAug2_ADavg = average_aquadopp(MAug2_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.7m (pump 1) and 0.2m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.52m (bin 42) and 0.02m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
MAug2_Sal_est = mean(MAug2_SP.PSAL(1,:))

[MAug2_ADavg] = ustar_from_aquadopp2(MAug2_ADavg,[0.52 0.11], MAug2_Sal_est); %bins adjusted 

clc
%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[MAug2_ADavg] = ustar_McGillis_Method(MAug2_ADavg, 0.70, 0.20, 42, 1);

%compare
close all
figure
hold on 
ustar_plot = plot(MAug2_ADavg.SDN, MAug2_ADavg.ustar, 'r');
ustar_WM_plot = plot(MAug2_ADavg.SDN, MAug2_ADavg.ustar_WM);
plot(MAug2_ADavg.SDN, zeros(size(MAug2_ADavg.SDN)),'k');
set(gca, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug2_ADavg.SDN(1) MAug2_ADavg.SDN(end)])
xlabel('Aug Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: MAug2_BMS *****

ADavg_vars = fieldnames(MAug2_ADavg);
for v = 1:length(ADavg_vars)
    MAug2_BMS.(ADavg_vars{v}) = (MAug2_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(MAug2_SP);
for v = 1:length(SP_vars)
    MAug2_BMS.(SP_vars{v}) = (MAug2_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(MAug2_BMS.SDN)
% min: 13-28-43-58 becuase SP was restarted in the field rather than on the
% minute


%% %% *************** Calculate Fluxes ****************

% actual pump heights  = 0.70m (pump 1) and 0.20m (pump 2) 
% 0.18m from substrate to ACDP head in Aug at U 
% adjusted height = 0.52 m (bin 42) and 0.02 m (too shallow) (bin 1)  
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[MAug2_BMS] = calc_NCC_3(MAug2_BMS,[0.52 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(MAug2_BMS.DOXY(1,:));
[MAug2_BMS] = calc_NCP_3(MAug2_BMS, [0.52 0.11],C1guess); %estimate C1 guess using median DOXY(1,:) value  

% Plot NCP and NCC
%close all
figure
hold on; box on; 
NEPplot = plot(MAug2_BMS.SDN, MAug2_BMS.NEP);
NECplot = plot(MAug2_BMS.SDN, MAug2_BMS.NEC);
plot(MAug2_BMS.SDN, zeros(size(MAug2_BMS.SDN)),'k');
set(gca, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug2_BMS.SDN(1) MAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');

%McGillis method flux calculations 
[MAug2_BMS] = calc_NCP_McGillis_Method(MAug2_BMS, 0.70, 0.20, MAug2_Sal_est);
[MAug2_BMS] = calc_NCC_McGillis_Method(MAug2_BMS, 0.70, 0.20, MAug2_Sal_est);

% Plot NCP and NCC
% close all
figure
hold on; box on; 
NEPplot = plot(MAug2_BMS.SDN, MAug2_BMS.NEP_WM);
NECplot = plot(MAug2_BMS.SDN, MAug2_BMS.NEC_WM);
plot(MAug2_BMS.SDN, zeros(size(MAug2_BMS.SDN)),'k');
set(gca, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug2_BMS.SDN(1) MAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 WM Fluxes');

%% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(MAug2_BMS.SDN, MAug2_BMS.NEP);
NEPplotWM = plot(MAug2_BMS.SDN, MAug2_BMS.NEP_WM);
plot(MAug2_BMS.SDN, zeros(size(MAug2_BMS.SDN)),'k');
set(gca, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug2_BMS.SDN(1) MAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(MAug2_BMS.SDN, MAug2_BMS.NEC);
NECplotWM = plot(MAug2_BMS.SDN, MAug2_BMS.NEC_WM);
plot(MAug2_BMS.SDN, zeros(size(MAug2_BMS.SDN)),'k');
set(gca, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug2_BMS.SDN(1) MAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
MAug2_BMS.DOXYstd = std(MAug2_BMS.DOXY);

% get DOXY std
MAug2_idoxystd = find(MAug2_BMS.DOXYstd > 2);% 58.4% of data is greater than 0.8 stdev
MAug2_ihighdoxystd2 = [];
for i = 1:length(MAug2_idoxystd)
    
    MAug2_ihighdoxystd2 = vertcat(MAug2_ihighdoxystd2,[MAug2_idoxystd(i)-1:1:MAug2_idoxystd(i)+1]');
end
% get unique IDs
MAug2_ihighdoxystd2 = unique(MAug2_ihighdoxystd2);
% remove 0's and out of index values
MAug2_ihighdoxystd2(MAug2_ihighdoxystd2==0) = [];
MAug2_ihighdoxystd2(MAug2_ihighdoxystd2> length(MAug2_BMS.SDN)) = [];

% make it into index
trex = false(size(MAug2_BMS.SDN));
trex(MAug2_ihighdoxystd2) = true;
MAug2_ihighdoxystd2 = trex;
clear trex;

MAug2_BMS.NEP_QC = MAug2_BMS.NEP;
MAug2_BMS.NEC_QC = MAug2_BMS.NEC;
MAug2_BMS.dDOXY_QC = MAug2_BMS.dDOXY; %DO gradient
MAug2_BMS.dTA_QC = MAug2_BMS.dTA;     %TA gradient
MAug2_BMS.NEP_WM_QC = MAug2_BMS.NEP_WM;
MAug2_BMS.NEC_WM_QC = MAug2_BMS.NEC_WM;

% set observations when DOXYstd>0.8 to NaN
MAug2_BMS.NEP_QC(MAug2_ihighdoxystd2) = NaN;
MAug2_BMS.NEC_QC(:,MAug2_ihighdoxystd2) = NaN;
MAug2_BMS.dDOXY_QC(MAug2_ihighdoxystd2) = NaN; %DO gradient
MAug2_BMS.dTA_QC(:,MAug2_ihighdoxystd2) = NaN;   %TA gradient
MAug2_BMS.NEP_WM_QC(MAug2_ihighdoxystd2) = NaN;
MAug2_BMS.NEC_WM_QC(:,MAug2_ihighdoxystd2) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(MAug2_BMS.SDN, MAug2_BMS.NEP, 'k');
NEPplotQC = plot(MAug2_BMS.SDN, MAug2_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(MAug2_BMS.SDN, zeros(size(MAug2_BMS.SDN)),'k');
set(gca, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug2_BMS.SDN(1) MAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');


%WM Plot
%close all
figure
hold on; box on;
NEPplot = plot(MAug2_BMS.SDN, MAug2_BMS.NEP_WM, 'k');
NEPplotQC = plot(MAug2_BMS.SDN, MAug2_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
plot(MAug2_BMS.SDN, zeros(size(MAug2_BMS.SDN)),'k');
set(gca, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MAug2_BMS.SDN(1) MAug2_BMS.SDN(end)])
xlabel('Aug Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Fluxes');


%% Bin data to hourly intervals

X = floor(nanmin(MAug2_BMS.SDN)):1/24:ceil(nanmax(MAug2_BMS.SDN));

MAug2_BMSbin.SDN = X;
% variables from 2 differnet heights
MAug2_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.DOXY(1,:), X);
MAug2_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.DOXY(2,:), X);

MAug2_BMSbin.pH(1,:)       = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.pH(1,:), X);
MAug2_BMSbin.pH(2,:)       = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.pH(2,:), X);

MAug2_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.PSAL(1,:), X);
MAug2_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.PSAL(2,:), X);

MAug2_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.O2SATPER(1,:), X);
MAug2_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.O2SATPER(2,:), X);

MAug2_BMSbin.Pres(1,:)     = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.Pres(1,:), X);
MAug2_BMSbin.Pres(2,:)     = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.Pres(2,:), X);

MAug2_BMSbin.DENS(1,:)     = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.DENS(1,:), X);
MAug2_BMSbin.DENS(2,:)     = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.DENS(2,:), X);

MAug2_BMSbin.PAR(1,:)      = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.PAR(1,:), X);
MAug2_BMSbin.PAR(2,:)      = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.PAR(2,:), X);

MAug2_BMSbin.bin_depth     = MAug2_BMS.bin_depth;

for i = 1:108
    MAug2_BMSbin.uv(i,:)   = bin_data_to_X_GF(MAug2_BMS.SDN,MAug2_BMS.uv(i,:), X);
end
% bin data hourly. Vector variables 
MAug2_BMSbin.PRES  = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.Pres, X);
MAug2_BMSbin.U0       = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.U0, X);
MAug2_BMSbin.DIR      = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.direction, X);
MAug2_BMSbin.ustar    = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.ustar, X);
MAug2_BMSbin.ustar_rm = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.ustar_runmean, X);
MAug2_BMSbin.dTA      = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.dTA, X);
MAug2_BMSbin.dTA_QC   = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.dTA_QC, X);
MAug2_BMSbin.dpH      = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.dpH, X);
MAug2_BMSbin.dDOXY    = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.dDOXY, X);
MAug2_BMSbin.dDOXY_QC = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.dDOXY_QC, X);
MAug2_BMSbin.NEP      = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.NEP, X);
MAug2_BMSbin.NEP_QC   = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.NEP_QC, X);
MAug2_BMSbin.NEC      = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.NEC, X);
MAug2_BMSbin.NEC_QC   = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.NEC_QC, X);
MAug2_BMSbin.NEP_WM   = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.NEP_WM, X);
MAug2_BMSbin.NEP_WM_QC= bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.NEP_WM_QC, X);
MAug2_BMSbin.NEC_WM   = bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.NEC_WM, X);
MAug2_BMSbin.NEC_WM_QC= bin_data_to_X_GF(MAug2_BMS.SDN, MAug2_BMS.NEC_WM_QC, X);

%% Plot Binned Fluxes 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP, 'b'); 
% NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC, 'r-'); 
% plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 Hourly Binned Fluxes FluxFit');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_QC, 'r-'); 
plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Hourly Binned Fluxes FluxFit QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_WM, 'r-'); 
% plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 WM Original Hourly Binned Fluxes');

% close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_WM_QC, 'r-'); 
% plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 WM Full QC Hourly Binned Fluxes');

%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
DOplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.dDOXY, 'c'); 
TAplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.dTA, 'k-'); 
plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug2 Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 Aug2 2020 Binned Gradiets');



%% ***QC*** Remove sections when velocity is too slow
ibad = MAug2_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow

% MAug2_BMSbin.NEP(ibad) = NaN;
% MAug2_BMSbin.NEC(ibad) = NaN;
MAug2_BMSbin.NEP_QC(ibad) = NaN;
MAug2_BMSbin.NEC_QC(ibad) = NaN;
MAug2_BMSbin.dDOXY_QC(ibad) = NaN;
MAug2_BMSbin.dTA_QC(ibad) = NaN;

% MAug2_BMSbin.NEP_WM(ibad) = NaN;
% MAug2_BMSbin.NEC_WM(ibad) = NaN;
MAug2_BMSbin.NEP_WM_QC(ibad) = NaN;
MAug2_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP, 'b'); 
% NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC, 'r-'); 
% plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_QC, 'r-'); 
plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_WM, 'r-'); 
% plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Aug Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 Aug 2020 Hourly Binned WM Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_WM_QC, 'r-'); 
plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Aug 2020 Hourly Binned WM Fluxes Full QC');

%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
figure
hold on; box on;
boxplot(MAug2_BMSbin.NEP_QC)
ylabel('NEP')
title('Aug2 NEP Boxplots')

figure
hold on; box on;
boxplot(MAug2_BMSbin.NEC_QC)
ylabel('NEC')
title('Aug2 NEC Boxplots')

figure
hold on; box on;
boxplot(MAug2_BMSbin.dDOXY_QC)
ylabel('DO')
title('Aug2 dDO Boxplots')

figure
hold on; box on;
boxplot(MAug2_BMSbin.dTA_QC)
ylabel('TA')
title('Aug2 dTA Boxplots')

% Plot Profiles at outliers - 
close all 
for i = 253   %1:length(MAug2_BMSbin.SDN)
    figure (i)
    scatter(MAug2_BMSbin.uv(1:42,i), MAug2_BMSbin.bin_depth(1:42))
    title(['Cudjoe Aug2 Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end
%outliers to be removed: 
% MAug2_BMSbin.NEP_QC() = NaN;
% MAug2_BMSbin.NEC_QC() = NaN;
% MAug2_BMSbin.dDOXY_QC() = NaN;
% MAug2_BMSbin.dTA_QC() = NaN;

close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_QC, 'r-'); 
plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned Fluxes Full QC');

%% Plot Diel Curves

% nepdbin = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP, 24);
% necdbin = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC, 24);
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

nepdbin_QC = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_QC, 24);
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
% nepdbin_WM = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_WM, 24);
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

nepdbin_WM_QC = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_WM_QC, 24);
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
% Extract daytime data using MAug2_BMSbin.PAR
MAug2_inight = MAug2_BMSbin.PAR(1,:) < 1; %find all nightime datapoints 

%create new arrays for daytime data
MAug2_BMSbin.SDN_day = MAug2_BMSbin.SDN;
MAug2_BMSbin.PAR_day = MAug2_BMSbin.PAR(1,:);

MAug2_BMSbin.NEP_day = MAug2_BMSbin.NEP;
MAug2_BMSbin.NEC_day = MAug2_BMSbin.NEC;
MAug2_BMSbin.NEP_day_QC = MAug2_BMSbin.NEP_QC;
MAug2_BMSbin.NEC_day_QC = MAug2_BMSbin.NEC_QC;

MAug2_BMSbin.dDOXY_day_QC = MAug2_BMSbin.dDOXY_QC;      %DO Gradient
MAug2_BMSbin.dTA_day_QC = MAug2_BMSbin.dTA_QC;          %TA Gradient

MAug2_BMSbin.NEP_WM_day = MAug2_BMSbin.NEP_WM;
MAug2_BMSbin.NEC_WM_day = MAug2_BMSbin.NEC_WM;
MAug2_BMSbin.NEP_WM_day_QC = MAug2_BMSbin.NEP_WM_QC;
MAug2_BMSbin.NEC_WM_day_QC = MAug2_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
MAug2_BMSbin.SDN_day(MAug2_inight) = NaN;
MAug2_BMSbin.PAR_day (MAug2_inight) = NaN;

MAug2_BMSbin.NEP_day(MAug2_inight) = NaN;
MAug2_BMSbin.NEC_day(MAug2_inight) = NaN;
MAug2_BMSbin.NEP_day_QC(MAug2_inight) = NaN;
MAug2_BMSbin.NEC_day_QC(MAug2_inight) = NaN;

MAug2_BMSbin.dDOXY_day_QC(MAug2_inight) = NaN;      %DO Gradient
MAug2_BMSbin.dTA_day_QC(MAug2_inight) = NaN;        %TA Gradient

MAug2_BMSbin.NEP_WM_day(MAug2_inight) = NaN;
MAug2_BMSbin.NEC_WM_day(MAug2_inight) = NaN;
MAug2_BMSbin.NEP_WM_day_QC(MAug2_inight) = NaN;
MAug2_BMSbin.NEC_WM_day_QC(MAug2_inight) = NaN;
%Plot to check only nighttime points removed
figure 
hold on
scatter(MAug2_BMSbin.SDN, MAug2_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MAug2_BMSbin.SDN_day, MAug2_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
MAug2_BMSbin.NEP_day(isnan(MAug2_BMSbin.NEP_day))=[];
MAug2_BMSbin.NEC_day(isnan(MAug2_BMSbin.NEC_day))=[];
MAug2_BMSbin.NEP_day_QC(isnan(MAug2_BMSbin.NEP_day_QC))=[];
MAug2_BMSbin.NEC_day_QC(isnan(MAug2_BMSbin.NEC_day_QC))=[];

MAug2_BMSbin.dDOXY_day_QC(isnan(MAug2_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
MAug2_BMSbin.dTA_day_QC(isnan(MAug2_BMSbin.dTA_day_QC))=[];       %TA Gradient

MAug2_BMSbin.NEP_WM_day(isnan(MAug2_BMSbin.NEP_WM_day))=[];
MAug2_BMSbin.NEC_WM_day(isnan(MAug2_BMSbin.NEC_WM_day))=[];
MAug2_BMSbin.NEP_WM_day_QC(isnan(MAug2_BMSbin.NEP_WM_day_QC))=[];
MAug2_BMSbin.NEC_WM_day_QC(isnan(MAug2_BMSbin.NEC_WM_day_QC))=[];

% create nighttime hours datasets
MAug2_BMSbin.SDN_night = MAug2_BMSbin.SDN;
MAug2_BMSbin.PAR_night = MAug2_BMSbin.PAR(1,:);

MAug2_BMSbin.NEP_night = MAug2_BMSbin.NEP;
MAug2_BMSbin.NEC_night = MAug2_BMSbin.NEC;
MAug2_BMSbin.NEP_night_QC = MAug2_BMSbin.NEP_QC;
MAug2_BMSbin.NEC_night_QC = MAug2_BMSbin.NEC_QC;

MAug2_BMSbin.dDOXY_night_QC = MAug2_BMSbin.dDOXY_QC;      %DO Gradient
MAug2_BMSbin.dTA_night_QC = MAug2_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
MAug2_BMSbin.SDN_night=MAug2_BMSbin.SDN_night(MAug2_inight);
MAug2_BMSbin.PAR_night=MAug2_BMSbin.PAR_night(MAug2_inight);

MAug2_BMSbin.NEP_night_QC=MAug2_BMSbin.NEP_night_QC(MAug2_inight);
MAug2_BMSbin.NEC_night_QC=MAug2_BMSbin.NEC_night_QC(MAug2_inight);

MAug2_BMSbin.dDOXY_night_QC=MAug2_BMSbin.dDOXY_night_QC(MAug2_inight);      %DO Gradient
MAug2_BMSbin.dTA_night_QC=MAug2_BMSbin.dTA_night_QC(MAug2_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(MAug2_BMSbin.SDN, MAug2_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MAug2_BMSbin.SDN_night, MAug2_BMSbin.PAR_night, 'r.'); % night plot

%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(MAug2_BMSbin.NEP_day,MAug2_BMSbin.NEC_day);
% MAug2_BMSbin.Reg_Line = m*MAug2_BMSbin.NEP_day + b;
% MAug2_BMSbin.Ratio = m;
% MAug2_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(MAug2_BMSbin.NEP_day,MAug2_BMSbin.NEC_day,'o')
% plot(MAug2_BMSbin.NEP_day,MAug2_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Aug 2020 Post-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug2_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug2_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(MAug2_BMSbin.NEP_day_QC,MAug2_BMSbin.NEC_day_QC);
MAug2_BMSbin.Reg_Line_QC = m_QC*MAug2_BMSbin.NEP_day_QC + b_QC;
MAug2_BMSbin.Ratio_QC = m_QC;
MAug2_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(MAug2_BMSbin.NEP_day_QC,MAug2_BMSbin.NEC_day_QC,'o')
plot(MAug2_BMSbin.NEP_day_QC,MAug2_BMSbin.Reg_Line_QC,'r')
%ylim([-50 50])
%xlim([-25 25])
xlabel('NCP');
ylabel('NCC');
title('Marker 32 Aug 2020 Post-Restoration NCC:NCP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug2_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug2_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(MAug2_BMSbin.NEP_WM_day,MAug2_BMSbin.NEC_WM_day);
% MAug2_BMSbin.Reg_Line_WM = m_WM*MAug2_BMSbin.NEP_WM_day + b_WM;
% MAug2_BMSbin.Ratio_WM = m_WM;
% MAug2_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(MAug2_BMSbin.NEP_WM_day,MAug2_BMSbin.NEC_WM_day,'o')
% plot(MAug2_BMSbin.NEP_WM_day,MAug2_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Aug 2020 Post-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug2_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug2_BMSbin.R2_WM)


[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(MAug2_BMSbin.NEP_WM_day_QC,MAug2_BMSbin.NEC_WM_day_QC);
MAug2_BMSbin.Reg_Line_WM_QC = m_WM_QC*MAug2_BMSbin.NEP_WM_day_QC + b_WM_QC;
MAug2_BMSbin.Ratio_WM_QC = m_WM_QC;
MAug2_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(MAug2_BMSbin.NEP_WM_day_QC,MAug2_BMSbin.NEC_WM_day_QC,'o')
% plot(MAug2_BMSbin.NEP_WM_day_QC,MAug2_BMSbin.Reg_Line_WM_QC,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 Aug 2020 Post-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug2_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug2_BMSbin.R2_WM_QC)


%% For NEC:NEP Regressions Using Gradients
close all 
clc

% multiply o2 gradient by -1 for O2 production
MAug2_BMSbin.dDOXY_Reg = -1.*MAug2_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
MAug2_BMSbin.dTA_Reg = 0.5.*MAug2_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(MAug2_BMSbin.dDOXY_day_QC);

figure 
hold on 
DOday = plot(1:Xlength, MAug2_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, MAug2_BMSbin.dDOXY_Reg);
xlabel('Aug2 Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Marker 32 Aug2 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, MAug2_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, MAug2_BMSbin.dTA_Reg);
xlabel('Aug2 Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Marker 32 Aug2 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(MAug2_BMSbin.dDOXY_Reg, MAug2_BMSbin.dTA_Reg);
MAug2_BMSbin.Reg_Line_G = m_G*MAug2_BMSbin.dDOXY_Reg + b_G;
MAug2_BMSbin.Ratio_G = m_G;
MAug2_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(MAug2_BMSbin.dDOXY_Reg,MAug2_BMSbin.dTA_Reg,'o')
plot(MAug2_BMSbin.dDOXY_Reg,MAug2_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Marker 32 Aug2 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MAug2_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MAug2_BMSbin.R2_G)

clc
disp('Finished with Aug 2 Metabolism Calculations');

save('MAug220_2.mat', 'MAug2_BMS', 'MAug2_BMSbin');

%% Plot Profiles 
% plot for profile within pump heights
% close all 
% for i =1:100 %length(MAug2_ADavg.SDN)
%     figure (i)
%     scatter(MAug2_ADavg.uv(1:108,i), MAug2_ADavg.bin_depth(1:108))
%     title(['Marker 32 Aug Velocity Profile Number ',num2str(i),])
%     xlabel('Velocity (m/s)');
%     ylabel('Height (m)');
% end

%% Subplots 
close all
clc

sgtitle('Marker 32 Aug2 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'08/12 12:00';'08/13 12:00';'08/14 12:00';'08/15 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
% legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'southwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(MAug2_BMSbin.SDN, zeros(size(MAug2_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MAug2_good_Xrange, 'XTick', MAug2_tick, 'xticklabel', MAug2_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Aug2 Days');
xticklabels({'08/12 12:00';'08/13 12:00';'08/14 12:00';'08/15 12:00'})
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MAug2_BMSbin.SDN, MAug2_BMSbin.NEC_QC, 24);
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
plot(MAug2_BMSbin.NEP_day_QC,MAug2_BMSbin.NEC_day_QC,'o')
plot(MAug2_BMSbin.NEP_day_QC,MAug2_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(MAug2_BMSbin.Ratio_QC,2);
str2 = num2str(MAug2_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(MAug2_BMSbin.dDOXY_Reg,MAug2_BMSbin.dTA_Reg,'o')
plot(MAug2_BMSbin.dDOXY_Reg,MAug2_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
title('dTA:dDO Ratio from Gradients');
str1 = num2str(MAug2_BMSbin.Ratio_G,2);
str2 = num2str(MAug2_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.075, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.075, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')




% Sept Data Analysis from Marker 32 Reef
% Michelle Platz - USF 
% 3/9/2021

% SeapHOx sensor deployed 9/2/2020 
    % SP datafile: 'M_902BMS.txt'
    % pump 1 height above benthos = 70 cm  
    % pump 2 height above benthos = 20 cm 
% ADP sensor deployed 9/02/2020 
    % ADCP datafile: 'M_9_02'
    % height from substrate to ADCP head = 18 cm

close all
clc
clear all
%% Initial look at data
% ***** create MSept_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
MSept_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('M_902BMS1.txt');

%calculate O2 saturation concentration using temperature and salinity
MSept_SPraw.DOXY = MSept_SPraw.O2SATPER.*calcO2sat(MSept_SPraw.MCAT_TC, MSept_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
MSept_SPraw.pHint_prelim = calc_dfet_pHint(MSept_SPraw.Vint, MSept_SPraw.DFET_TC, -0.4);

% ***** create MSept_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
MSept_SP = parse_to_pumpheights_ARM_2pump_Mar19(MSept_SPraw, [datenum('09-02-2020 09:00:00'), datenum('09-28-2020 13:29:30')]);

% Calculate Gradients 
MSept_SP = calc_TA_gradientV2(MSept_SP, 2368.31, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

MSept_SP.dDOXY = MSept_SP.DOXY(1,:) - MSept_SP.DOXY(2,:); %Oxygen Gradient
MSept_SP.dpH = MSept_SP.pH(1,:) - MSept_SP.pH(2,:); %pH Gradient 
MSept_SP.dTA = MSept_SP.TAtop - MSept_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
MSept_DateString = {'09/03/2020 12:00:00';'09/04/2020 12:00:00';'09/05/2020 12:00:00';...
    '09/06/2020 12:00:00';'09/07/2020 12:00:00';'09/08/2020 12:00:00';'09/09/2020 12:00:00';'09/10/2020 12:00:00';'09/11/2020 12:00:00';'09/12/2020 12:00:00';...
    '09/13/2020 12:00:00';'09/14/2020 12:00:00';'09/15/2020 12:00:00';'09/16/2020 12:00:00';'09/17/2020 12:00:00';'09/18/2020 12:00:00';'09/19/2020 12:00:00';...
    '09/20/2020 12:00:00';'09/21/2020 12:00:00';'09/22/2020 12:00:00';'09/23/2020 12:00:00';'09/24/2020 12:00:00';'09/25/2020 12:00:00';'09/26/2020 12:00:00';...
    '09/27/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
MSept_tick = datenum(MSept_DateString,formatIn);

MSept_Xrange = [datenum('09-02-2020 09:00:00'), datenum('09-28-2020 13:29:30')];

figure
hold on; box on;
plot(MSept_SP.SDN, MSept_SP.dDOXY); %oxygen gradient 
plot(MSept_SP.SDN, MSept_SP.dTA); %TA gradient
plot(MSept_SP.SDN, zeros(size(MSept_SP.SDN))); %zero line
set(gca, 'xlim', MSept_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32Sept 2020 Unbinned DO and TA Gradients');


%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
MSept_SiteChar.SDN = MSept_SPraw.SDN;
MSept_SiteChar.TC = MSept_SPraw.OPT_TC;
MSept_SiteChar.PAR = MSept_SPraw.PAR;
MSept_SiteChar.PSAL = MSept_SPraw.PSAL;
MSept_SiteChar.Pres = MSept_SPraw.Pres;

%Plot pressure data to see when surface interval observations are
close all
figure 
hold on; 
plot(MSept_SiteChar.SDN, MSept_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
MSept_SiteChar.SDN = MSept_SPraw.SDN(121:74960);
MSept_SiteChar.TC = MSept_SPraw.OPT_TC(121:74960);
MSept_SiteChar.PAR = MSept_SPraw.PAR(121:74960);
MSept_SiteChar.PSAL = MSept_SPraw.PSAL(121:74960);
MSept_SiteChar.Pres = MSept_SPraw.Pres(121:74960);

%extract full length of ADCP datafile  
MSept_ADfull=aquadoppraw2mat('M_9_02', 70, [datenum('09-02-2020 16:00:00'), datenum('11-14-2020 13:49:25')]);

% add AD variables to Site Char 
MSept_SiteChar.AD_SDN = MSept_ADfull.SDN;
MSept_SiteChar.AD_Pres = MSept_ADfull.Pres;
MSept_SiteChar.AD_TC = MSept_ADfull.TC;
MSept_SiteChar.bin_depth = MSept_ADfull.bin_depth;
MSept_SiteChar.u = MSept_ADfull.u;
MSept_SiteChar.v = MSept_ADfull.v;
MSept_SiteChar.w = MSept_ADfull.w;
MSept_SiteChar.uv = MSept_ADfull.uv;
MSept_SiteChar.direction = MSept_ADfull.direction;

%Plot pressure data to see when surface interval observations are
close all
figure 
hold on; 
plot(MSept_SiteChar.AD_SDN, MSept_SiteChar.AD_Pres)

%clip ends of data to remove surfave interval observations 
% MSept_SiteChar.AD_SDN = MSept_ADfull.SDN(71:end);
% MSept_SiteChar.AD_Pres = MSept_ADfull.Pres(71:end);
% MSept_SiteChar.AD_TC = MSept_ADfull.TC(71:end);
% MSept_SiteChar.bin_depth = MSept_ADfull.bin_depth;
% MSept_SiteChar.u = MSept_ADfull.u(:,71:end);
% MSept_SiteChar.v = MSept_ADfull.v(:,71:end);
% MSept_SiteChar.w = MSept_ADfull.w(:,71:end);
% MSept_SiteChar.uv = MSept_ADfull.uv(:,71:end);
% MSept_SiteChar.direction = MSept_ADfull.direction(:,71:end);

% find U0 
MSept_z1 = 0.70;
MSept_z2 = 0.20;
MSept_ADheight = 0.18;
MSept_ADbin_depth_1m = 1-(MSept_ADheight);% = 0.82
MSept_i1m = find(MSept_SiteChar.bin_depth==(0.82));
MSept_SiteChar.U0 = MSept_SiteChar.uv(MSept_i1m,:);

% save data in separate datastructure
save('MSept20_SiteChar_2.mat', 'MSept_SiteChar' )


%% Constrain Xrange from graph results and extract good gradient data - 
close all 

MSept_good_Xrange = [datenum('09-02-2020 16:00:00'), datenum('09-11-2020 12:00:00')]; 

% plot to check range is correct
figure
hold on; box on;
plot(MSept_SP.SDN, MSept_SP.dDOXY); %oxygen gradient 
plot(MSept_SP.SDN, MSept_SP.dTA); %TA gradient
plot(MSept_SP.SDN, zeros(size(MSept_SP.SDN))); %zero line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32Sept 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all
M_Sept_BMS_idx_start = find(MSept_SP.SDN==datenum('09-02-2020 16:00:00'))
M_Sept_BMS_idx_end = find(MSept_SP.SDN==datenum('09-11-2020 12:00:00'))

% Create new data vectors of just the good data
M_Sept_BMS_good_data = M_Sept_BMS_idx_start:M_Sept_BMS_idx_end;
Initial_data_points = length(M_Sept_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(MSept_SP);
for v = 1:length(vars)
    MSept_SP.(vars{v}) = (MSept_SP.(vars{v})(:,M_Sept_BMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: MSept_AD *****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
MSept_AD=aquadoppraw2mat('M_9_02', 70, [datenum('09-02-2020 16:00:00'), datenum('09-11-2020 12:15:00')]);

%averages data to the middle of the minute interval spacified 
MSept_ADavg = average_aquadopp(MSept_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.7m (pump 1) and 0.2m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.52m (bin 42) and 0.02m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
MSept_Sal_est = mean(MSept_SP.PSAL(1,3:end));

[MSept_ADavg] = ustar_from_aquadopp2(MSept_ADavg,[0.52 0.11], MSept_Sal_est); %bins adjusted 

clc
%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[MSept_ADavg] = ustar_McGillis_Method(MSept_ADavg, 0.70, 0.20, 42, 1);

%compare ustar calculation methods
close all
figure
hold on 
ustar_plot = plot(MSept_ADavg.SDN, MSept_ADavg.ustar, 'r');
ustar_WM_plot = plot(MSept_ADavg.SDN, MSept_ADavg.ustar_WM);
plot(MSept_ADavg.SDN, zeros(size(MSept_ADavg.SDN)),'k');
set(gca, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MSept_ADavg.SDN(1) MSept_ADavg.SDN(end)])
xlabel('Sept Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Marker 32Sept 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: MSept_BMS *****

ADavg_vars = fieldnames(MSept_ADavg);
for v = 1:length(ADavg_vars)
    MSept_BMS.(ADavg_vars{v}) = (MSept_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
SP_vars = fieldnames(MSept_SP);
for v = 1:length(SP_vars)
    MSept_BMS.(SP_vars{v}) = (MSept_SP.(SP_vars{v}));
end
% check that SDN is on 15 min interval
datestr(MSept_BMS.SDN)
% min: 13-28-43-58 becuase SP was restarted in the field rather than on the
% minute
%% %% *************** Calculate Fluxes ****************

% actual pump heights  = 0.70m (pump 1) and 0.20m (pump 2) 
% 0.18m from substrate to ACDP head in Sept at U 
% adjusted height = 0.52 m (bin 42) and 0.02 m (too shallow) (bin 1)  
clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[MSept_BMS] = calc_NCC_3(MSept_BMS,[0.52 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(MSept_BMS.DOXY(1,:));
[MSept_BMS] = calc_NCP_3(MSept_BMS, [0.52 0.11],C1guess); %estimate C1 guess using median DOXY(1,:) value  

% Plot NCP and NCC
%close all
figure
hold on; box on; 
NEPplot = plot(MSept_BMS.SDN, MSept_BMS.NEP);
NECplot = plot(MSept_BMS.SDN, MSept_BMS.NEC);
plot(MSept_BMS.SDN, zeros(size(MSept_BMS.SDN)),'k');
set(gca, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MSept_BMS.SDN(1) MSept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32Sept 2020 Fluxes');

%McGillis method flux calculations 
[MSept_BMS] = calc_NCP_McGillis_Method(MSept_BMS, 0.70, 0.20, MSept_Sal_est);
[MSept_BMS] = calc_NCC_McGillis_Method(MSept_BMS, 0.70, 0.20, MSept_Sal_est);

% Plot NCP and NCC
% close all
figure
hold on; box on; 
NEPplot = plot(MSept_BMS.SDN, MSept_BMS.NEP_WM);
NECplot = plot(MSept_BMS.SDN, MSept_BMS.NEC_WM);
plot(MSept_BMS.SDN, zeros(size(MSept_BMS.SDN)),'k');
set(gca, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MSept_BMS.SDN(1) MSept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32Sept 2020 WM Fluxes');

%% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(MSept_BMS.SDN, MSept_BMS.NEP);
NEPplotWM = plot(MSept_BMS.SDN, MSept_BMS.NEP_WM);
plot(MSept_BMS.SDN, zeros(size(MSept_BMS.SDN)),'k');
set(gca, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MSept_BMS.SDN(1) MSept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('Marker 32Sept 2020 Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(MSept_BMS.SDN, MSept_BMS.NEC);
NECplotWM = plot(MSept_BMS.SDN, MSept_BMS.NEC_WM);
plot(MSept_BMS.SDN, zeros(size(MSept_BMS.SDN)),'k');
set(gca, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MSept_BMS.SDN(1) MSept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('Marker 32Sept 2020 Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
MSept_BMS.DOXYstd = std(MSept_BMS.DOXY);

% get DOXY std
MSept_idoxystd = find(MSept_BMS.DOXYstd > 2);% 58.4% of data is greater than 0.8 stdev
MSept_ihighdoxystd = [];
for i = 1:length(MSept_idoxystd)
    
    MSept_ihighdoxystd = vertcat(MSept_ihighdoxystd,[MSept_idoxystd(i)-1:1:MSept_idoxystd(i)+1]');
end
% get unique IDs
MSept_ihighdoxystd = unique(MSept_ihighdoxystd);
% remove 0's and out of index values
MSept_ihighdoxystd(MSept_ihighdoxystd==0) = [];
MSept_ihighdoxystd(MSept_ihighdoxystd> length(MSept_BMS.SDN)) = [];

% make it into index
trex = false(size(MSept_BMS.SDN));
trex(MSept_ihighdoxystd) = true;
MSept_ihighdoxystd = trex;
clear trex;

%create _QC datasets to preserve original data
MSept_BMS.NEP_QC = MSept_BMS.NEP;
MSept_BMS.NEC_QC = MSept_BMS.NEC;
MSept_BMS.dDOXY_QC = MSept_BMS.dDOXY; %DO gradient
MSept_BMS.dTA_QC = MSept_BMS.dTA;     %TA gradient
MSept_BMS.NEP_WM_QC = MSept_BMS.NEP_WM;
MSept_BMS.NEC_WM_QC = MSept_BMS.NEC_WM;

% set observations when DOXYstd>X to NaN - only remove NaNs from QC datasets
MSept_BMS.NEP_QC(MSept_ihighdoxystd) = NaN;
MSept_BMS.NEC_QC(:,MSept_ihighdoxystd) = NaN;
MSept_BMS.dDOXY_QC(MSept_ihighdoxystd) = NaN; %DO gradient
MSept_BMS.dTA_QC(:,MSept_ihighdoxystd) = NaN;   %TA gradient
MSept_BMS.NEP_WM_QC(MSept_ihighdoxystd) = NaN;
MSept_BMS.NEC_WM_QC(:,MSept_ihighdoxystd) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(MSept_BMS.SDN, MSept_BMS.NEP, 'k');
NEPplotQC = plot(MSept_BMS.SDN, MSept_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(MSept_BMS.SDN, zeros(size(MSept_BMS.SDN)),'k');
set(gca, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MSept_BMS.SDN(1) MSept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32Sept 2020 Fluxes');


%WM Plot
%close all
figure
hold on; box on;
NEPplot = plot(MSept_BMS.SDN, MSept_BMS.NEP_WM, 'k');
NEPplotQC = plot(MSept_BMS.SDN, MSept_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
plot(MSept_BMS.SDN, zeros(size(MSept_BMS.SDN)),'k');
set(gca, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MSept_BMS.SDN(1) MSept_BMS.SDN(end)])
xlabel('Sept Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32Sept 2020 Fluxes');


%% Bin data to hourly intervals

X = floor(nanmin(MSept_BMS.SDN)):1/24:ceil(nanmax(MSept_BMS.SDN));

MSept_BMSbin.SDN = X;
% variables from 2 differnet heights
MSept_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.DOXY(1,:), X);
MSept_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.DOXY(2,:), X);

MSept_BMSbin.pH(1,:)       = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.pH(1,:), X);
MSept_BMSbin.pH(2,:)       = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.pH(2,:), X);

MSept_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.PSAL(1,:), X);
MSept_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.PSAL(2,:), X);

MSept_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.O2SATPER(1,:), X);
MSept_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.O2SATPER(2,:), X);

MSept_BMSbin.Pres(1,:)     = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.Pres(1,:), X);
MSept_BMSbin.Pres(2,:)     = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.Pres(2,:), X);

MSept_BMSbin.DENS(1,:)     = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.DENS(1,:), X);
MSept_BMSbin.DENS(2,:)     = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.DENS(2,:), X);

MSept_BMSbin.PAR(1,:)      = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.PAR(1,:), X);
MSept_BMSbin.PAR(2,:)      = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.PAR(2,:), X);

MSept_BMSbin.bin_depth     = MSept_BMS.bin_depth;

for i = 1:108
    MSept_BMSbin.uv(i,:)   = bin_data_to_X_GF(MSept_BMS.SDN,MSept_BMS.uv(i,:), X);
end

% bin data hourly. Vector variables 
MSept_BMSbin.PRES  = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.Pres, X);
MSept_BMSbin.U0       = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.U0, X);
MSept_BMSbin.DIR      = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.direction, X);
MSept_BMSbin.ustar    = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.ustar, X);
MSept_BMSbin.ustar_rm = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.ustar_runmean, X);
MSept_BMSbin.dTA      = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.dTA, X);
MSept_BMSbin.dTA_QC   = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.dTA_QC, X);
MSept_BMSbin.dpH      = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.dpH, X);
MSept_BMSbin.dDOXY    = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.dDOXY, X);
MSept_BMSbin.dDOXY_QC = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.dDOXY_QC, X);
MSept_BMSbin.NEP      = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.NEP, X);
MSept_BMSbin.NEP_QC   = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.NEP_QC, X);
MSept_BMSbin.NEC      = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.NEC, X);
MSept_BMSbin.NEC_QC   = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.NEC_QC, X);
MSept_BMSbin.NEP_WM   = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.NEP_WM, X);
MSept_BMSbin.NEP_WM_QC= bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.NEP_WM_QC, X);
MSept_BMSbin.NEC_WM   = bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.NEC_WM, X);
MSept_BMSbin.NEC_WM_QC= bin_data_to_X_GF(MSept_BMS.SDN, MSept_BMS.NEC_WM_QC, X);


%% Plot Binned Fluxes 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP, 'b'); 
% NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC, 'r-'); 
% plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32Sept 2020 Hourly Binned Fluxes FluxFit');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC_QC, 'r-'); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32Sept 2020 Hourly Binned Fluxes FluxFit QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC_WM, 'r-'); 
% plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32Sept 2020 WM Original Hourly Binned Fluxes');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC_WM_QC, 'r-'); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32Sept 2020 WM Full QC Hourly Binned Fluxes');

%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
% DOplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dDOXY, 'c'); 
% TAplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dTA, 'k-'); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32Sept 2020 Binned Gradiets');



%% ***QC*** Remove sections when velocity is too slow
ibad = MSept_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow

% MSept_BMSbin.NEP(ibad) = NaN;
% MSept_BMSbin.NEC(ibad) = NaN;
MSept_BMSbin.NEP_QC(ibad) = NaN;
MSept_BMSbin.NEC_QC(ibad) = NaN;
MSept_BMSbin.dDOXY_QC(ibad) = NaN;
MSept_BMSbin.dTA_QC(ibad) = NaN;

% MSept_BMSbin.NEP_WM(ibad) = NaN;
% MSept_BMSbin.NEC_WM(ibad) = NaN;
MSept_BMSbin.NEP_WM_QC(ibad) = NaN;
MSept_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP, 'b'); 
% NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC, 'r-'); 
% plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32Sept 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC_QC, 'r-'); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32Sept 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC_WM, 'r-'); 
% plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('Sept Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32Sept 2020 Hourly Binned WM Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC_WM_QC, 'r-'); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32Sept 2020 Hourly Binned WM Fluxes Full QC');


%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
figure
hold on; box on;
boxplot(MSept_BMSbin.NEP_QC)
ylabel('NEP')
title('Sept NEP Boxplots')

figure
hold on; box on;
boxplot(MSept_BMSbin.NEC_QC)
ylabel('NEC')
title('Sept NEC Boxplots')

figure
hold on; box on;
boxplot(MSept_BMSbin.dDOXY_QC)
ylabel('DO')
title('Sept dDO Boxplots')

figure
hold on; box on;
boxplot(MSept_BMSbin.dTA_QC)
ylabel('TA')
title('Sept dTA Boxplots')

% Outliers: 

% 24: dTA and NEC outlier 
    %OUTLIER 
    
%106: dTA outlier (value: 2.2)
    %    '06-Sep-2020 09:00:00'
    %BAD PROFILE - OUTLIER 

%218: nec OUTLIER: -8.8
    %    '11-Sep-2020 01:00:00'
%217: NEC outlier: -8.3
%221: NEC outlier: -8
%214:: NEC outlier: -7.8
%228: NEC outlier: -7.5
%166: NEC outlier: -7.5
    %BAD PROFILE - OUTLIER 

datestr(MSept_BMSbin.SDN(166))

% Plot Profiles at outliers - 
close all 
for i = 166   %1:length(MSept_BMSbin.SDN)
    figure (i)
    scatter(MSept_BMSbin.uv(1:108,i), MSept_BMSbin.bin_depth(1:108))
    title(['Cudjoe Sept Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end
% outliers to be removed: 24 106 166
% MSept_BMSbin.NEP_QC(24) = NaN;
% MSept_BMSbin.NEC_QC(24) = NaN;
% MSept_BMSbin.dDOXY_QC(24) = NaN;
% MSept_BMSbin.dTA_QC(24) = NaN;
% 
% MSept_BMSbin.NEP_QC(106) = NaN;
% MSept_BMSbin.NEC_QC(106) = NaN;
% MSept_BMSbin.dDOXY_QC(106) = NaN;
% MSept_BMSbin.dTA_QC(106) = NaN;
% 
% MSept_BMSbin.NEP_QC(166) = NaN;
% MSept_BMSbin.NEC_QC(166) = NaN;
% MSept_BMSbin.dDOXY_QC(166) = NaN;
% MSept_BMSbin.dTA_QC(166) = NaN;


close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC_QC, 'r-'); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Oct Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 Oct 2020 Hourly Binned Fluxes Full QC');
%% Plot Diel Curves

% nepdbin = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEP, 24);
% necdbin = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEC, 24);
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

nepdbin_QC = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEC_QC, 24);
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
% nepdbin_WM = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEC_WM, 24);
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

nepdbin_WM_QC = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEC_WM_QC, 24);
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
% 
% 

%% Extract Daytime data for Ratios
clc
% Extract daytime data using MSept_BMSbin.PAR
MSept_inight = MSept_BMSbin.PAR(1,:) < 5; %find all nightime datapoints 

%create new arrays for daytime data
MSept_BMSbin.SDN_day = MSept_BMSbin.SDN;
MSept_BMSbin.PAR_day = MSept_BMSbin.PAR(1,:);

MSept_BMSbin.NEP_day = MSept_BMSbin.NEP;
MSept_BMSbin.NEC_day = MSept_BMSbin.NEC;
MSept_BMSbin.NEP_day_QC = MSept_BMSbin.NEP_QC;
MSept_BMSbin.NEC_day_QC = MSept_BMSbin.NEC_QC;

MSept_BMSbin.dDOXY_day_QC = MSept_BMSbin.dDOXY_QC;      %DO Gradient
MSept_BMSbin.dTA_day_QC = MSept_BMSbin.dTA_QC;          %TA Gradient

MSept_BMSbin.NEP_WM_day = MSept_BMSbin.NEP_WM;
MSept_BMSbin.NEC_WM_day = MSept_BMSbin.NEC_WM;
MSept_BMSbin.NEP_WM_day_QC = MSept_BMSbin.NEP_WM_QC;
MSept_BMSbin.NEC_WM_day_QC = MSept_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
MSept_BMSbin.SDN_day(MSept_inight) = NaN;
MSept_BMSbin.PAR_day (MSept_inight) = NaN;

MSept_BMSbin.NEP_day(MSept_inight) = NaN;
MSept_BMSbin.NEC_day(MSept_inight) = NaN;
MSept_BMSbin.NEP_day_QC(MSept_inight) = NaN;
MSept_BMSbin.NEC_day_QC(MSept_inight) = NaN;

MSept_BMSbin.dDOXY_day_QC(MSept_inight) = NaN;      %DO Gradient
MSept_BMSbin.dTA_day_QC(MSept_inight) = NaN;        %TA Gradient

MSept_BMSbin.NEP_WM_day(MSept_inight) = NaN;
MSept_BMSbin.NEC_WM_day(MSept_inight) = NaN;
MSept_BMSbin.NEP_WM_day_QC(MSept_inight) = NaN;
MSept_BMSbin.NEC_WM_day_QC(MSept_inight) = NaN;

%Plot to check only nighttime points removed
figure 
hold on
scatter(MSept_BMSbin.SDN, MSept_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MSept_BMSbin.SDN_day, MSept_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
MSept_BMSbin.NEP_day(isnan(MSept_BMSbin.NEP_day))=[];
MSept_BMSbin.NEC_day(isnan(MSept_BMSbin.NEC_day))=[];
MSept_BMSbin.NEP_day_QC(isnan(MSept_BMSbin.NEP_day_QC))=[];
MSept_BMSbin.NEC_day_QC(isnan(MSept_BMSbin.NEC_day_QC))=[];

MSept_BMSbin.dDOXY_day_QC(isnan(MSept_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
MSept_BMSbin.dTA_day_QC(isnan(MSept_BMSbin.dTA_day_QC))=[];       %TA Gradient

MSept_BMSbin.NEP_WM_day(isnan(MSept_BMSbin.NEP_WM_day))=[];
MSept_BMSbin.NEC_WM_day(isnan(MSept_BMSbin.NEC_WM_day))=[];
MSept_BMSbin.NEP_WM_day_QC(isnan(MSept_BMSbin.NEP_WM_day_QC))=[];
MSept_BMSbin.NEC_WM_day_QC(isnan(MSept_BMSbin.NEC_WM_day_QC))=[];

% create nighttime hours datasets
MSept_BMSbin.SDN_night = MSept_BMSbin.SDN;
MSept_BMSbin.PAR_night = MSept_BMSbin.PAR(1,:);

MSept_BMSbin.NEP_night = MSept_BMSbin.NEP;
MSept_BMSbin.NEC_night = MSept_BMSbin.NEC;
MSept_BMSbin.NEP_night_QC = MSept_BMSbin.NEP_QC;
MSept_BMSbin.NEC_night_QC = MSept_BMSbin.NEC_QC;

MSept_BMSbin.dDOXY_night_QC = MSept_BMSbin.dDOXY_QC;      %DO Gradient
MSept_BMSbin.dTA_night_QC = MSept_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
MSept_BMSbin.SDN_night=MSept_BMSbin.SDN_night(MSept_inight);
MSept_BMSbin.PAR_night=MSept_BMSbin.PAR_night(MSept_inight);

MSept_BMSbin.NEP_night_QC=MSept_BMSbin.NEP_night_QC(MSept_inight);
MSept_BMSbin.NEC_night_QC=MSept_BMSbin.NEC_night_QC(MSept_inight);

MSept_BMSbin.dDOXY_night_QC=MSept_BMSbin.dDOXY_night_QC(MSept_inight);      %DO Gradient
MSept_BMSbin.dTA_night_QC=MSept_BMSbin.dTA_night_QC(MSept_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(MSept_BMSbin.SDN, MSept_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MSept_BMSbin.SDN_night, MSept_BMSbin.PAR_night, 'r.'); % night plot

%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(MSept_BMSbin.NEP_day,MSept_BMSbin.NEC_day);
% MSept_BMSbin.Reg_Line = m*MSept_BMSbin.NEP_day + b;
% MSept_BMSbin.Ratio = m;
% MSept_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(MSept_BMSbin.NEP_day,MSept_BMSbin.NEC_day,'o')
% plot(MSept_BMSbin.NEP_day,MSept_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32Sept 2020 Pre-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MSept_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MSept_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(MSept_BMSbin.NEP_day_QC,MSept_BMSbin.NEC_day_QC);
MSept_BMSbin.Reg_Line_QC = m_QC*MSept_BMSbin.NEP_day_QC + b_QC;
MSept_BMSbin.Ratio_QC = m_QC;
MSept_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(MSept_BMSbin.NEP_day_QC,MSept_BMSbin.NEC_day_QC,'o')
plot(MSept_BMSbin.NEP_day_QC,MSept_BMSbin.Reg_Line_QC,'r')
%ylim([-50 50])
%xlim([-25 25])
xlabel('NCP');
ylabel('NCC');
title('Marker 32Sept 2020 Pre-Restoration NCC:NCP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MSept_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MSept_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(MSept_BMSbin.NEP_WM_day,MSept_BMSbin.NEC_WM_day);
% MSept_BMSbin.Reg_Line_WM = m_WM*MSept_BMSbin.NEP_WM_day + b_WM;
% MSept_BMSbin.Ratio_WM = m_WM;
% MSept_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(MSept_BMSbin.NEP_WM_day,MSept_BMSbin.NEC_WM_day,'o')
% plot(MSept_BMSbin.NEP_WM_day,MSept_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32Sept 2020 Pre-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MSept_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MSept_BMSbin.R2_WM)


[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(MSept_BMSbin.NEP_WM_day_QC,MSept_BMSbin.NEC_WM_day_QC);
MSept_BMSbin.Reg_Line_WM_QC = m_WM_QC*MSept_BMSbin.NEP_WM_day_QC + b_WM_QC;
MSept_BMSbin.Ratio_WM_QC = m_WM_QC;
MSept_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(MSept_BMSbin.NEP_WM_day_QC,MSept_BMSbin.NEC_WM_day_QC,'o')
% plot(MSept_BMSbin.NEP_WM_day_QC,MSept_BMSbin.Reg_Line_WM_QC,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32Sept 2020 Pre-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MSept_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MSept_BMSbin.R2_WM_QC)


%% For NEC:NEP Regressions Using Gradients
% close all 
clc

% multiply o2 gradient by -1 for O2 production
MSept_BMSbin.dDOXY_Reg = -1.*MSept_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
MSept_BMSbin.dTA_Reg = 0.5.*MSept_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(MSept_BMSbin.dDOXY_day_QC);
figure 
hold on 
DOday = plot(1:Xlength, MSept_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, MSept_BMSbin.dDOXY_Reg);
xlabel('Sept Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Marker 32Sept 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, MSept_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, MSept_BMSbin.dTA_Reg);
xlabel('Sept Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Marker 32Sept 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(MSept_BMSbin.dDOXY_Reg, MSept_BMSbin.dTA_Reg);
MSept_BMSbin.Reg_Line_G = m_G*MSept_BMSbin.dDOXY_Reg + b_G;
MSept_BMSbin.Ratio_G = m_G;
MSept_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(MSept_BMSbin.dDOXY_Reg,MSept_BMSbin.dTA_Reg,'o')
plot(MSept_BMSbin.dDOXY_Reg,MSept_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Marker 32Sept 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MSept_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MSept_BMSbin.R2_G)

clc
disp('Finished with Sept Metabolism Calculations');

save('MSept20_2.mat', 'MSept_BMS', 'MSept_BMSbin');

%% Plot Profiles 
% plot for profile within pump heights
% close all 
% for i =1:100 %length(MSept_ADavg.SDN)
%     figure (i)
%     scatter(MSept_ADavg.uv(1:108,i), MSept_ADavg.bin_depth(1:108))
%     title(['Marker 32Sept Velocity Profile Number ',num2str(i),])
%     xlabel('Velocity (m/s)');
%     ylabel('Height (m)');
% end


%% Subplots 
close all
clc

sgtitle('Marker 32 September 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
% DOplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dDOXY, 'c'); 
% TAplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.dTA, 'k-'); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'09/03 12:00';'09/04 12:00';'09/05 12:00';...
    '09/06 12:00';'09/07 12:00';'09/08 12:00';'09/09 12:00';'09/10 12:00';'09/11 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
% legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northwest');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(MSept_BMSbin.SDN, MSept_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(MSept_BMSbin.SDN, zeros(size(MSept_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MSept_good_Xrange, 'XTick', MSept_tick, 'xticklabel', MSept_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('Sept Days');
xticklabels({'09/03 12:00';'09/04 12:00';'09/05 12:00';...
    '09/06 12:00';'09/07 12:00';'09/08 12:00';'09/09 12:00';'09/10 12:00';'09/11 12:00'})
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MSept_BMSbin.SDN, MSept_BMSbin.NEC_QC, 24);
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
plot(MSept_BMSbin.NEP_day_QC,MSept_BMSbin.NEC_day_QC,'o')
plot(MSept_BMSbin.NEP_day_QC,MSept_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(MSept_BMSbin.Ratio_QC,2);
str2 = num2str(MSept_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(MSept_BMSbin.dDOXY_Reg,MSept_BMSbin.dTA_Reg,'o')
plot(MSept_BMSbin.dDOXY_Reg,MSept_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
title('dTA:dDO Ratio from Gradients');
str1 = num2str(MSept_BMSbin.Ratio_G,2);
str2 = num2str(MSept_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')

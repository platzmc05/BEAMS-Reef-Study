% July Analysis from Marker 32 Reef
% Michelle Platz - USF 
% 3/10/2021

% SeapHOx sensor deployed 7/01/2020 at 4pm EST
    % pump 1 height above benthos = 72 cm  
    % pump 2 height above benthos = 23 cm 
% ADP sensor deployed 7/03/2020 at 4pm EST
    % height from substrate to ADCP head = 18 cm

close all
clc
clear all

MJuly_z1 = 0.72;
MJuly_z2 = 0.23;
MJuly_ADheight = 0.18;
MJuly_ADbin_depth_1m = 1-MJuly_ADheight;
%% Initial look at data
% ***** create MJuly_SPraw data structure ***** observations every 30 seconds
%Parse SeapHOx data from datafile by variable 
MJuly_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('M_701BMS.txt');

%calculate O2 saturation concentration using temperature and salinity
MJuly_SPraw.DOXY = MJuly_SPraw.O2SATPER.*calcO2sat(MJuly_SPraw.MCAT_TC, MJuly_SPraw.PSAL)./100;

%calculate pH from durafet using internal reference electrode and Nernst equation 
MJuly_SPraw.pHint_prelim = calc_dfet_pHint(MJuly_SPraw.Vint, MJuly_SPraw.DFET_TC, -0.4);

% ***** create MJuly_SP data structure *****  observations every 15 mins
% sort data into respective pump heights
% daterange start must be first obs. of pump 1 cycle: pump 1/obs. 1
% daterange end must be end of pump 2 cycle: pump 2/obs.30
MJuly_SP = parse_to_pumpheights_ARM_2pump_Mar19(MJuly_SPraw, [datenum('07-01-2020 16:00:00'), datenum('07-26-2020 10:29:30')]);

% Calculate Gradients 
MJuly_SP = calc_TA_gradientV2(MJuly_SP, 2368.31, [0.8:0.1:1.2], 1, 2);
% Top TA is TA0 (estimated from average of discrete samples)
% calcualtes TA2, which is based on the Barnes equations.
% Q values tested: [0.8, 0.9, 1, 1.1, 1.2]

MJuly_SP.dDOXY = MJuly_SP.DOXY(1,:) - MJuly_SP.DOXY(2,:); %Oxygen Gradient
MJuly_SP.dpH = MJuly_SP.pH(1,:) - MJuly_SP.pH(2,:); %pH Gradient 
MJuly_SP.dTA = MJuly_SP.TAtop - MJuly_SP.TAbtm(3,:); % TA gradient - assuming Q=1

%% Plot Unbinned Gradients to determine good data Xrange
close all
clc
% Create Datestring for Plots
MJuly_DateString = {'07/01/2020 12:00:00';'07/02/2020 12:00:00';'07/03/2020 12:00:00';'07/04/2020 12:00:00';...
    '07/05/2020 12:00:00';'07/06/2020 12:00:00';'07/07/2020 12:00:00';'07/08/2020 12:00:00';'07/09/2020 12:00:00';...
    '07/10/2020 12:00:00';'07/11/2020 12:00:00';'07/12/2020 12:00:00';'07/13/2020 12:00:00';'07/14/2020 12:00:00';...
    '07/15/2020 12:00:00';'07/16/2020 12:00:00';'07/17/2020 12:00:00';'07/18/2020 12:00:00';'07/19/2020 12:00:00';...
    '07/20/2020 12:00:00';'07/21/2020 12:00:00';'07/22/2020 12:00:00';'07/23/2020 12:00:00';'07/24/2020 12:00:00';...
    '07/25/2020 12:00:00'};

formatIn = 'mm/dd/yyyy HH:MM:SS';
MJuly_tick = datenum(MJuly_DateString,formatIn);

MJuly_Xrange = [datenum('07-01-2020 16:00:00'), datenum('07-26-2020 10:29:30')];

figure
hold on; box on;
plot(MJuly_SP.SDN, MJuly_SP.dDOXY); %oxygen gradient 
plot(MJuly_SP.SDN, MJuly_SP.dTA); %TA gradient
plot(MJuly_SP.SDN, zeros(size(MJuly_SP.SDN))); %zero line
set(gca, 'xlim', MJuly_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 July 2020 Unbinned DO and TA Gradients');

%% Save Full Site Characterization Datasets
% save from SPraw to get 30 sec measurement intervals
MJuly_SiteChar.SDN = MJuly_SPraw.SDN;
MJuly_SiteChar.TC = MJuly_SPraw.OPT_TC;
MJuly_SiteChar.PAR = MJuly_SPraw.PAR;
MJuly_SiteChar.PSAL = MJuly_SPraw.PSAL;
MJuly_SiteChar.Pres = MJuly_SPraw.Pres;

close all
figure 
hold on; 
plot(MJuly_SiteChar.SDN, MJuly_SiteChar.Pres)

%clip ends of data to remove surfave interval observations 
MJuly_SiteChar.SDN = MJuly_SPraw.SDN(188:71636);
MJuly_SiteChar.TC = MJuly_SPraw.OPT_TC(188:71636);
MJuly_SiteChar.PAR = MJuly_SPraw.PAR(188:71636);
MJuly_SiteChar.PSAL = MJuly_SPraw.PSAL(188:71636);
MJuly_SiteChar.Pres = MJuly_SPraw.Pres(188:71636);

%extract full length of ADCP datafile  
MJuly_ADfull=aquadoppraw2mat('M_7_03', 70, [datenum('07-03-2020 16:00:00'), datenum('07-26-2020 16:00:00')]);

% add AD variables to Site Char 
MJuly_SiteChar.AD_SDN = MJuly_ADfull.SDN;
MJuly_SiteChar.AD_Pres = MJuly_ADfull.Pres;
MJuly_SiteChar.AD_TC = MJuly_ADfull.TC;
MJuly_SiteChar.bin_depth = MJuly_ADfull.bin_depth;
MJuly_SiteChar.u = MJuly_ADfull.u;
MJuly_SiteChar.v = MJuly_ADfull.v;
MJuly_SiteChar.w = MJuly_ADfull.w;
MJuly_SiteChar.uv = MJuly_ADfull.uv;
MJuly_SiteChar.direction = MJuly_ADfull.direction;

%Plot to see when surface interval observations are
close all
figure 
hold on; 
plot(MJuly_SiteChar.AD_SDN, MJuly_SiteChar.AD_Pres)

%clip ends of data to remove surfave interval observations 
MJuly_SiteChar.AD_SDN = MJuly_ADfull.SDN(1:65636);
MJuly_SiteChar.AD_Pres = MJuly_ADfull.Pres(1:65636);
MJuly_SiteChar.AD_TC = MJuly_ADfull.TC(1:65636);
MJuly_SiteChar.bin_depth = MJuly_ADfull.bin_depth;
MJuly_SiteChar.u = MJuly_ADfull.u(:,1:65636);
MJuly_SiteChar.v = MJuly_ADfull.v(:,1:65636);
MJuly_SiteChar.w = MJuly_ADfull.w(:,1:65636);
MJuly_SiteChar.uv = MJuly_ADfull.uv(:,1:65636);
MJuly_SiteChar.direction = MJuly_ADfull.direction(:,1:65636);

% find U0 
MJuly_z1 = 0.72;
MJuly_z2 = 0.23;
MJuly_ADheight = 0.18;
MJuly_ADbin_depth_1m = 1-MJuly_ADheight;
MJuly_i1m = find(MJuly_SiteChar.bin_depth==MJuly_ADbin_depth_1m);
MJuly_SiteChar.U0 = MJuly_SiteChar.uv(72,:);

% save data in separate datastructure
save('MJuly20_SiteChar_2.mat', 'MJuly_SiteChar')



%% Constrain Xrange from graph results and extract good gradient data - 
close all 

MJuly_good_Xrange = [datenum('07-01-2020 16:00:00'), datenum('07-15-2020 14:00:00')]; 

% plot to check range is correct
figure
hold on; box on;
plot(MJuly_SP.SDN, MJuly_SP.dDOXY); %oxygen gradient 
plot(MJuly_SP.SDN, MJuly_SP.dTA); %TA gradient
plot(MJuly_SP.SDN, zeros(size(MJuly_SP.SDN))); %zero line
set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 July 2020 Unbinned DO and TA Gradients');

%% Create good dataframe

close all

M_July_BMS_idx_start = find(MJuly_SP.SDN==datenum('07-01-2020 16:00:00'))
M_July_BMS_idx_end = find(MJuly_SP.SDN==datenum('07-15-2020 13:00:00'))

% Create new data vectors of just the good data
M_July_BMS_good_data = M_July_BMS_idx_start:M_July_BMS_idx_end;
Initial_data_points = length(M_July_BMS_good_data)

%% Extract good data for all SeapHOx Parameters

clc

vars = fieldnames(MJuly_SP);
for v = 1:length(vars)
    MJuly_SP.(vars{v}) = (MJuly_SP.(vars{v})(:,M_July_BMS_good_data));
end
    
%% *************** ADCP DATA ****************
% ***** create new data structure: MJuly_AD *****
%***** ADCP deployed later in the week on 7/3/2020 @16:00*****

clc
close all 
% data points every 30 seconds
% pull only good dataframe identified above
MJuly_AD=aquadoppraw2mat('M_7_03', 70, [datenum('07-03-2020 18:00:00'), datenum('07-15-2020 13:15:00')]);

%averages data to the middle of the minute interval spacified 
MJuly_ADavg = average_aquadopp(MJuly_AD, 15.1);

%% Calc ustar 
% calculates ustar from current profiles 
% actual heights  = 0.72m (pump 1) and 0.23m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.54m (bin 44) and 0.05m (bin 1)  - bins from which to pull ADCP data 
% salinity - estimated from mean of SP_Sal data over observation period -

clc
% already removed data outside data frame so can take average of whole set
MJuly_Sal_est = mean(MJuly_SP.PSAL(1,3:end));

[MJuly_ADavg] = ustar_from_aquadopp2(MJuly_ADavg,[0.54 0.11], MJuly_Sal_est); %bins adjusted 

%[ADavg] = ustar_McGillis_Method(ADavg, ztop, zbtm, bintop, binbtm)
[MJuly_ADavg] = ustar_McGillis_Method(MJuly_ADavg, 0.72, 0.23, 44, 1);

%compare ustar calculation methods
close all
figure
hold on 
ustar_plot = plot(MJuly_ADavg.SDN, MJuly_ADavg.ustar, 'r');
ustar_WM_plot = plot(MJuly_ADavg.SDN, MJuly_ADavg.ustar_WM);
plot(MJuly_ADavg.SDN, zeros(size(MJuly_ADavg.SDN)),'k');
set(gca, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MJuly_ADavg.SDN(1) MJuly_ADavg.SDN(end)])
xlabel('July Days');
ylabel('ustar values');
legend([ustar_plot ustar_WM_plot], {'ustar plot','ustar WM plot'}, 'location', 'northeast');
title('Marker 32 July 2020 Ustar Values');


%% Combine SP and AD data into one data structure  
%***** create new data structure: MJuly_BMS *****

ADavg_vars = fieldnames(MJuly_ADavg);
for v = 1:length(ADavg_vars)
    MJuly_BMS.(ADavg_vars{v}) = (MJuly_ADavg.(ADavg_vars{v}));
end

% SP second to override SDN
%Becuase ADCP was deployed late, need to extract coresponding SP data
% [datenum('07-03-2020 18:00:00'), datenum('07-15-2020 13:00:00')]
SP_short_start = find(MJuly_SP.SDN==datenum('07-03-2020 18:00:00'))
SP_vars = fieldnames(MJuly_SP);
for v = 1:length(SP_vars)
    MJuly_BMS.(SP_vars{v}) = (MJuly_SP.(SP_vars{v})(:,201:1333));
end
% check that SDN is on 15 min interval
datestr(MJuly_BMS.SDN)

%% %% *************** Calculate Fluxes ****************

% actual heights  = 0.72m (pump 1) and 0.23m (pump 2) 
% 0.18m from substrate to ACDP head - 
% adjusted height = 0.54m (bin 44) and 0.05m (bin 1)  - bins from which to pull ADCP data 

clc

%NCC - calculates TA flux and NCC from ustar and TA concetration gradients
[MJuly_BMS] = calc_NCC_3(MJuly_BMS,[0.54 0.11]);

%NCP - calculates DO flux and NCP from ustar and DO concetration gradients
C1guess = median(MJuly_BMS.DOXY(1,:))
[MJuly_BMS] = calc_NCP_3(MJuly_BMS, [0.54 0.11],179); %estimate C1 guess using median DOXY(1,:) value  

% Plot NCP and NCC
%close all
figure
hold on; box on; 
NEPplot = plot(MJuly_BMS.SDN, MJuly_BMS.NEP);
NECplot = plot(MJuly_BMS.SDN, MJuly_BMS.NEC);
plot(MJuly_BMS.SDN, zeros(size(MJuly_BMS.SDN)),'k');
set(gca, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MJuly_BMS.SDN(1) MJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32 July 2020 Fluxes');

%McGillis method flux calculations 
[MJuly_BMS] = calc_NCP_McGillis_Method(MJuly_BMS, 0.72, 0.23, MJuly_Sal_est);
[MJuly_BMS] = calc_NCC_McGillis_Method(MJuly_BMS, 0.72, 0.23, MJuly_Sal_est);

% Plot NCP and NCC
% close all
figure
hold on; box on; 
NEPplot = plot(MJuly_BMS.SDN, MJuly_BMS.NEP_WM);
NECplot = plot(MJuly_BMS.SDN, MJuly_BMS.NEC_WM);
plot(MJuly_BMS.SDN, zeros(size(MJuly_BMS.SDN)),'k');
set(gca, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MJuly_BMS.SDN(1) MJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCP or NCC [mmol/m2/hr]');
legend([NEPplot NECplot], {'NCP','NCC'}, 'location', 'northeast');
title('Marker 32 July 2020 WM Fluxes');

%% Compare Flux_fit vs WM Plots 

% NCP Plot  
close all
figure
hold on; box on; 
NEPplot = plot(MJuly_BMS.SDN, MJuly_BMS.NEP);
NEPplotWM = plot(MJuly_BMS.SDN, MJuly_BMS.NEP_WM);
plot(MJuly_BMS.SDN, zeros(size(MJuly_BMS.SDN)),'k');
set(gca, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MJuly_BMS.SDN(1) MJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotWM], {'NEP','NEP WM'}, 'location', 'northeast');
title('M32 July 2020 Fluxes');

%NCC plot 
figure
hold on; box on; 
NECplot = plot(MJuly_BMS.SDN, MJuly_BMS.NEC);
NECplotWM = plot(MJuly_BMS.SDN, MJuly_BMS.NEC_WM);
plot(MJuly_BMS.SDN, zeros(size(MJuly_BMS.SDN)),'k');
set(gca, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MJuly_BMS.SDN(1) MJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCC [mmol/m2/hr]');
legend([NECplot NECplotWM], {'NEC','NEC WM'}, 'location', 'northeast');
title('M32 July 2020 Fluxes');


%% ***QC*** Find when the stdev of DO is > 2 umol/kg at a given pump height, 
%indicates boundary layer was non-steady state and therefore unfit for gradient flux analysis 
clc
%calculate standard deviation of each DOXY observation
MJuly_BMS.DOXYstd = std(MJuly_BMS.DOXY);

% get DOXY std
MJuly_idoxystd = find(MJuly_BMS.DOXYstd > 2);% 58.4% of data is greater than 0.8 stdev
MJuly_ihighdoxystd = [];
for i = 1:length(MJuly_idoxystd)
    
    MJuly_ihighdoxystd = vertcat(MJuly_ihighdoxystd,[MJuly_idoxystd(i)-1:1:MJuly_idoxystd(i)+1]');
end
% get unique IDs
MJuly_ihighdoxystd = unique(MJuly_ihighdoxystd);
% remove 0's and out of index values
MJuly_ihighdoxystd(MJuly_ihighdoxystd==0) = [];
MJuly_ihighdoxystd(MJuly_ihighdoxystd> length(MJuly_BMS.SDN)) = [];

% make it into index
trex = false(size(MJuly_BMS.SDN));
trex(MJuly_ihighdoxystd) = true;
MJuly_ihighdoxystd = trex;
clear trex;

%create _QC datasets to preserve original data
MJuly_BMS.NEP_QC = MJuly_BMS.NEP;
MJuly_BMS.NEC_QC = MJuly_BMS.NEC;
MJuly_BMS.dDOXY_QC = MJuly_BMS.dDOXY; %DO gradient
MJuly_BMS.dTA_QC = MJuly_BMS.dTA;     %TA gradient
MJuly_BMS.NEP_WM_QC = MJuly_BMS.NEP_WM;
MJuly_BMS.NEC_WM_QC = MJuly_BMS.NEC_WM;

% set observations when DOXYstd>X to NaN - only remove NaNs from QC datasets
MJuly_BMS.NEP_QC(MJuly_ihighdoxystd) = NaN;
MJuly_BMS.NEC_QC(:,MJuly_ihighdoxystd) = NaN;
MJuly_BMS.dDOXY_QC(MJuly_ihighdoxystd) = NaN; %DO gradient
MJuly_BMS.dTA_QC(:,MJuly_ihighdoxystd) = NaN;   %TA gradient
MJuly_BMS.NEP_WM_QC(MJuly_ihighdoxystd) = NaN;
MJuly_BMS.NEC_WM_QC(:,MJuly_ihighdoxystd) = NaN;

% plot to see what got removed
close all
figure
hold on; box on;
NEPplot = plot(MJuly_BMS.SDN, MJuly_BMS.NEP, 'k');
NEPplotQC = plot(MJuly_BMS.SDN, MJuly_BMS.NEP_QC, 'r', 'linewidth', 1.5);
plot(MJuly_BMS.SDN, zeros(size(MJuly_BMS.SDN)),'k');
set(gca, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MJuly_BMS.SDN(1) MJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32 July 2020 Fluxes');


%WM Plot
%close all
figure
hold on; box on;
NEPplot = plot(MJuly_BMS.SDN, MJuly_BMS.NEP_WM, 'k');
NEPplotQC = plot(MJuly_BMS.SDN, MJuly_BMS.NEP_WM_QC, 'r', 'linewidth', 1.5);
plot(MJuly_BMS.SDN, zeros(size(MJuly_BMS.SDN)),'k');
set(gca, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlim([MJuly_BMS.SDN(1) MJuly_BMS.SDN(end)])
xlabel('July Days');
ylabel('NCP [mmol/m2/hr]');
legend([NEPplot NEPplotQC], {'NEP removed','NEP QC'}, 'location', 'northeast');
title('Marker 32 July 2020 Fluxes');


%% Bin data to hourly intervals

X = floor(nanmin(MJuly_BMS.SDN)):1/24:ceil(nanmax(MJuly_BMS.SDN));

MJuly_BMSbin.SDN = X;
% variables from 2 differnet heights
MJuly_BMSbin.DOXY(1,:)     = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.DOXY(1,:), X);
MJuly_BMSbin.DOXY(2,:)     = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.DOXY(2,:), X);

MJuly_BMSbin.pH(1,:)       = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.pH(1,:), X);
MJuly_BMSbin.pH(2,:)       = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.pH(2,:), X);

MJuly_BMSbin.PSAL(1,:)     = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.PSAL(1,:), X);
MJuly_BMSbin.PSAL(2,:)     = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.PSAL(2,:), X);

MJuly_BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.O2SATPER(1,:), X);
MJuly_BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.O2SATPER(2,:), X);

MJuly_BMSbin.Pres(1,:)     = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.Pres(1,:), X);
MJuly_BMSbin.Pres(2,:)     = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.Pres(2,:), X);

MJuly_BMSbin.DENS(1,:)     = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.DENS(1,:), X);
MJuly_BMSbin.DENS(2,:)     = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.DENS(2,:), X);

MJuly_BMSbin.PAR(1,:)      = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.PAR(1,:), X);
MJuly_BMSbin.PAR(2,:)      = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.PAR(2,:), X);

MJuly_BMSbin.bin_depth     = MJuly_BMS.bin_depth;

for i = 1:108
    MJuly_BMSbin.uv(i,:)   = bin_data_to_X_GF(MJuly_BMS.SDN,MJuly_BMS.uv(i,:), X);
end

% bin data hourly. Vector variables 
MJuly_BMSbin.PRES  = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.Pres, X);
MJuly_BMSbin.U0       = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.U0, X);
MJuly_BMSbin.DIR      = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.direction, X);
MJuly_BMSbin.ustar    = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.ustar, X);
MJuly_BMSbin.ustar_rm = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.ustar_runmean, X);
MJuly_BMSbin.dTA      = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.dTA, X);
MJuly_BMSbin.dTA_QC   = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.dTA_QC, X);
MJuly_BMSbin.dpH      = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.dpH, X);
MJuly_BMSbin.dDOXY    = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.dDOXY, X);
MJuly_BMSbin.dDOXY_QC = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.dDOXY_QC, X);
MJuly_BMSbin.NEP      = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.NEP, X);
MJuly_BMSbin.NEP_QC   = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.NEP_QC, X);
MJuly_BMSbin.NEC      = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.NEC, X);
MJuly_BMSbin.NEC_QC   = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.NEC_QC, X);
MJuly_BMSbin.NEP_WM   = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.NEP_WM, X);
MJuly_BMSbin.NEP_WM_QC= bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.NEP_WM_QC, X);
MJuly_BMSbin.NEC_WM   = bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.NEC_WM, X);
MJuly_BMSbin.NEC_WM_QC= bin_data_to_X_GF(MJuly_BMS.SDN, MJuly_BMS.NEC_WM_QC, X);


%% Plot Binned Fluxes 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP, 'b'); 
% NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC, 'r-'); 
% plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 July 2020 Hourly Binned Fluxes FluxFit');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_QC, 'r-'); 
plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 July 2020 Hourly Binned Fluxes FluxFit QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_WM, 'r-'); 
% plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 July 2020 WM Original Hourly Binned Fluxes');

% close all
clc
figure
hold on; box on;
NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_WM_QC, 'b'); 
NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_WM_QC, 'r-'); 
plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 July 2020 WM Full QC Hourly Binned Fluxes');

%% Plot Binned Gradients 
close all
figure
hold on; box on;
DOplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
% DOplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dDOXY, 'c'); 
% TAplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dTA, 'k-'); 
plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 July 2020 Binned Gradiets');



%% ***QC*** Remove sections when velocity is too slow
ibad = MJuly_BMSbin.U0 < 0.03; % when velociy at 1m above substrate is too slow

% MJuly_BMSbin.NEP(ibad) = NaN;
% MJuly_BMSbin.NEC(ibad) = NaN;
MJuly_BMSbin.NEP_QC(ibad) = NaN;
MJuly_BMSbin.NEC_QC(ibad) = NaN;

MJuly_BMSbin.dDOXY_QC(ibad) = NaN;
MJuly_BMSbin.dTA_QC(ibad) = NaN;

% MJuly_BMSbin.NEP_WM(ibad) = NaN;
% MJuly_BMSbin.NEC_WM(ibad) = NaN;
MJuly_BMSbin.NEP_WM_QC(ibad) = NaN;
MJuly_BMSbin.NEC_WM_QC(ibad) = NaN;

% Plot to see what was removed 
close all
% clc
% figure
% hold on; box on;
% NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP, 'b'); 
% NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC, 'r-'); 
% plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 July 2020 Hourly Binned Fluxes');

clc
figure
hold on; box on;
NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_QC, 'r-'); 
plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 July 2020 Hourly Binned Fluxes Full QC');


% WM plots
% figure
% hold on; box on;
% NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_WM, 'b'); 
% NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_WM, 'r-'); 
% plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 July 2020 Hourly Binned WM Fluxes');

% clc
% figure
% hold on; box on;
% NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_WM_QC, 'b'); 
% NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_WM_QC, 'r-'); 
% plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
% set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
% datetick('x', 'dd', 'keeplimits', 'keepticks');
% xlabel('July Days');
% ylabel('NEP or NEC');
% legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
% title('Marker 32 July 2020 Hourly Binned WM Fluxes Full QC');


%% Boxplots - Identify remaining outliers
% Fluxfit Calcs
figure
hold on; box on;
boxplot(MJuly_BMSbin.NEP_QC)
ylabel('NEP')
title('July NEP Boxplots')

figure
hold on; box on;
boxplot(MJuly_BMSbin.NEC_QC)
ylabel('NEC')
title('July NEC Boxplots')

figure
hold on; box on;
boxplot(MJuly_BMSbin.dDOXY_QC)
ylabel('DO')
title('July dDO Boxplots')

figure
hold on; box on;
boxplot(MJuly_BMSbin.dTA_QC)
ylabel('TA')
title('July dTA Boxplots')

% Outliers to investigate

% 66: NEP outlier (value: -63)NEC outlier (value: 9.5)and dTA
    %OUTLIER
    
% 67: NEP outlier (value: -23)

% 280: NEP outlier (value: -43)NEC outlier (value: -7.7)and dTA
    %OUTLIER

% 114: NEP outlier (value: -30)
    %OUTLIER - 5pm - value unlikely

% 115: NEP outlier (value: -23)
    %OUTLIER- 4pm - value unlikely

% 284: NEP outlier (value: -23) NEC outlier (value: -5.8)and dTA
%    '14-Jul-2020 19:00:00'
    %OUTLIER
% 283: NEC outlier (value: -5.8)and dTA
    %OUTLIER

% 37: NEC outlier (value: 8) and dTA
    %OUTLIER
    
% 38: NEC outlier (value: 16) and dTA
    %OUTLIER
    
% 86: NEC outlier (value: 6.7)
% 77: NEC outlier (value: 6)and dTA
% 254: NEC outlier (value: 6)and dTA
% 85: NEC outlier (value: 5.6)
% 43: NEC outlier (value: 5.6)
% 202: NEC outlier (value: -8)
    %OUTLIER

% 238: NEC outlier (value: -5.8)
    %OUTLIER - bad profile

% 187: dTA outlier (value: 

datestr(MJuly_BMSbin.SDN(114))

% Plot Profiles at outliers -  
close all 
for i = 67  %1:length(MAug2_BMSbin.SDN)
    figure (i)
    scatter(MJuly_BMSbin.uv(1:108,i), MJuly_BMSbin.bin_depth(1:108))
    title(['Marker 32 July Velocity Profile Number ',num2str(i),])
    xlabel('Velocity (m/s)');
    ylabel('Height (m)');
end


%outliers to be removed: 37 38 66 202 238 280 283 284 
% MJuly_BMSbin.NEP_QC(37) = NaN;
% MJuly_BMSbin.NEC_QC(37) = NaN;
% MJuly_BMSbin.dDOXY_QC(37) = NaN;
% MJuly_BMSbin.dTA_QC(37) = NaN;
% 
% MJuly_BMSbin.NEP_QC(38) = NaN;
% MJuly_BMSbin.NEC_QC(38) = NaN;
% MJuly_BMSbin.dDOXY_QC(38) = NaN;
% MJuly_BMSbin.dTA_QC(38) = NaN;
% 
% MJuly_BMSbin.NEP_QC(66) = NaN;
% MJuly_BMSbin.NEC_QC(66) = NaN;
% MJuly_BMSbin.dDOXY_QC(66) = NaN;
% MJuly_BMSbin.dTA_QC(66) = NaN;
% 
% MJuly_BMSbin.NEP_QC(202) = NaN;
% MJuly_BMSbin.NEC_QC(202) = NaN;
% MJuly_BMSbin.dDOXY_QC(202) = NaN;
% MJuly_BMSbin.dTA_QC(202) = NaN;
% 
% MJuly_BMSbin.NEP_QC(238) = NaN;
% MJuly_BMSbin.NEC_QC(238) = NaN;
% MJuly_BMSbin.dDOXY_QC(238) = NaN;
% MJuly_BMSbin.dTA_QC(238) = NaN;
% 
% MJuly_BMSbin.NEP_QC(280) = NaN;
% MJuly_BMSbin.NEC_QC(280) = NaN;
% MJuly_BMSbin.dDOXY_QC(280) = NaN;
% MJuly_BMSbin.dTA_QC(280) = NaN;
% 
% MJuly_BMSbin.NEP_QC(283) = NaN;
% MJuly_BMSbin.NEC_QC(283) = NaN;
% MJuly_BMSbin.dDOXY_QC(283) = NaN;
% MJuly_BMSbin.dTA_QC(283) = NaN;
% 
% MJuly_BMSbin.NEP_QC(284) = NaN;
% MJuly_BMSbin.NEC_QC(284) = NaN;
% MJuly_BMSbin.dDOXY_QC(284) = NaN;
% MJuly_BMSbin.dTA_QC(284) = NaN;
% 
% MJuly_BMSbin.NEP_QC(114) = NaN;
% MJuly_BMSbin.NEC_QC(114) = NaN;
% MJuly_BMSbin.dDOXY_QC(114) = NaN;
% MJuly_BMSbin.dTA_QC(114) = NaN;
% 
% MJuly_BMSbin.NEP_QC(115) = NaN;
% MJuly_BMSbin.NEC_QC(115) = NaN;
% MJuly_BMSbin.dDOXY_QC(115) = NaN;
% MJuly_BMSbin.dTA_QC(115) = NaN;

close all
clc
figure
subplot(2,1,1)
hold on; box on;
DOplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dDOXY_QC, 'b', 'linewidth', 1.5); 
TAplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dTA_QC, 'r-', 'linewidth', 1.5); 
plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('dDO or dTA');
legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Marker 32 July 2020 Binned Gradiets');

subplot(2,1,2)
hold on; box on;
NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_QC, 'b'); 
NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_QC, 'r-'); 
plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MJuly_good_Xrange, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
ylabel('NEP or NEC');
legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'northeast');
title('Marker 32 July 2020 Hourly Binned Fluxes Full QC');


%% Plot Diel Curves

% nepdbin = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP, 24);
% necdbin = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin, 'bo', 'markersize', 3);
% plot(1:24, necdbin, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin,1), 'bo-');
% plot(1:24, nanmedian(necdbin,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('July Diel Plot 2020');
% xlabel('hour of day');

nepdbin_QC = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_QC, 24);
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
title('July Diel Plot 2020  QC');
xlabel('hour of day');

% WM data
% nepdbin_WM = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_WM, 24);
% necdbin_WM = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_WM, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM,1), 'ro-')
% ylabel(['NEP or \color{red}NEC']);
% title('July 2020 WM');
% xlabel('hour of day');

nepdbin_WM_QC = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_WM_QC, 24);
necdbin_WM_QC = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_WM_QC, 24);
% figure
% hold on; box on;
% plot(1:24, zeros(size(1:24)), 'k:');
% plot(1:24, nepdbin_WM_QC, 'bo', 'markersize', 3);
% plot(1:24, necdbin_WM_QC, 'ro', 'markersize', 3);
% plot(1:24, nanmedian(nepdbin_WM_QC,1), 'bo-');
% plot(1:24, nanmedian(necdbin_WM_QC,1), 'ro-')
% xlim([1 24]);
% ylabel(['NEP or \color{red}NEC']);
% title('July 2020 WM Full QC');
% xlabel('hour of day');



%% Extract Daytime data for Ratios
clc
% Extract daytime data using MJuly_BMSbin.PAR
MJuly_inight = MJuly_BMSbin.PAR(1,:) < 1; %find all nightime datapoints 

%create new arrays for daytime data
MJuly_BMSbin.SDN_day = MJuly_BMSbin.SDN;
MJuly_BMSbin.PAR_day = MJuly_BMSbin.PAR(1,:);

MJuly_BMSbin.NEP_day = MJuly_BMSbin.NEP;
MJuly_BMSbin.NEC_day = MJuly_BMSbin.NEC;
MJuly_BMSbin.NEP_day_QC = MJuly_BMSbin.NEP_QC;
MJuly_BMSbin.NEC_day_QC = MJuly_BMSbin.NEC_QC;

MJuly_BMSbin.dDOXY_day_QC = MJuly_BMSbin.dDOXY_QC;      %DO Gradient
MJuly_BMSbin.dTA_day_QC = MJuly_BMSbin.dTA_QC;          %TA Gradient

MJuly_BMSbin.NEP_WM_day = MJuly_BMSbin.NEP_WM;
MJuly_BMSbin.NEC_WM_day = MJuly_BMSbin.NEC_WM;
MJuly_BMSbin.NEP_WM_day_QC = MJuly_BMSbin.NEP_WM_QC;
MJuly_BMSbin.NEC_WM_day_QC = MJuly_BMSbin.NEC_WM_QC;

%set all nightime values to NaN
MJuly_BMSbin.SDN_day(MJuly_inight) = NaN;
MJuly_BMSbin.PAR_day (MJuly_inight) = NaN;

MJuly_BMSbin.NEP_day(MJuly_inight) = NaN;
MJuly_BMSbin.NEC_day(MJuly_inight) = NaN;
MJuly_BMSbin.NEP_day_QC(MJuly_inight) = NaN;
MJuly_BMSbin.NEC_day_QC(MJuly_inight) = NaN;

MJuly_BMSbin.dDOXY_day_QC(MJuly_inight) = NaN;      %DO Gradient
MJuly_BMSbin.dTA_day_QC(MJuly_inight) = NaN;        %TA Gradient

MJuly_BMSbin.NEP_WM_day(MJuly_inight) = NaN;
MJuly_BMSbin.NEC_WM_day(MJuly_inight) = NaN;
MJuly_BMSbin.NEP_WM_day_QC(MJuly_inight) = NaN;
MJuly_BMSbin.NEC_WM_day_QC(MJuly_inight) = NaN;

%Plot to check only nighttime points removed
figure 
hold on
scatter(MJuly_BMSbin.SDN, MJuly_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MJuly_BMSbin.SDN_day, MJuly_BMSbin.PAR_day, 'r.'); % day plot

%Remove NaN values from fluxes
MJuly_BMSbin.NEP_day(isnan(MJuly_BMSbin.NEP_day))=[];
MJuly_BMSbin.NEC_day(isnan(MJuly_BMSbin.NEC_day))=[];
MJuly_BMSbin.NEP_day_QC(isnan(MJuly_BMSbin.NEP_day_QC))=[];
MJuly_BMSbin.NEC_day_QC(isnan(MJuly_BMSbin.NEC_day_QC))=[];

MJuly_BMSbin.dDOXY_day_QC(isnan(MJuly_BMSbin.dDOXY_day_QC))=[];   %DO Gradient
MJuly_BMSbin.dTA_day_QC(isnan(MJuly_BMSbin.dTA_day_QC))=[];       %TA Gradient

MJuly_BMSbin.NEP_WM_day(isnan(MJuly_BMSbin.NEP_WM_day))=[];
MJuly_BMSbin.NEC_WM_day(isnan(MJuly_BMSbin.NEC_WM_day))=[];
MJuly_BMSbin.NEP_WM_day_QC(isnan(MJuly_BMSbin.NEP_WM_day_QC))=[];
MJuly_BMSbin.NEC_WM_day_QC(isnan(MJuly_BMSbin.NEC_WM_day_QC))=[];

% create nighttime hours datasets
MJuly_BMSbin.SDN_night = MJuly_BMSbin.SDN;
MJuly_BMSbin.PAR_night = MJuly_BMSbin.PAR(1,:);

MJuly_BMSbin.NEP_night = MJuly_BMSbin.NEP;
MJuly_BMSbin.NEC_night = MJuly_BMSbin.NEC;
MJuly_BMSbin.NEP_night_QC = MJuly_BMSbin.NEP_QC;
MJuly_BMSbin.NEC_night_QC = MJuly_BMSbin.NEC_QC;

MJuly_BMSbin.dDOXY_night_QC = MJuly_BMSbin.dDOXY_QC;      %DO Gradient
MJuly_BMSbin.dTA_night_QC = MJuly_BMSbin.dTA_QC;          %TA Gradient

% extract nighttime hours
MJuly_BMSbin.SDN_night=MJuly_BMSbin.SDN_night(MJuly_inight);
MJuly_BMSbin.PAR_night=MJuly_BMSbin.PAR_night(MJuly_inight);

MJuly_BMSbin.NEP_night_QC=MJuly_BMSbin.NEP_night_QC(MJuly_inight);
MJuly_BMSbin.NEC_night_QC=MJuly_BMSbin.NEC_night_QC(MJuly_inight);

MJuly_BMSbin.dDOXY_night_QC=MJuly_BMSbin.dDOXY_night_QC(MJuly_inight);      %DO Gradient
MJuly_BMSbin.dTA_night_QC=MJuly_BMSbin.dTA_night_QC(MJuly_inight);        %TA Gradient


%Plot to check only nighttime points removed
figure 
hold on
scatter(MJuly_BMSbin.SDN, MJuly_BMSbin.PAR(1,:), 'o');% day/night plot
scatter(MJuly_BMSbin.SDN_night, MJuly_BMSbin.PAR_night, 'r.'); % night plot


%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 

close all 
clc

% [m,b,r,sm,sb]=lsqfitgm(MJuly_BMSbin.NEP_day,MJuly_BMSbin.NEC_day);
% MJuly_BMSbin.Reg_Line = m*MJuly_BMSbin.NEP_day + b;
% MJuly_BMSbin.Ratio = m;
% MJuly_BMSbin.R2 = r;
% % plot
% figure
% hold on; box on;
% plot(MJuly_BMSbin.NEP_day,MJuly_BMSbin.NEC_day,'o')
% plot(MJuly_BMSbin.NEP_day,MJuly_BMSbin.Reg_Line,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 July 2020 Pre-Restoration NCC:NCP Ratio FluxFit');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MJuly_BMSbin.Ratio)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MJuly_BMSbin.R2)


[m_QC,b_QC,r_QC,sm_QC,sb_QC]=lsqfitgm(MJuly_BMSbin.NEP_day_QC,MJuly_BMSbin.NEC_day_QC);
MJuly_BMSbin.Reg_Line_QC = m_QC*MJuly_BMSbin.NEP_day_QC + b_QC;
MJuly_BMSbin.Ratio_QC = m_QC;
MJuly_BMSbin.R2_QC = r_QC;
% plot
figure
hold on; box on;
plot(MJuly_BMSbin.NEP_day_QC,MJuly_BMSbin.NEC_day_QC,'o')
plot(MJuly_BMSbin.NEP_day_QC,MJuly_BMSbin.Reg_Line_QC,'r')
%ylim([-50 50])
%xlim([-25 25])
xlabel('NCP');
ylabel('NCC');
title('Marker 32 July 2020 Pre-Restoration NCC:NCP Ratio FluxFit QC');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MJuly_BMSbin.Ratio_QC)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MJuly_BMSbin.R2_QC)


% WM Ratios 
% [m_WM,b_WM,r_WM,sm_WM,sb_WM]=lsqfitgm(MJuly_BMSbin.NEP_WM_day,MJuly_BMSbin.NEC_WM_day);
% MJuly_BMSbin.Reg_Line_WM = m_WM*MJuly_BMSbin.NEP_WM_day + b_WM;
% MJuly_BMSbin.Ratio_WM = m_WM;
% MJuly_BMSbin.R2_WM = r_WM;
% % plot
% figure
% hold on; box on;
% plot(MJuly_BMSbin.NEP_WM_day,MJuly_BMSbin.NEC_WM_day,'o')
% plot(MJuly_BMSbin.NEP_WM_day,MJuly_BMSbin.Reg_Line_WM,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 July 2020 Pre-Restoration NCC:NCP Ratio');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MJuly_BMSbin.Ratio_WM)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MJuly_BMSbin.R2_WM)


[m_WM_QC,b_WM_QC,r_WM_QC,sm_WM_QC,sb_WM_QC]=lsqfitgm(MJuly_BMSbin.NEP_WM_day_QC,MJuly_BMSbin.NEC_WM_day_QC);
MJuly_BMSbin.Reg_Line_WM_QC = m_WM_QC*MJuly_BMSbin.NEP_WM_day_QC + b_WM_QC;
MJuly_BMSbin.Ratio_WM_QC = m_WM_QC;
MJuly_BMSbin.R2_WM_QC = r_WM_QC;
% plot
% figure
% hold on; box on;
% plot(MJuly_BMSbin.NEP_WM_day_QC,MJuly_BMSbin.NEC_WM_day_QC,'o')
% plot(MJuly_BMSbin.NEP_WM_day_QC,MJuly_BMSbin.Reg_Line_WM_QC,'r')
% %ylim([-50 50])
% %xlim([-25 25])
% xlabel('NCP');
% ylabel('NCC');
% title('Marker 32 July 2020 Pre-Restoration NCC:NCP Ratio WM Data Full QC');
% annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MJuly_BMSbin.Ratio_WM_QC)
% annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MJuly_BMSbin.R2_WM_QC)


%% For NEC:NEP Regressions Using Gradients
% close all 
clc

% multiply o2 gradient by -1 for O2 production
MJuly_BMSbin.dDOXY_Reg = -1.*MJuly_BMSbin.dDOXY_day_QC;
% divide TA data by 2 for alkalinity anomaly 
MJuly_BMSbin.dTA_Reg = 0.5.*MJuly_BMSbin.dTA_day_QC;

% plot to see changes - NaNs (nightime points) have already been removed
Xlength = length(MJuly_BMSbin.dDOXY_day_QC);
figure 
hold on 
DOday = plot(1:Xlength, MJuly_BMSbin.dDOXY_day_QC);
DOreg = plot(1:Xlength, MJuly_BMSbin.dDOXY_Reg);
xlabel('July Days');
ylabel('DO Gradient');
legend([DOday DOreg], {'Daytime DO','Flipped DO'}, 'location', 'northeast');
title('Marker 32 July 2020 Hourly Binned Daytime DO Gradients');

figure 
hold on 
DOday = plot(1:Xlength, MJuly_BMSbin.dTA_day_QC);
DOreg = plot(1:Xlength, MJuly_BMSbin.dTA_Reg);
xlabel('July Days');
ylabel('TA Gradient');
legend([DOday DOreg], {'Daytime TA','Regression TA'}, 'location', 'northeast');
title('Marker 32 July 2020 Hourly Binned Daytime TA Gradients');

% Regression using gradient data:
[m_G,b_G,r_G,sm_G,sb_G]=lsqfitgm(MJuly_BMSbin.dDOXY_Reg, MJuly_BMSbin.dTA_Reg);
MJuly_BMSbin.Reg_Line_G = m_G*MJuly_BMSbin.dDOXY_Reg + b_G;
MJuly_BMSbin.Ratio_G = m_G;
MJuly_BMSbin.R2_G = r_G;
% plot
figure
hold on; box on;
plot(MJuly_BMSbin.dDOXY_Reg,MJuly_BMSbin.dTA_Reg,'o')
plot(MJuly_BMSbin.dDOXY_Reg,MJuly_BMSbin.Reg_Line_G ,'r')
xlabel('NCP');
ylabel('NCC');
title('Marker 32 July 2020 NCC:NCP Ratio from Gradients');
annotation('textbox', [0.7, 0.25, 0.17, 0.1], 'String', "NCC:NCP =" + MJuly_BMSbin.Ratio_G)
annotation('textbox', [0.7, 0.15, 0.17, 0.1], 'String', "R2 =" + MJuly_BMSbin.R2_G)

clc
disp('Finished with July Metabolism Calculations');

save('MJuly20_2.mat', 'MJuly_BMS', 'MJuly_BMSbin');

%% Plot Profiles 
% plot for profile within pump heights
% close all 
% for i =1:100 %length(MJuly_ADavg.SDN)
%     figure (i)
%     scatter(MJuly_ADavg.uv(1:108,i), MJuly_ADavg.bin_depth(1:108))
%     title(['Marker 32 July Velocity Profile Number ',num2str(i),])
%     xlabel('Velocity (m/s)');
%     ylabel('Height (m)');
% end


%% Subplots 
close all
clc
MJuly_Xrange_short = [datenum('07-03-2020 18:00:00'), datenum('07-15-2020 13:00:00')];


sgtitle('Marker 32 July 2020 Results')
subplot(3,3,[1,2,3]); %Binned Gradient Plot 
hold on; box on;
DOplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dDOXY_QC, 'b-.', 'linewidth', 1.5); 
TAplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dTA_QC, 'r-.', 'linewidth', 1.5); 
% DOplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dDOXY, 'c'); 
% TAplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.dTA, 'k-'); 
plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MJuly_Xrange_short, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xticklabels({'07/03 12:00';'07/04 12:00';...
    '07/05 12:00';'07/06 12:00';'07/07 12:00';'07/08 12:00';'07/09 12:00';...
    '07/10 12:00';'07/11 12:00';'07/12 12:00';'07/13 12:00';'07/14 12:00';...
    '07/15 12:00'; '07/16 12:00'})
ylabel('\color{blue}dDO \color{black}or \color{red}dTA');
%legend([DOplot TAplot], {'DO Gradient','TA Gradient'}, 'location', 'northeast');
title('Hourly Binned Gradiets');

subplot(3,3,[4,5,6]); %Binned Flux Plot 
hold on; box on;
NEPplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_QC, 'b', 'linewidth', 1.5); 
NECplot = plot(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_QC, 'r-', 'linewidth', 1.5); 
plot(MJuly_BMSbin.SDN, zeros(size(MJuly_BMSbin.SDN)), 'k'); %Zero Line
set(gca, 'xlim', MJuly_Xrange_short, 'XTick', MJuly_tick, 'xticklabel', MJuly_tick, 'XGrid', 'on');
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('July Days');
xticklabels({'07/03 12:00';'07/04 12:00';...
    '07/05 12:00';'07/06 12:00';'07/07 12:00';'07/08 12:00';'07/09 12:00';...
    '07/10 12:00';'07/11 12:00';'07/12 12:00';'07/13 12:00';'07/14 12:00';...
    '07/15 12:00'; '07/16 12:00'})
ylabel('\color{blue}NEP \color{black}or \color{red}NEC');
%legend([NEPplot NECplot], {'NEP','NEC'}, 'location', 'southwest');
title('Hourly Binned Fluxes');

subplot(3,3,7); % Diel Composite Plot 
hold on 
nepdbin_QC = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEP_QC, 24);
necdbin_QC = parse_to_diel(MJuly_BMSbin.SDN, MJuly_BMSbin.NEC_QC, 24);
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
plot(MJuly_BMSbin.NEP_day_QC,MJuly_BMSbin.NEC_day_QC,'o')
plot(MJuly_BMSbin.NEP_day_QC,MJuly_BMSbin.Reg_Line_QC,'r')
xlabel('NEP');
ylabel('NEC');
title('NEC:NEP Ratio from Fluxes');
str1 = num2str(MJuly_BMSbin.Ratio_QC,2);
str2 = num2str(MJuly_BMSbin.R2_QC,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.412, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.412, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')


subplot(3,3,9); % Ratio Plot using gradietns 
hold on; box on;
plot(MJuly_BMSbin.dDOXY_Reg,MJuly_BMSbin.dTA_Reg,'o')
plot(MJuly_BMSbin.dDOXY_Reg,MJuly_BMSbin.Reg_Line_G ,'r')
xlabel('dDO');
ylabel('dTA');
ylim([-2 2])
title('dTA:dDO Ratio from Gradients');
str1 = num2str(MJuly_BMSbin.Ratio_G,2);
str2 = num2str(MJuly_BMSbin.R2_G,2);
%            [left to right, up, box length, box width]    
annotation('textbox', [0.693, 0.284, 0.0735, 0.03], 'String', "NEC:NEP =" + str1, 'HorizontalAlignment', 'left')
annotation('textbox', [0.693, 0.254, 0.0735, 0.03], 'String', "R^2 =" + str2, 'HorizontalAlignment', 'left')

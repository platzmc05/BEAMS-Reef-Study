% Intitial look at August SeapHOx data from Marker 32 Reef
% Michelle Platz
% 8/11/2020

close all
clc

%Parses SeapHOx data from datafile by variable 
    %37 data columns in U630 dataset
MAug20_SPraw = parse_pHOxGFdata_ARM_V3_Mar19('M_812BMS.txt');
    % First Aug Datafile: M_805BMS.txt
    % Second Aug Datafile: M_812BMS.txt
% make sure to check number of columns in the dataset and parse data
% accordingly - 3 extra columns in March dataset 

%%

%calculates O2 saturation concentration 
MAug20_SPraw.DOXY = MAug20_SPraw.O2SATPER.*calcO2sat(MAug20_SPraw.MCAT_TC, MAug20_SPraw.PSAL)./100;

 %ipHbad = SPraw.SDN > datenum('Apr-11-2018 17:45:00');
 %SPraw.pHint_prelim(ipHbad) = NaN;
 %SPraw.pHext_prelim(ipHbad) = NaN;

%calculates pH from durafet using internal reference electrode and Nernst
%equation 
MAug20_SPraw.pHint_prelim = calc_dfet_pHint(MAug20_SPraw.Vint, MAug20_SPraw.DFET_TC, -0.4);

%%
% *** end times in the daterange need to be the end of a pump cycle to avoid an
% error***
MAug20_SP = parse_to_pumpheights_ARM_2pump_Mar19(MAug20_SPraw, [datenum('08-12-2020 14:00:00'), datenum('09-1-2020 08:00:30')]);

%for Aug18: [datenum('Aug-8-2018 10:30:00'), datenum('Aug-22-2018 09:29:30')]
%for Mar19: [datenum('03-16-2019 10:21:30') datenum('03-25-2019 12:51:00')]
%for July20: U = [datenum('June-30-2020 16:00:00'), datenum('July-26-2020 09:00:00')]
            %M = [datenum('07-10-2020 16:00:00'), datenum('07-26-2020 10:29:30')]
%for Aug20: U = [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 12:00:30')]
            %M = [datenum('08-05-2020 12:00:00'), datenum('08-12-2020 09:00:30')]

MAug20_SP = calc_TA_gradientV2(MAug20_SP, 2304.87, [0.8:0.1:1.2], 1, 2);

%% Calculate Gradients 
MAug20_SP.dDOXY = MAug20_SP.DOXY(1,:) - MAug20_SP.DOXY(2,:); %Oxygen Gradient
MAug20_SP.dpH = MAug20_SP.pH(1,:) - MAug20_SP.pH(2,:); %pH Gradient 
MAug20_SP.dTA = MAug20_SP.TAtop - MAug20_SP.TAbtm(3,:); % - assuming Q=1

%% Bin Data to Hours
MAug20_Xrange = [datenum('08-12-2020 14:00:00'), datenum('08-15-2020 13:16:30')];
%hourly bins 
MAug20_SP.SDNbin = MAug20_Xrange(1):1/24:MAug20_Xrange(end);
MAug20_SP.dDOXYbin = bin_data_to_X(MAug20_SP.SDN, MAug20_SP.dDOXY, MAug20_SP.SDNbin);
MAug20_SP.dTAbin   = bin_data_to_X(MAug20_SP.SDN, MAug20_SP.dTA, MAug20_SP.SDNbin);
MAug20_SP.dpHbin   = bin_data_to_X(MAug20_SP.SDN, MAug20_SP.dpH, MAug20_SP.SDNbin);

MAug20_SP.DOXY_top_bin  = bin_data_to_X(MAug20_SP.SDN, MAug20_SP.DOXY(1,:), MAug20_SP.SDNbin);
MAug20_SP.DOXY_btm_bin  = bin_data_to_X(MAug20_SP.SDN, MAug20_SP.DOXY(2,:), MAug20_SP.SDNbin);
MAug20_SP.TA_top_bin    = bin_data_to_X(MAug20_SP.SDN, MAug20_SP.TAtop, MAug20_SP.SDNbin);
MAug20_SP.TA_btm_bin    = bin_data_to_X(MAug20_SP.SDN, MAug20_SP.TAbtm(3,:), MAug20_SP.SDNbin);
MAug20_SP.PARbin        = bin_data_to_X(MAug20_SP.SDN, MAug20_SP.PAR, MAug20_SP.SDNbin);

%% Create Datestring for Plots

% MAug20_DateString = {'08/05/2020 12:00:00';'08/06/2020 12:00:00';'08/07/2020 12:00:00';'08/08/2020 12:00:00';'08/09/2020 12:00:00';...
%     '08/10/2020 12:00:00';'08/11/2020 12:00:00'};

MAug20_DateString = {'08/12/2020 12:00:00';'08/13/2020 12:00:00';'08/14/2020 12:00:00';'08/15/2020 12:00:00';'08/16/2020 12:00:00';...
    '08/17/2020 12:00:00';'08/18/2020 12:00:00';'08/19/2020 12:00:00';'08/20/2020 12:00:00';'08/21/2020 12:00:00';'08/22/2020 12:00:00';...
    '08/23/2020 12:00:00';'08/24/2020 12:00:00';'08/25/2020 12:00:00';'08/26/2020 12:00:00';'08/27/2020 12:00:00';'08/28/2020 12:00:00';...
    '08/29/2020 12:00:00';'08/30/2020 12:00:00';'08/31/2020 12:00:00'};

    
formatIn = 'mm/dd/yyyy HH:MM:SS';
MAug20_tick = datenum(MAug20_DateString,formatIn);

%% Plot Unbinned pH Gradient

figure
hold on; box on;

plot(MAug20_SP.SDN, MAug20_SP.dpH); %pH gradient
plot(MAug20_SP.SDN, zeros(size(MAug20_SP.SDN))); %zero line

set(gca, 'xlim', MAug20_Xrange, 'XTick', MAug20_tick, 'xticklabel', MAug20_tick, 'XGrid', 'on', 'ylim', [-0.06 0.06]);
datetick('x', 'dd', 'keeplimits', 'keepticks');
xlabel('August Days');
ylabel('\Delta pH');
title('Marker 32 Aug 2020 Unbinned pH Gradient');

%% Plot Unbinned DO and TA Gradients

figure
hold on; box on;

plot(MAug20_SP.SDN, MAug20_SP.dDOXY); %oxygen gradient 
plot(MAug20_SP.SDN, MAug20_SP.dTA); %TA gradient
plot(MAug20_SP.SDN, zeros(size(MAug20_SP.SDN))); %zero line

set(gca, 'xlim', MAug20_Xrange, 'XTick', MAug20_tick, 'xticklabel', MAug20_tick, 'XGrid', 'on', 'ylim', [-20 20]);
datetick('x', 'dd', 'keeplimits', 'keepticks');
%xlim([737503.910590278 737507.743923611])
xlabel('August Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend('\DeltaO_2', '\DeltaTA', 'location', 'northeast');
title('Marker 32 Aug 2020 Unbinned DO and TA Gradients');

%% Plot Binned Gradients

%UAug20_Xrange2 = [datenum('06-30-2020 15:00:00'), datenum('07-3-2020 16:00:00')];

figure
hold on; box on;

DOplot = plot(MAug20_SP.SDNbin, MAug20_SP.dDOXYbin, 'b'); %binned DO Gradient 
TAplot = plot(MAug20_SP.SDNbin, MAug20_SP.dTAbin, 'r-'); %binned TA Gradient 
plot(MAug20_SP.SDNbin, zeros(size(MAug20_SP.SDNbin)), 'k'); %Zero Line

set(gca, 'xlim', MAug20_Xrange, 'XTick', MAug20_tick, 'xticklabel', MAug20_tick, 'XGrid', 'on', 'ylim', [-4 2]);
datetick('x', 'dd', 'keeplimits', 'keepticks');

ylim([-10 10])
%xlim([737972.7083333334 737975.4583333334])
xlabel('Aug Days');
ylabel('\DeltaO_2 or \DeltaTA');
legend([DOplot TAplot], {'\DeltaO_2','\DeltaTA'}, 'location', 'northeast');

title('Marker 32 August 2020 Hourly Binned Concentration Gradient');



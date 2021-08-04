%% Reef Study Full Datasets
% Michelle Platz - USF 
% 3/11/2021

clc
close all 
clear all

%Pull complete SP Datasets: 5 each reef
% U_630BMS
% U_805BMS
% U_812BMS
% U_902BMS
% U_928BMS
% 
% M_701BMS
% M_805BMS
% M_812BMS
% M_902BMS
% M1001BMS





%pull complete ADCP Datasets: 3 each reef - update date ranges

UJuly_AD=aquadoppraw2mat('U_6_30', 70, [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 10:15:00')]);
UAug_AD=aquadoppraw2mat('U_8_05', 70, [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 10:15:00')]);
USept_AD=aquadoppraw2mat('U_9_02', 70, [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 10:15:00')]);
MJuly_AD=aquadoppraw2mat('M_6_30', 70, [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 10:15:00')]);
MAug_AD=aquadoppraw2mat('M_8_05', 70, [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 10:15:00')]);
MSept_AD=aquadoppraw2mat('M_9_02', 70, [datenum('08-05-2020 10:00:00'), datenum('08-12-2020 10:15:00')]);





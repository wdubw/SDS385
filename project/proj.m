clc
clear
%% Read Data
% Training data: forecasts from period 1 of 104 events
p1 = xlsread('wofcdata.xlsx','Period_1');
% ID of forecastors
ID = xlsread('wofcdata.xlsx','Period_1','A2:A1234');
% Test dast: forecasts from period 2 of 90 events
p2 = xlsread('wofcdata.xlsx','Period_2');
% results of total 194 events: column 1 contains the ID's of the events and
% column 5 contains the results
result = xlsread('wofcdata.xlsx','IFP_description','A2:E195');
% for training data
rs = result(:,5); % outcomes of events

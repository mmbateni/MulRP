function [observation,stationname,nstations,threshold,years_sim] = input_paras()
%
% this function load disaggregated sub daily data from multiple files ['S',num2str(i),'disagg.mat']
% an shows statistics of those observational subdaily data
%
%% value of the following varibles must be given 
% filenamein=input('Enter an input filename (string):');
filenamein='input_data';%demo_1,demo, demo1
% threshold=input('Enter a daily precipitation threshold:'); % precipitation threshold in 'mm' 
threshold=0.5;%0.5;
% years_sim=input('Enter the number of years to generate:'); % the number of years to be generated
years_sim=500;
load(filenamein)
nstations=size(sitename,2);
% define the stations from S1 to Sn
for i=1:nstations
    stationname{i,1}=['S',num2str(i)]; 
end
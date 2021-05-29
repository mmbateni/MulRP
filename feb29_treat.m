% treat the Feb 29th
% the input data should include six columns with the order of year, month, 
% day, and precip

function [data1]=feb29_treat(data)
nn=find(data(:,2)==2 & data(:,3)==29);
% evenly distribute the precip in Feb 29th to Feb 28 and Mar 1st
% the temperature in Feb 29th is simply removed
for i=1:length(nn)
    data(nn(i)-1,4)=data(nn(i)-1,4)+data(nn(i),4)/2;
    data(nn(i)+1,4)=data(nn(i)+1,4)+data(nn(i),4)/2;
end
data(nn,:)=[];
data1=data;

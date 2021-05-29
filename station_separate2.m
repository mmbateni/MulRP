function [generationOcc]=station_separate2(data,years_sim,nstations)
% data is a structure of 12*1, which each elemet is a matrix of
% (length_month*n_simyear)*n_stations
%'generationOcc' has  long elements for each station
%% combine seasonal precip
days=years_sim*365; % the number of days to be generated
multiexp_pre=zeros(days,nstations);
multiexp_pre=addymd(multiexp_pre);
%
nn=find(multiexp_pre(:,2)==12); % dec
multiexp_pre(nn,4:end)=data(12).occ;
nn=find(multiexp_pre(:,2)==1); % jan
multiexp_pre(nn,4:end)=data(1).occ;
nn=find(multiexp_pre(:,2)==2); % feb
multiexp_pre(nn,4:end)=data(2).occ;
nn=find(multiexp_pre(:,2)==3); % mar
multiexp_pre(nn,4:end)=data(3).occ;
nn=find(multiexp_pre(:,2)==4); % apr
multiexp_pre(nn,4:end)=data(4).occ;
nn=find(multiexp_pre(:,2)==5); % may
multiexp_pre(nn,4:end)=data(5).occ;
nn=find(multiexp_pre(:,2)==6); % jan
multiexp_pre(nn,4:end)=data(6).occ;
nn=find(multiexp_pre(:,2)==7); % feb
multiexp_pre(nn,4:end)=data(7).occ;
nn=find(multiexp_pre(:,2)==8); % mar
multiexp_pre(nn,4:end)=data(8).occ;
nn=find(multiexp_pre(:,2)==9); % apr
multiexp_pre(nn,4:end)=data(9).occ;
nn=find(multiexp_pre(:,2)==10); % may
multiexp_pre(nn,4:end)=data(10).occ;
nn=find(multiexp_pre(:,2)==11); % may
multiexp_pre(nn,4:end)=data(11).occ;
%
%% put precip in a file for each station
for i=1:nstations
    S=zeros(days,3);
    S=addymd(S);
    S(:,4)=multiexp_pre(:,i+3); % generated precip event
    eval(['generationOcc.S',num2str(i),'=S;'])
end
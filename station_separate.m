function [generation]=station_separate(data,years_sim,nstations,length_month)
% separate the generated precipitation and temperature by station
% data is a structure of 4*1, which each elemet is a matrix of
% (length_season*n_simyear)*n_stations
%'generation' has  a long elements for each station
%% combine seasonal precip
days=years_sim*365; % the number of days to be generated
multiexp_pre=zeros(days,nstations);
multiexp_pre=addymd(multiexp_pre);
%
nn=find(multiexp_pre(:,2)==12); % dec
multiexp_pre(nn,4:end)=data(1,1).pre(1:years_sim*length_month(1,12),:);
nn=find(multiexp_pre(:,2)==1); % jan
multiexp_pre(nn,4:end)=data(1,1).pre(years_sim*length_month(1,12)+1 ...
    :years_sim*length_month(1,12)+years_sim*length_month(1,1),:);
nn=find(multiexp_pre(:,2)==2); % feb
multiexp_pre(nn,4:end)=data(1,1).pre(years_sim*length_month(1,12) ...
    +years_sim*length_month(1,1)+1:end,:);
%
nn=find(multiexp_pre(:,2)==3); % mar
multiexp_pre(nn,4:end)=data(1,2).pre(1:years_sim*length_month(1,3),:);
nn=find(multiexp_pre(:,2)==4); % apr
multiexp_pre(nn,4:end)=data(1,2).pre(years_sim*length_month(1,3)+1 ...
    :years_sim*length_month(1,3)+years_sim*length_month(1,4),:);
nn=find(multiexp_pre(:,2)==5); % may
multiexp_pre(nn,4:end)=data(1,2).pre(years_sim*length_month(1,3) ...
    +years_sim*length_month(1,4)+1:end,:);
%
nn=find(multiexp_pre(:,2)==6); % jun
multiexp_pre(nn,4:end)=data(1,3).pre(1:years_sim*length_month(1,6),:);
nn=find(multiexp_pre(:,2)==7); % jul
multiexp_pre(nn,4:end)=data(1,3).pre(years_sim*length_month(1,6)+1 ...
    :years_sim*length_month(1,6)+years_sim*length_month(1,7),:);
nn=find(multiexp_pre(:,2)==8); % aug
multiexp_pre(nn,4:end)=data(1,3).pre(years_sim*length_month(1,6) ...
    +years_sim*length_month(1,7)+1:end,:);
%
nn=find(multiexp_pre(:,2)==9); % sep
multiexp_pre(nn,4:end)=data(1,4).pre(1:years_sim*length_month(1,9),:);
nn=find(multiexp_pre(:,2)==10); % oct
multiexp_pre(nn,4:end)=data(1,4).pre(years_sim*length_month(1,9)+1 ...
    :years_sim*length_month(1,9)+years_sim*length_month(1,10),:);
nn=find(multiexp_pre(:,2)==11); %nor
multiexp_pre(nn,4:end)=data(1,4).pre(years_sim*length_month(1,9) ...
    +years_sim*length_month(1,10)+1:end,:);
multiexp_pre(:,1:3)=[];
%
%% put precip in a file for each station
for i=1:nstations
    S=zeros(days,3);
    S=addymd(S);
    S(:,4)=multiexp_pre(:,i); % generated precip
    eval(['generation.S',num2str(i),'=S;'])
end

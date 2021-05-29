% This function establishes link between occurrence index and average precip 
% amounts (and stdev when using gamma distribution) for each station and 
% construct the multi-exponential (or multi-gamma) distribution for 
% each station 
function [phat_multiexp,dist,se]=multi_distribution_par(amounts,stationname,...
    begin_year,end_year,occurrences_gen,nstations,seasons,begin_season,length_season,...
    threshold,Tolerance,maxiter)%observation,years_sim,dist, input:occ_index
warning('off','all');
minyear=min(begin_year);
maxyear=max(end_year);
raincum=[];
% generate season precip occrrence (0 & 1s) Dim=(length_month*years_sim)*nstations (the spatial correlation is induced)
occurencesDJF=[occurrences_gen(1,12).occ; occurrences_gen(1,1).occ; occurrences_gen(1,2).occ];
occurencesMAM=[occurrences_gen(1,3).occ; occurrences_gen(1,4).occ; occurrences_gen(1,5).occ];
occurencesJJA=[occurrences_gen(1,6).occ; occurrences_gen(1,7).occ; occurrences_gen(1,8).occ];
occurencesSON=[occurrences_gen(1,9).occ; occurrences_gen(1,10).occ; occurrences_gen(1,11).occ];
%
wDJF=corrcoef(occurencesDJF);
wMAM=corrcoef(occurencesMAM);
wJJA=corrcoef(occurencesJJA);
wSON=corrcoef(occurencesSON);
se=zeros(nstations,nstations,4);
figure('Name',['Historic.CDF.amounts'],'NumberTitle','off');
ind_plot=0;
% % preliminary treatment create matrices 'shortname'occ for occurences, with
% % NaN as missing data and 'shortname'amounts for amounts with NaN as
% % missing data, all using a threshold sepcified above 
% for i=1:nstations
%     eval(['mat=' char(stationname(i)) ';']);
%     amounts=mat;
%     occ=zeros(size(mat));
%     occ(amounts>threshold)=1;
%     amounts(mat<=threshold)=0;
%     eval([char(stationname(i)) 'amounts=amounts;']);
%     eval([char(stationname(i)) 'occ=occ;']);    
% end
% %
%% rearrange all data, create precipitation matrices of equal size for all
% stations containing NaN for missing values
%

for i=1:nstations
    eval(['amount=' 'amounts.',char(stationname(i)),';' ]);  
    firstrows=begin_year(i)-minyear;
    lastrows=maxyear-end_year(i);
    amnts=[NaN*ones(firstrows,365); amount; NaN*ones(lastrows,365)];
    eval([char(stationname(i)) 'amounts=amnts;' ]);
end
for ii=1:nstations  % put december in first
    eval(['SS',num2str(ii),'amounts=S',num2str(ii),'amounts;']);
    eval(['SS',num2str(ii),'amounts(:,1:31)=S',num2str(ii),'amounts(:,335:365);']);
    eval(['SS',num2str(ii),'amounts(:,32:365)=S',num2str(ii),'amounts(:,1:334);']);
    eval(['S',num2str(ii),'amounts=SS',num2str(ii),'amounts;']);
end
%%%%	
for ijk=1:nstations
	  station_num=ijk;  % pick station to analyze       
    for qq=1:4   % loop for each season
        [iii,~]=size(S1amounts);
        mm=iii*length_season(qq);
        precipdata=zeros(mm,nstations);
	
        % select the month block for each station, rearrange it as a column and
        % produce a matrix containing all precip data for each station and
        % each season
        for i=1:nstations
            eval(['m1=' char(stationname(i)) 'amounts(:,begin_season(qq):begin_season(qq)+length_season(qq)-1);' ]);
            [c1,c2]=size(m1);
            c=c1*c2;
            precipdata(:,i)=reshape(m1,c1*c2,1);
        end		
		    % explore the relationship of mean precip for all data for station 1
		    occ_coef=-999*ones(c,1);
        % weight matrix
        if qq==1
            tempo=wDJF;
        elseif qq==2
            tempo=wMAM;
        elseif qq==3
            tempo=wJJA;
        else
            tempo=wSON;
        end       
        for i=1:c   % for each day
            if (precipdata(i,station_num))>0
                k=find(precipdata(i,:)>=0);  % find existing data
                kk=precipdata(i,:)>0;
                occ=zeros(1,nstations);
                occ(kk)=1;
                occ=occ(k);  % only take existing data             

                 tempo2=tempo(station_num,:);
                 occ_coef(i)=occ*tempo2(k)'/sum(tempo2(k));
            end
        end
        ind=precipdata(:,station_num)>0;%& moran>=-1
        raincum=precipdata(ind,station_num)';
        wer=precipdata>0;
        precipdata_season=precipdata(wer);
%% *****************************************************************
    % finding best fitted Extended Generalized Pareto distribution  
    % to the rainfall values for each station pooled in four seasons  
##    opt = optimset('GradObj','off',...
##   'Hessian','off',...
##   'LargeScale','off',...
##   'Display','off',...
##   'MaxIter', maxiter,...
##   'MaxFunEvals', 1000000, ...
##   'TolFun', Tolerance ,...
##   'TolX', 1.000e-006,...
##   'FunValCheck','off',...
##   'DerivativeCheck','off',...
##   'Diagnostics','off',...
##   'Algorithm','interior-point');
%% *****************************************************************
% starting values for estimation of precipitation depth distribution
% at each station for each season
        [DistName,ParameterValues,~]=fitdistall(raincum',threshold,maxiter,Tolerance,precipdata_season);
        phat_value=ParameterValues;
        dist_fit=DistName;
%       parinit=[1,2.5,0.08];
%       [phat_value,~,~,~,~,~,hessian] = fmincon(@(parinit) egp_ll(raincum,threshold,parinit),parinit,[0,0,0;0,-1,0;0,0,-1],[0;-0.5;0],[],[],[],[],[],opt);
%       se(:,ijk,qq) = sqrt(diag(inv(hessian)));%se=1;
        eval(['raincum' char(seasons(qq)) '=raincum;' ]);
%       histogram(raincum);
%       [N,EDGES] = histcounts(raincum);
        [N,Cbins] = hist(raincum);
        N1=cumsum(N)/length(raincum);
%       B=mean([EDGES(1:end-1);EDGES(2:end)]);
        B=Cbins;
        ind_plot=ind_plot+1;
        subplot(nstations,4,ind_plot);
        [xb,yb]=stairs(B,N1);
        plot(xb,yb);
        hold on;
        y = linspace(0.05,0.95);
        x=distinv(y,dist_fit,phat_value,0);%pd OR DistName,ParameterValues
%       x=egpinv(y,phat_value,0);%0 could be threshold
        plot(x,y,'k');
        title([char(seasons(qq)) char(stationname(station_num))]);
% *****************************************************************         
        eval(['phat_value' char(seasons(qq)) '=phat_value;' ]);
        eval(['dist' char(seasons(qq)) '=dist_fit;' ]);
        eval(['Amounts_nzero.' char(stationname(station_num)) char(seasons(qq)) '=raincum;']);
    end
        eval(['phat_multiexp.',char(stationname(station_num)),'.phatDJF=phat_valueDJF;']);%'DistParams' ,char(stationname(station_num)),char(seasons(qq))
        eval(['phat_multiexp.',char(stationname(station_num)),'.phatMAM=phat_valueMAM;']);
        eval(['phat_multiexp.',char(stationname(station_num)),'.phatJJA=phat_valueJJA;']);
        eval(['phat_multiexp.',char(stationname(station_num)),'.phatSON=phat_valueSON;']);
        eval(['dist.',char(stationname(station_num)),'.DistDJF=distDJF;']);%'DistParams' ,char(stationname(station_num)),char(seasons(qq))
        eval(['dist.',char(stationname(station_num)),'.DistMAM=distMAM;']);
        eval(['dist.',char(stationname(station_num)),'.DistJJA=distJJA;']);
        eval(['dist.',char(stationname(station_num)),'.DistSON=distSON;']);
end
% assignin('base','Amounts_nzero',Amounts_nzero);





function [daily_obs,X_total,n_subdaily,XB] = descrip_stats(stationname,nstations,threshold)
%
% this function load disaggregated sub daily data from multiple files ['S',num2str(i),'disagg.mat']
% an shows statistics of those observational subdaily data
%%
daily_obs=cell(1,nstations);
X_total=cell(1,nstations);
n_subdaily=zeros(1,nstations);
Stats= {'Skewness_nonzero';'Skewness_all';'Mean_nonzero';'Mean_all';'STD_nonzero';'STD_all'};
XB=cell(1,nstations);
for i=1:nstations
    Xt=eval(['csvread(' '"S',num2str(i),'disagg.csv")']);%['load ' 'S',num2str(i),'disagg_1.mat/csv']
    Xt(any(isnan(Xt), 2), :) = [];
    daily_obs{i}=sum(Xt,2);
    X_total{i}=Xt;
    B = Xt ;
for ik = 1:size(Xt,1)
    abs_sum=norm(Xt(ik,:),1);
    if abs_sum>0
       B(ik,:) = Xt(ik,:)./abs_sum ;
    end
end
XB{i}=B;
%% Statstics of Historical values
figure('Name',['Historic6hr.Nonzero.' char(stationname(i))],'NumberTitle','off');
n_subdaily(i)=size(Xt,2);
for ij=1:n_subdaily(i)
    subplot(1,n_subdaily(i),ij);
    Xt_dis=Xt(:,ij);
    Xt_dis_nz = Xt_dis(Xt_dis>0);%>threshold
    boxplot(Xt_dis_nz);
%   boxplot(Xt_dis)
    skew_nz(ij)=skewness(Xt_dis_nz);
    skew_tot(ij)=skewness(Xt_dis);
    mean_nz(ij)=mean(Xt_dis_nz);
    mean_tot(ij)=mean(Xt_dis);
    std_nz(ij)=std(Xt_dis_nz);
    std_tot(ij)=std(Xt_dis);
end
data_table=[skew_nz;skew_tot;mean_nz;mean_tot;std_nz;std_tot];
#T = table(data_table,'RowNames',Stats);
T = table(Stats,data_table);
title=strcat(['Historic.subdaily.' char(stationname(i))]);
disp (title)
prettyprint (T)
#figure('Name',['Historic.subdaily.' char(stationname(i))],'NumberTitle','off');
#uitable('Data',T{:,:},'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
end
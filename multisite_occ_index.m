function [seasonal_cor,corr_amounts_rand_index,corr_amounts_gen_index,...
    precip_gen_index]=multisite_occ_index(data,phat_multiexp,dist,corr_amounts,...
    stationname,occurrences_gen,nstations,season,begin_season,length_season,...
    threshold,years_sim,Tolerance,maxiter,begin_year,end_year,Amounts_nzero,seed)%site_bins,dist
%% generate precipitation amounts based on the occurrence index
% calculate seasonal correlation for precip amount
[corr_amountsS]=seasonal_corr(data,stationname,nstations,threshold,...
                begin_season,length_season,season,begin_year,end_year);
f1=figure('Name',['Simulated.CDF.amounts'],'NumberTitle','off');
f2=figure('Name',['QQPlot.YSimulated.VS.XHistoric'],'NumberTitle','off');
ind_plot1=0;
ind_plot2=0;
randn ("seed",seed)% for Octave
for ijk=1:4
    C=corr_amountsS(1,ijk).cor;  
	% obtain correlation matrix of generated precip occurrence 
    if ijk == 1 
	    occurrences=[occurrences_gen(1,12).occ; occurrences_gen(1,1).occ; occurrences_gen(1,2).occ];
    end
    
    if ijk == 2  
	    occurrences=[occurrences_gen(1,3).occ; occurrences_gen(1,4).occ; occurrences_gen(1,5).occ];
    end
    
    if ijk == 3  
	    occurrences=[occurrences_gen(1,6).occ; occurrences_gen(1,7).occ; occurrences_gen(1,8).occ];
    end
    
    if ijk == 4  
	    occurrences=[occurrences_gen(1,9).occ; occurrences_gen(1,10).occ; occurrences_gen(1,11).occ]; 
    end
    seasonal_cor(ijk).cor=C;
    seasonal_cor(ijk).season=char(season(ijk));  
%   generate random number		
    random=zeros(nstations,length_season(ijk)*years_sim);
    for i=1:nstations
        random(i,:)=randn(1,length_season(ijk)*years_sim);
    end
%   automatic determination of correlation matrix of random numbers 
%   tttt=cputime;
    Char_season=char(season(ijk));
    n=length_season(ijk)*years_sim;
    [M,K,precip]= spatial_iterate_amo_index(C,phat_multiexp,dist,nstations,random,...
     n,threshold,occurrences,Char_season,Tolerance,maxiter);
%   time_in_second=cputime-tttt;
%   correlation of random numbers needed
    corr_amounts_rand_index(ijk).cor=M;
    corr_amounts_rand_index(ijk).season=char(season(ijk));
%   resulting correlation of amount
    corr_amounts_gen_index(ijk).cor=K;
    corr_amounts_gen_index(ijk).season=char(season(ijk));
%   generated precip amount
    precip_gen_index(ijk).pre=precip;
    precip_gen_index(ijk).season=char(season(ijk)); 
%
%% produce graphics
%
    for u=1:nstations
        precip_station=precip(:,u);
        [N,Cbins] = hist(precip_station);
%       [N,EDGES] = histcounts(precip_station);
        N1=cumsum(N)./length(precip_station);
%       B=mean([EDGES(1:end-1);EDGES(2:end)]);
        B=Cbins;
        ind_plot1=ind_plot1+1;
        figure(f1);
        subplot(4,nstations,ind_plot1);
        stairs(B,N1);
        title([char(season(ijk)) char(stationname(u))]);
        eval(['Simulated_Amounts_nzero.',char(stationname(u)),char(season(ijk)),'=precip_station;']);
    end
    for v=1:nstations
        str1=[char(stationname(v)),char(season(ijk))];
        str2=['Amounts_nzero.' str1];
        x_obs=eval(str2);
        y_sim=precip(:,v);
        y_sim_s=zeros(length(x_obs),1);
        for ki=1:100
            y_sim_short=y_sim(randperm(length(y_sim),length(x_obs)));
            y_sim_short=sort(y_sim_short);
            y_sim_s=y_sim_s+y_sim_short;
        end
        y_sim=(y_sim_s/100)';
        ind_plot2=ind_plot2+1;
        figure(f2);
        subplot(4,nstations,ind_plot2);
        qqplot(x_obs,y_sim);
        axis([0 100 0 100])
        title([char(season(ijk)) char(stationname(v))]);
    end
end   
%
%% produce graphics
%
% correlation of generated precip amount    
figure('Name','amount correlation of historical and generated precip','NumberTitle','off');
for ii=1:4
    C=seasonal_cor(1,ii).cor;
    CC=reshape(triu(C,1),1,[]);
    i=CC~=0;
    CC=CC(i);
    K=corr_amounts_gen_index(1,ii).cor;
    KK=reshape(triu(K,1),1,[]);
    i=KK~=0;
    KK=KK(i);
    subplot(2,2,ii)
    plot(CC,KK,'o',[0 1],[0 1]);
    title(['all stations for season ' char(season(ii))]);
    xlabel('observed correlation');
    ylabel('generated correlation');    
    axis([0 1 0 1])
    axis square
    set(gca,'FontSize',12)
    set(gcf,'Color',[1 1 1])
end  
figure('Name','Cormatrix Historical Data Values ','NumberTitle','off');
for ij=1:4
    obs_Corr_amount=corr_amountsS(1,ij).cor; 
    subplot(2,2,ij);
    plot_corrmat(obs_Corr_amount,['Corr.amount.Seasonal.' char(season(ij))])
end
figure('Name','Cormatrix Simulated Data Values ','NumberTitle','off');
for ij=1:4
    sim_Corr_amount=corr_amounts_gen_index(ij).cor;
    subplot(2,2,ij);
    plot_corrmat(sim_Corr_amount,['corr.Amounts(nzero).' char(season(ij))])
end
end
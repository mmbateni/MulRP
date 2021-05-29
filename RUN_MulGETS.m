
% *************************************************************************
% joint Return Period calculation using Multisite daily precipitation generator (MulRP)
% *************************************************************************
%
% MulRP is OCTAVE-based multisite data generator for generating
% spatially correlated daily precipitation series. The algorithm originates from the 
% Wilks approach proposed in 1998. In 2007, Francois Brissette et al. 
% presented an algorithm for efficient generation of multisite
% precipitation data following the Wilks approach. In this work,
% the distribution representing the precipitation intensity at each location,
% is the E-GPD distribution. This distribution was first proposed by Papastathopoulos and Tawn (2013),
% who referred to it as an extended GP-Type III distribution, 
% and it has since been shown to adequately model the whole range of 
% precipitation intensities (Naveau et al., 2016).
% Compared to other heavy-tailed distributions applied to daily precipitation amounts
% (e.g., mixtures of GPD and gamma distribution; see Vrac and Naveau, 2007),
% the E-GPD is parsimonious and provides a very good compromise 
% between flexibility and stability, which is an essential feature for extrapolation.
%This distribution can be described by a smooth transition between a gamma-like
% distribution and a heavy-tailed Generalized Pareto distribution (GPD).
% Here, we used E-GPD with (Naveau et al,2016)Model (1) with G(v)=v^k. for
% that type of GPD, ? > 0 is a scale parameter, and ? > 0. 
% In this work, the estimation of the ? parameter is boundedbelow by 0. 
% When ? < 0, the E-GPD distribution has an upper bound. 
% As shown by many recent studies (e.g., Serinaldi and Kilsby, 2014a), 
% negative estimates of ? are usually due to parameter uncertainty and are not realistic. 
% The two remaining parameters of the E-GPD, the scale parameter ? 
% and the parameter of the transformation ?, are estimated from 
% the observations available at that station.  As local estimations of the GPD tail
% exhibit a lack of robustness, we propose estimating the ? parameter of the E-GPD
% using data from all stations.

% 
% ****************************
% Input data
% ****************************
% The input data consists of daily precipitation for multisite,
% precipitation data shall be separated by stations with a matlab
% structure named "observation", the data of each station shall be provided
% with the order of year, month, day and precipitation
% Missing data should be assigned as NaN.
%
% ****************************
% Output data
% ****************************
% The output consists of daily precipitation, Tmax and Tmin, the generated 
% meteorological time series is separated by stations with a matlab structure
% named "generation", the order of each station is year, month, day  
%  and precipitation
% 
% ****************************
% Generation procedure
% **************************** 
% the generation of precipitation includes five steps
% step 1: determination of weather generator parameters: p00, p10 and 
%         precip distribution function on a monthly scale
% step 2: computation of correlation matrices of precipatation occurrence 
%         and amounts 
% step 3: generate spatially correlated precipitation occurrence
% step 4: establish link between occurrence index and average 
%         precip amounts for each station and construct the multi-exponential
%         or multi-gamma distribution for each station.generate 
%         precipitation amounts based on the occurrence index of generated occurrence
% step 5: Rearranging output
% 

%
% References:
% (1) Wilks, D. S., 1998. Multisite generalization of a daily stochastic
% precipitation generation model. J. Hydrol. 210, 178-191.
% (2) Brissette, F.P., Khalili, M., Leconte, R., 2007. Efficient stochastic
% generation of multi-site sythetic precipitation data. J. Hydrol. 345,
% 121-133.
% (3) Naveau, P., R. Huser, P. Ribereau, and A. Hannart (2016), Modeling jointly low,
% moderate, and heavy rainfall intensities without a threshold selection,
% Water Resour. Res., 52, 2753-2769, doi:10.1002/2015WR018552

%% For setting up OCTAVE
% install required packages
##pkg install https://github.com/apjanke/octave-tablicious/releases/download/v0.3.5/tablicious-0.3.5.tar.gz
##pkg install -forge struct
##pkg install -forge optim
##pkg install -forge statistics
##pkg install -forge nan
% load required packages
pkg load tablicious
pkg load optim
pkg load statistics
pkg load nan
%
clear all
close all
display('**********************************************************START********************************************************')
display('*******************************************************basic input*****************************************************')

%% basic inputs and declare important variables
% ******************************************************************************
filenameout='out';
% filenameout=input('Enter an output filename (string):');
% if graph==1, the correlation matrices of generated data will be plotted 
% over the obseved coresponding; the seasonal dependency of mean (and 
% standard deviation of) precipitation and occurrence index will also be
% plotted for each station
% ******************************************************************************
% declare several important variables as inputs
[observation,stationname,nstations,threshold,years_sim]=input_paras();
months={'jan' ;'feb'; 'mar'; 'apr';'may';'jun';'jul';'aug';'sep';'oct';'nov';'dec'};
length_month=[31 28 31 30 31 30 31 31 30 31 30 31]; % days of each month
begin_month=[1 32 60 91 121 152 182 213 244 274 305 335]; % the first julin day of each month
seasons={'DJF' ;'MAM'; 'JJA'; 'SON';};
begin_season=[1 91 183 275];
length_season=[90 92 92 91];
Tolerance=0.001;% convergence threshold
maxiter=500; % convergence iteration
seed=234;
%% basic processing of inputs
% find the start year and the end year for each station
begin_year=zeros(1,nstations);
end_year=begin_year;
for i=1:nstations
    begin_year(1,i)=eval(['observation.S',num2str(i),'(1,1);']); % begin year of the dataset
    end_year(1,i)=eval(['observation.S',num2str(i),'(end,1);']); % end year of the dataset
end
%
% extract precip for each station
for i=1:nstations
    S=eval(['observation.S',num2str(i),';']); 
    [data1]=feb29_treat(S); % treat Feb 29th, the model does not account for Feb 29th
    pre=data1(:,4);
    pre=reshape(pre,365,[]);
    pre=pre';
    eval(['data.',char(stationname(i)) '=pre;']);
end
%% ************************************************************************
%         START MULTISITE PRECIPITATION GENERATION (SIX STEPS)
% *************************************************************************
%% step 1: determination of transition probabilities: p00, p10  
% and precip distribution function on a monthly scale
[trans_prob,amounts,occs]=trans_proba(data,stationname,nstations,threshold,...
    begin_month,length_month,months);
%% step 2: computation of correlation matrices of precipatation occurrence 
% and amounts 
[corr_occ,corr_amounts]=corr_precip(occs,amounts,stationname,nstations,...
    begin_month,length_month,months,begin_year,end_year);
%% step 3: generate spatially correlated precipitation occurrence
% the step 3 includes two sub-steps: 
% (1) automatic determination of correlation matrix of random number using
% algorithm .Diagonalization and repalcement of negative eigenvalues
% if matrix is non-positive definite 
% (2) generate precipitation occurrence using first-order Markov chain with
% Cholesky factorization 
[corr_occ_rand,corr_occ_gen,occurrences_gen,normal_rand]=multisite_occ_generation(corr_occ,...
    trans_prob,stationname,nstations,length_month,years_sim,months,Tolerance,maxiter,seed);
%% step 4: establish link between occurrence index (eq.10) and average 
% precipitation amounts for each station and construct the amount
% distribution for each station (eq.11)
% generate precipitation amounts based on the occurrence index and save the results
[phat_multiexp,dist,se]=multi_distribution_par(amounts,stationname,begin_year,end_year,...
    occurrences_gen,nstations,seasons,begin_season,length_season,threshold,Tolerance,maxiter);
[seasonal_cor,corr_amounts_rand_index,corr_amounts_gen_index,precip_gen_index]=multisite_occ_index(...
    data,phat_multiexp,dist,corr_amounts,stationname,occurrences_gen,nstations,seasons,...
    begin_season,length_season,threshold,years_sim,Tolerance,maxiter,begin_year,end_year,Amounts_nzero,seed);
%% step 5: Rearranging Data
% separate the generated precipitation by month and drawing histograms
    [amounts_gen,amounts_gen_imonth]=month_separate(precip_gen_index,years_sim,nstations,length_month);
% separate the generated precipitation by stations
    [generation]=station_separate(precip_gen_index,years_sim,nstations,length_month);
    [generationOcc]=station_separate2(occurrences_gen,years_sim,nstations);
    occ_generation=cell2mat(struct2cell(generationOcc));
    pre_generation=cell2mat(struct2cell(generation));
    out_event=reshape(occ_generation(:,4),[],nstations);
    out=reshape(pre_generation(:,4),[],nstations);
% save the generated final results
    csvwrite('generated_seq.csv',out)
% *************************************************************************
%                   END MULTISITE PRECIPITATION GENERATION
% *************************************************************************
% *************************************************************************
% *************************************************************************
%

% *************************************************************************
% Non-parametric Temporal Disaggregation of precipitation data
% *************************************************************************
%
% Here, we disaggregate daily precipitation data into the 4 values of
% 6-hours for each day. A modified K?nearest neighbor (K?NN)?based method
% proposed by Nowak et al, 2010 is used. Morever, a Three-day Rainfall Pattern 
% classification is adopted. herefore, when daily rainfall data is disaggregated
% into sub-daily rainfall, statistical properties such as mean, standard deviation, etc.,
% should be preserved. However, the modi?ed KNNR in Nowak et al, 2010 is not applicable
% to disaggregate daily rainfall into sub-daily because it is not guaranteed that the
% statistical properties will be reproduced. Therefore, sub-daily rainfall distribution
% or some information about sub-daily rainfall depth are required.the seven rainfall patters 
% are presented. When the amount of rainfall in D day is in the middle, 
% rainfall patterns “decrease” from previous day (D1) to subsequent day (D + 1)
% in Types 1, 3, and 4. Rainfall patterns “increase” from previous day (D1) 
% to subsequent day (D + 1) in types 2, 5, and 6. Lastly, Type 7 shows the “zero” amount of rainfalls
% in previous (D1) and subsequent (D + 1) days.
% When the hourly rainfall is retransformed, small number of negative values could
% be generated. In this study, the sum of the negative values was distributed over
% the positive hourly rainfall data according to the ratio of rainfall depth. Therefore,
% in?uence of negative value was removed.
% %
% References:
% (1)Nowak, K., Prairie, J., Rajagopalan, B., & Lall, U. (2010). A nonparametric
% stochastic approach for multisite disaggregation of annual to daily streamflow.
% Water Resources Research, 46(8).
% (2) Park, H., & Chung, G. (2020). A Nonparametric Stochastic Approach for
% Disaggregation of Daily to Hourly Rainfall Using 3-Day Rainfall Patterns.
% Water, 12(8), 2306.

%% ************************************************************************
%         START DISAGGREGATION OF  PRECIPITATION DATA (FOUR STEPS)
% *************************************************************************
clear all
close all
%% step 1: LOADING ORIGINAL SUB-DAILY DATA 
[~,stationname,nstations,threshold,~]=input_paras();
[daily_obs,X_total,n_sub,obs_normalized] = descrip_stats(stationname,nstations,threshold);
ind_daily_obs= cellfun(@(x) x>threshold,daily_obs,'UniformOutput', false);
%% step 2: LOADING SIMULATED (AGGREGATED) DAILY DATA 
len_Zsim_i=zeros(1,nstations);
out_sim=csvread('generated_seq.csv');
 for i=1:nstations
     Zsim_i=(out_sim(:,i))';
     Zsim{i}=Zsim_i;
     len_Zsim_i(i)=length(Zsim_i);
     diff_Zsim_i=diff(Zsim_i);
     diff_Zsim{i}=diff_Zsim_i;
 end  
 %% step 2: IDENTIFICATION OF Three-day Rainfall Patterns
 for i=1:nstations
    ind_Zsim{i}=(pattern3((Zsim{1,i})));
    ind_daily{i}=(pattern3((daily_obs{1,i})));
 end
 %% step 3: DISAGGREGATING DATA OF EACH STATION
Stats= {'Skewness_nonzero';'Skewness_all';'Mean_nonzero';'Mean_all';'STD_nonzero';'STD_all'};
historic_6h=cell(1,nstations);
sim_6h=cell(1,nstations);
 for i=1:nstations
            Zsim_disagg=zeros(n_sub(i),len_Zsim_i(i));
            Zobs_disagg=zeros(n_sub(i),length((daily_obs{1,i})));
            ZsimBin_disagg=zeros(n_sub(i),len_Zsim_i(i));
            ZobsBin_disagg=zeros(n_sub(i),length((daily_obs{1,i})));
           
       for ik=0:7
           ind_pattern1=find(ind_daily{i}==ik);
           ind_pattern2=find(ind_Zsim{1,i}==ik);
           Xt=(X_total{1,i});
           daily_obs_t=(daily_obs{1,i});
%          ind_daily_obs_t=cell2mat(ind_daily_obs{i});
           XNormt=(obs_normalized{1,i});
           Xtu=(Xt(ind_pattern1,:))';
%
           XNormtu=(XNormt(ind_pattern1,:))';
           XNormtu_mean=mean(XNormtu,2);
%
           daily_obs_tu=daily_obs_t(ind_pattern1); 
%          ind_daily_obs_tu=ind_daily_obs_t(ind_pattern1);
           Zsimu=Zsim{1,i}(ind_pattern2);         
           [x_s] = disagg_knn(Xtu,Zsimu);
           [x_h] = disagg_knn(Xtu,daily_obs_tu);
           Zsim_disagg(:,ind_pattern2)=x_s;
           Zobs_disagg(:,ind_pattern1)=x_h;
           ZsimBin_disagg(:,ind_pattern2)=repmat(XNormtu_mean,1,length(ind_pattern2));
           ZobsBin_disagg(:,ind_pattern1)=repmat(XNormtu_mean,1,length(ind_pattern1));
       end
       Zsim_disagg=Zsim_disagg';
       Zobs_disagg=Zobs_disagg';
       ZsimBin_disagg=ZsimBin_disagg';
       ZobsBin_disagg=ZobsBin_disagg';
       sim_6h{1,i}=Zsim_disagg;
       historic_6h{1,i}=Zobs_disagg;
       figure('Name',['Disaggregated6hr.Nonzero.' char(stationname(i))],'NumberTitle','off');
       for kj=1:n_sub(i)
       subplot(1,n_sub(i),kj);
       Zsim_dis=Zsim_disagg(:,kj);
       Zsim_dis_nz = Zsim_dis(Zsim_dis>0);
       boxplot(Zsim_dis_nz);
       skew_nz_gen(kj)=skewness(Zsim_dis_nz);
       skew_tot_gen(kj)=skewness(Zsim_dis);
       mean_nz_gen(kj)=mean(Zsim_dis_nz);
       mean_tot_gen(kj)=mean(Zsim_dis);
       std_nz_gen(kj)=std(Zsim_dis_nz);
       std_tot_gen(kj)=std(Zsim_dis);
%      boxplot(Zsim_dis)
       end
       figure('Name',['DisaggregatedHistoric6hr.Nonzero.' char(stationname(i))],'NumberTitle','off');
       for kk=1:n_sub(i)
       subplot(1,n_sub(i),kk);
       Zobs_dis=Zobs_disagg(:,kk);
       Zobz_dis_nz = Zobs_dis(Zobs_dis>0);
       boxplot(Zobz_dis_nz);
       skew_nz_obs(kk)=skewness(Zobz_dis_nz);
       skew_tot_obs(kk)=skewness(Zobs_dis);
       mean_nz_obs(kk)=mean(Zobz_dis_nz);
       mean_tot_obs(kk)=mean(Zobs_dis);
       std_nz_obs(kk)=std(Zobz_dis_nz);
       std_tot_obs(kk)=std(Zobs_dis);
%      boxplot(Zsim_dis)
       end
    data_table_sim=[skew_nz_gen;skew_tot_gen;mean_nz_gen;mean_tot_gen;std_nz_gen;std_tot_gen];
%    T = table(data_table_sim,'RowNames',Stats);
    T1 = table(Stats,data_table_sim);
    title_t1=strcat(['Simulated.subdaily.' char(stationname(i))]);
    disp (title_t1)
%   figure('Name',['Simulated.subdaily.' char(stationname(i))],'NumberTitle','off');
%   uitable('Data',T{:,:},'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    prettyprint (T1)
    
    data_table_obs=[skew_nz_obs;skew_tot_obs;mean_nz_obs;mean_tot_obs;std_nz_obs;std_tot_obs];
%    T = table(data_table_obs,'RowNames',Stats);
    T2 = table(Stats,data_table_obs);
    title_t2=strcat(['Observation.subdaily.' char(stationname(i))]);
    disp (title_t2)
%   figure('Name',['Observation.subdaily.' char(stationname(i))],'NumberTitle','off');
%   uitable('Data',T{:,:},'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    prettyprint (T2)
end
%% step 4: SAVING THE DISAGGREGATED RESULTS
    for i=1:nstations
      fname=[fullfile('D:\MulRP',['Station' num2str(i,'%4d')],'6h.csv')];  % build a file name
      csvwrite(fname,cell2mat(sim_6h(1,i)))                % open file handle with it% write as unformatted % close this one...
    end
% *************************************************************************
%                   END  DISAGGREGATION OF  PRECIPITATION DATA
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
%% ************************************************************************
%         MULTIVARIATE RETURN PERIOD
% *************************************************************************
% ## Research Question
% 
%   *|*RP _river_*=10 | *RP _river_*=50 | *RP _river_*=100
% --- | --- | ---| ---
% *RP _drain_*=10 | `?` |`?`|`?`
% *RP _drain_*=50| `?` |`?`|`?`
% *RP _drain_*=100| `?` |`?`|`?`
% ## Joint Return Periods 
% 
% ### **OR Case**
% 
% The probability associated with the exceedance of either one of 
% the two variables, that is, the bivariate event EOR(x,y)=[X>x?Y>y] (OR case)
% Due to lack of data on flood water depths, we tried to work with precipitation data
% at upstream and the Monza city. 
% Upper-tail dependence structure of these two data could be generalized to 
% the dependence structure of water depth due to riverine and flash floods.
% We used **upper-tail dependence index** between upstream and at-site precipitation,
%  in [0,1], and set  the estimated value of the **upper-tail dependence index equal to the
%  **Kendall's tau** dependence measure between depth 
%  of water from riverine flood and from flash flood.
% Here, all examined copulas have only one parameter. They are asymetric copulas 
% (with more dependence at upper-tail) including **Flipped Clayton**, **Joe**,
% and **Gumbel**.
% %
% Reference:
% 1)Shiau, J.T. Fitting Drought Duration and Severity with Two-Dimensional Copulas.
% Water Resour Manage 20, 795–815 (2006). https://doi.org/10.1007/s11269-005-9008-9
% 
% 2)Salvadori, G., Durante, F., & De Michele, C. (2011). On the return period and design
% in a multivariate framework. https://www.hydrol-earth-syst-sci.net/15/3293/2011/

%% ************************************************************************
%         MULTIVARIATE RETURN PERIOD CALCULATOR (... STEPS)
% *************************************************************************
clear all
close all
[~,~,nstations,threshold,~]=input_paras();
for i=1:nstations
    fname=[fullfile('D:\MulRP',['Station' num2str(i,'%4d')],'6h.csv')];  % build a file name
    sim_6h{1,i}=csvread(fname);                % open file handle with it% 
end
%% step 1: CALCULATING THE TAIL DEPENDENCE INDEX
sim_6h_upper=(sim_6h{1,2}+sim_6h{1,3})';
sim_6h_monza=(sim_6h{1,1})';
sim_6h_tot=[sim_6h_monza(:),sim_6h_upper(:)];
clearvars sim_6h_monza sim_6h_upper % to save available memory
%tail dep index, Fixed Threshold
corr_prc_fx = fitLambda(sim_6h_tot, 0.95);%95 percent higher values
corr_prc_fx = corr_prc_fx(2,1);
%tail dep index, Optimum Threshold (Frahm et al (2005))
[corr_prc,tr_opt,wb] = optimal_tdc(sim_6h_tot,"upper");
%Kendall's tau. Kendall's_tau=asin(rho)/(0.5*pi)
corr_prc = (asin(corr_prc))/(0.5*pi); 
tdcs=taildepplot(sim_6h_tot,"upper");
chiplot(sim_6h_tot);
%% step 2: CALCULATING THE MULTIVARIATE RETURN PERIOD BASED ON COPULA VALUES
RP = [10,20,50,100];
pr = 1-(1./RP);
[u1,u2] = meshgrid(pr,pr);
% y = copulacdf('Clayton',[u1(:),u2(:)],corr_prc);
%  FOR FLIPPED/ORIGINAL CLAYTON COPULA : ?=2?/(1-?)
y = fclaytoncdf([u1(:),u2(:)], (2*corr_prc)/(1-corr_prc));
RP_copula_fclayton   =  round(100*(1./(1-y)))/100;
%  FOR GUMBEL COPULA : ?=1/(1-?)
y = copulacdf('Gumbel',[u1(:),u2(:)],1/(1-corr_prc));
RP_copula_gumbel   =  round(100*(1./(1-y)))/100;
%  FOR INDEPENDENT CASE
y = u1(:).*u2(:);
RP_ind   =  round(100*(1./(1-y)))/100;
%  FOR JOE COPULA : ?=1/(1-?)
PAR = joetau2par(corr_prc);
y = joecdf([u1(:),u2(:)],PAR);
RP_copula_joe   =  round(100*(1./(1-y)))/100;
% *************************************************************************
%                   END  MULTIVARIATE RETURN PERIOD CALCULATOR
% *************************************************************************
% *************************************************************************
% *************************************************************************
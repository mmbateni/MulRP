%
% generate Tmax and Tmin using spatial correlated random numbers with
% scheme of WeGETS
% 

function [gTmax,gTmin]=multisite_temp_generation(observation,stationname,par,...
    occurrences_gen,nor_ran_tmax,nor_ran_tmin,begin_year,end_year,...
    begin_month,length_month,nstations,years_sim,graph)

% basic treatment for multisite generated precip occurrence
n=years_sim*365;
multi_occ=zeros(n,nstations);
multi_occ=addymd(multi_occ);
for i=1:12
    mm=find(multi_occ(:,2)==i);
    multi_occ(mm,4:end)=occurrences_gen(1,i).occ(:,1:nstations);
end

%% generate Tmax and Tmin using normally distributed random number
t=1:n;
T=365/(2*pi);
gTmax=zeros(n,nstations);
gTmin=zeros(n,nstations);
for ii=1:nstations
    X=multi_occ(:,ii+3)';

    %  time series of min temperature and max temperature 
    %
    %  we first generate the fourier estimates of average and standard deviations
    %  of the time series
    %
    ay=zeros(4,n);
    sy=zeros(4,n);
    for j=1:4
       ay(j,:)=par(ii).aC0(j)+par(ii).aC1(j)*sin(t/T+par(ii).aD1(j))+par(ii).aC2(j)*sin(2*t/T+par(ii).aD2(j));
       sy(j,:)=par(ii).sC0(j)+par(ii).sC1(j)*sin(t/T+par(ii).sD1(j))+par(ii).sC2(j)*sin(2*t/T+par(ii).sD2(j));
    end
    %
    %  procedure is started by assuming that the residuals are equal to 0
    %  random component eps is normally distributed N(0,1)
    %
    res=[0; 0];
    ksi=zeros(2,n);  % fbri7

    eps=zeros(2,n);
    eps(1,:)=nor_ran_tmax(:,ii);
    eps(2,:)=nor_ran_tmin(:,ii);

    for i=1:n   
       res=par(ii).A*res+par(ii).B*eps(:,i);
       ksi(:,i)=res;
    end
    %
    % the means and standard deviations obtained by Fourier time series are conditioned
    % on the wet or dry status of the day determined by using the Markov chain model
    %
    v=ones(2,1);
    XX=kron(v,X);

    cay=XX.*ay(1:2,:)+(1-XX).*ay(3:4,:);
    csy=XX.*sy(1:2,:)+(1-XX).*sy(3:4,:);

    %
    % the Tmax and Tmin are generated conditioned on each other (Jie Chen modified)
    % the smaller standard deviation of Tmax or Tmin is used as a base, and the
    % other parameter is generated conditioned on the chosen parameter. If the
    % standard deviation of Tmax is larger than or equal to the standard
    % deviation of Tmin, daily temperatures are generated by:
    % Tmin=Mean(min)+Std(min)*rand
    % Tmax=Tmin+(Mean(max)-Mean(min))+(Std(max)^2-Std(min)^2)^0.5*rand
    % If the standard deviation of Tmax is less than those of Tmin, daily
    % temperatures are genareted by:
    % Tmax=Mean(max)+Std(max)*rand
    % Tmin=Tmax-(Mean(max)-Mean(min))-(Std(min)^2-Std(max)^2)^0.5*rand

    for i=1:length(ksi)
        if csy(1,i)>=csy(2,i)
            Xp(2,i)=ksi(2,i)*csy(2,i)+cay(2,i);
            Xp(1,i)=ksi(1,i)*(csy(1,i)^2-csy(2,i)^2)^(1/2)+(cay(1,i)-cay(2,i))+Xp(2,i);
        else
            Xp(1,i)=ksi(1,i)*csy(1,i)+cay(1,i);
            Xp(2,i)=ksi(2,i)*(csy(2,i)^2-csy(1,i)^2)^(1/2)-(cay(1,i)-cay(2,i))+Xp(1,i);
        end
    end

    % range control the Tmin, insure that Tmin is always small than Tmax
    for i=1:length(ksi)
        if Xp(1,i)<=Xp(2,i)
            Xp(2,i)=Xp(1,i)-abs(Xp(1,i))*0.2;
        end
    end 
    Tmax=Xp(1,:);
    Tmin=Xp(2,:);

    gTmax(:,ii)=Tmax';
    gTmin(:,ii)=Tmin';

end

%% produce the graph of correlation matrix
%
% at the annual scale
% correlation of observed temperature
if graph==1
    [corr_tmax,corr_tmin]=corr_temp(observation,stationname,nstations,begin_year,end_year);

    % correlation multisite weather generator generated temperature
    gTmax_cor=corr(gTmax);
    gTmin_cor=corr(gTmin);

    % plot the annual correlation matrix
    % extract upper triangular part for correlation matrix
    Tmax_cor=triu(corr_tmax);
    Tmin_cor=triu(corr_tmin);
    gTmax_cor=triu(gTmax_cor);
    gTmin_cor=triu(gTmin_cor);

    % plot the results
    figure
    subplot(1,2,1)
    plot(Tmax_cor,gTmax_cor,'bo')
    hold on
    plot(0:1,0:1,'-')
    xlim([0.8 1])
    ylim([0.8,1])
    xlabel('Observed correlation')
    ylabel('Generated correlation')
    title('Multi-site Tmax')

    subplot(1,2,2)
    plot(Tmin_cor,gTmin_cor,'bo')
    hold on
    plot(0:1,0:1,'-')
    xlim([0.8 1])
    ylim([0.8,1])
    xlabel('Observed correlation')
    ylabel('Generated correlation')
    title('Multi-site Tmin')
end

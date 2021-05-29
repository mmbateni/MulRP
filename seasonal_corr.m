function [corr_amnts]=seasonal_corr(data,stationname,nstations,threshold,...
    begin_season,length_season,season,begin_year,end_year)
%% calculate the spatial correlation matrices for precip amounts for each season  
% extract precip for each station
for i=1:nstations
    pre=data.(char(stationname(i)));
    P=zeros(size(pre));
    P(:,1:31)=pre(:,335:365); % move DEC at the begin for calculating seasonal correlation
    P(:,32:365)=pre(:,1:334);
    eval([char(stationname(i)) '=P;']);   
end
%% preliminary treatment create matrices 'stationname'occ for occurences,
% with NaN as missing data and 'stationname'amounts for amounts with NaN as
% missing data, all using a threshold sepcified above
for i=1:nstations
    eval(['mat=' char(stationname(i)) ';']);    
    amounts=mat;
    occ=zeros(size(mat));
    j=find(mat<=threshold);
    amounts(j)=0;
    eval([char(stationname(i)) 'amounts=amounts;']);    
end
%
% calculate seasonal correlation matrices
%
for i=1:nstations
    for j=1:nstations
        eval(['station_i_amounts=' char(stationname(i)) 'amounts;' ]);
        eval(['station_j_amounts=' char(stationname(j)) 'amounts;' ]);        
        % find intersection of year period between stations       
        yearstart=max([begin_year(i) begin_year(j)]);
        yearend=min([end_year(i) end_year(j)]);
        years_i=[begin_year(i):end_year(i)];
        years_j=[begin_year(j):end_year(j)];
        ifind=find(years_i>=yearstart & years_i<=yearend);
        jfind=find(years_j>=yearstart & years_j<=yearend);
        nyears=length(ifind);
        station_i_amounts=station_i_amounts(ifind,:);
        station_j_amounts=station_j_amounts(jfind,:);
        % calculate seasonal correlations for precip amounts
        for iseason=1:4
           if nyears~=0   % if there is an intersection between station data
               amounts_i=station_i_amounts(:,begin_season(iseason):begin_season(iseason)+length_season(iseason)-1);
               amounts_j=station_j_amounts(:,begin_season(iseason):begin_season(iseason)+length_season(iseason)-1);
               amounts_i=reshape(amounts_i',length_season(iseason)*nyears,1);
               amounts_j=reshape(amounts_j',length_season(iseason)*nyears,1);
               
               ij=find(amounts_i>=0 & amounts_j>=0 );  % find all elements where there is precip data in both
               corr_val=corrcoef(amounts_i(ij),amounts_j(ij));               
               if numel(corr_val)>1%for Matlab
                  corr_val=corr_val(1,2); 
               end
           else
               corr_val=-999;
           end
           corr_amnts(iseason).cor(i,j)=corr_val;
           corr_amnts(iseason).season=char(season(iseason));
        end
        
    end
end
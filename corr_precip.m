% calculate the spatial correlation matrices for precip occurences and
% amounts for each month  
function [corr_occ,corr_amounts]=corr_precip(occs,amounts,stationname,nstations,...
    begin_month,length_month,months,begin_year,end_year)
warning('off','all');
% calculate monthly correlation matrices
corr_occ_mat=zeros(nstations,nstations,12);
corr_amounts_mat=zeros(nstations,nstations,12);
for i=1:nstations
    for j=1:nstations
%       station_i_occ=occs.(char(stationname(i)));
%       station_j_occ=occs.(char(stationname(j)));
%       station_i_amounts= amounts.(char(stationname(i)));
%       station_j_amounts= amounts.(char(stationname(j)));
        eval(['station_i_occ=' 'occs.',char(stationname(i)),';' ]);
        eval(['station_i_amounts=' 'amounts.',char(stationname(i)),';']);
        eval(['station_j_occ=' 'occs.',char(stationname(j)),';']);
        eval(['station_j_amounts=' 'amounts.',char(stationname(j)),';' ]);        
%       find intersection of year period between stations       
        yearstart=max([begin_year(i) begin_year(j)]);
        yearend=min([end_year(i) end_year(j)]);
        years_i=[begin_year(i):end_year(i)];
        years_j=[begin_year(j):end_year(j)];
        ifind=find(years_i>=yearstart & years_i<=yearend);
        jfind=find(years_j>=yearstart & years_j<=yearend);
        nyears=length(ifind);
        station_i_occ=station_i_occ(ifind,:);
        station_i_amounts=station_i_amounts(ifind,:);
        station_j_occ=station_j_occ(jfind,:);
        station_j_amounts=station_j_amounts(jfind,:);
        % calculate monthly correlations for both occurences and amounts
        for imonth=1:12
           % start with occurence
           if nyears~=0   % if there is an intersection between station data
               occ_i=station_i_occ(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
               occ_j=station_j_occ(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
               occ_i=reshape(occ_i',length_month(imonth)*nyears,1);
               occ_j=reshape(occ_j',length_month(imonth)*nyears,1);
               ij=find(occ_i>=0 & occ_j>=0 );  % find all elements where there is no missing values in both
               corr=corrcoef(occ_i(ij),occ_j(ij));
               corr=corr(1,2);
           else
               corr=-999;
           end
           if isnan(corr)
               corr=prod(corr_occ(imonth).cor(i-1,1:j));
           end
           eval(['corr_occ_' char(months(imonth)) '(i,j)=corr;']);
           corr_occ(imonth).cor(i,j)=corr;
           corr_occ(imonth).month=char(months(imonth));
           corr_occ_mat(i,j,imonth)=corr;
           % do amounts
           if nyears~=0   % if there is an intersection between station data
               amounts_i=station_i_amounts(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
               amounts_j=station_j_amounts(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
               amounts_i=reshape(amounts_i',length_month(imonth)*nyears,1);
               amounts_j=reshape(amounts_j',length_month(imonth)*nyears,1);
               ij=find(amounts_i>=0 & amounts_j>=0 );  % find all elements where there is precip data in both
               corr=corrcoef(amounts_i(ij),amounts_j(ij));               
               corr=corr(1,2);
           else
               corr=-999;
           end
         corr_amounts(imonth).cor(i,j)=corr;
         corr_amounts(imonth).month=char(months(imonth));
         corr_amounts_mat(i,j,imonth)=corr;
        end    
    end  
end
    figure('Name','Historical Data Occurence','NumberTitle','off');
for imonth=1:12
    subplot(4,3,imonth);
    plot_corrmat(corr_occ_mat(:,:,imonth),['Corr.Occ.' char(months(imonth))])
end
    figure('Name','Historical Data Amounts','NumberTitle','off');
for imonth=1:12
    subplot(4,3,imonth);
    plot_corrmat(corr_amounts_mat(:,:,imonth),['corr.Amounts(nzero).' char(months(imonth))])
end
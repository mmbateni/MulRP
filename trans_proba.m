%% calculate monthly p00 and p10 parameters for precip occurrence, and the
% mean and stdev of the distribution for precip amounts
function [trans_prob,amounts,occs]=trans_proba(data,stationname,nstations,...
    threshold,begin_month,length_month,months)
% preliminary treatment create matrices 'shortname'occ for occurences, with
% NaN as missing data using the threshold specified above
% fitting 3rd parameter of the distribution to all data of stations for
% each month
f = fieldnames(data)';
c = cell(length(f),1);
s = cell2struct(c,f);
amounts = s;
occs = s;
for i=1:nstations
    eval(['mat=' 'data.',char(stationname(i)) ';' ]);    
    amount=mat;
    occ=zeros(size(amount));    
    j=amount>threshold;
    occ(j)=1;
    j=amount<=threshold;
    amount(j)=0;    
    amounts.(char(stationname(i)))=amount;
    occs.(char(stationname(i)))=occ;
%   amounts_tot(1:size(amounts,1),:,i)=amount;
end
%
%% calculate transition matrices and the distribution
% 
Stats = {'p00';'p10'};
for i=1:nstations
    amount=amounts.(char(stationname(i)));   
    occ=occs.(char(stationname(i)));
    p_trans=zeros(2,12);
    figure('Name',['Historic. dist. of wet days ',char(stationname(i))],'NumberTitle','off');
    for imonth=1:12
        % occurence
        occmonth=occ(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
        amountsmonth=amount(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
        [p00,p10,~] = transition(occmonth);
        p_trans(:,imonth)=[p00;p10];
        eval(['trans_prob_' char(months(imonth)) '(1,i)=p00;']);
        eval(['trans_prob_' char(months(imonth)) '(2,i)=p10;']);
        trans_prob(imonth).prob(1,i)=p00;
        trans_prob(imonth).prob(2,i)=p10;
        trans_prob(imonth).month=char(months(imonth));
        qq=amountsmonth>0;
        subplot(3,4,imonth)
        hist(amountsmonth(qq));     
    end
%            T = table(p_trans(:,1),p_trans(:,2),p_trans(:,3),p_trans(:,4),...
%            p_trans(:,5),p_trans(:,6),p_trans(:,7),p_trans(:,8),p_trans(:,9),...
%            p_trans(:,10),p_trans(:,11),p_trans(:,12),'VariableNames',cellstr(char(months)),...
%            'RowNames',Stats);
             TP_months=[p_trans(:,1),p_trans(:,2),p_trans(:,3),p_trans(:,4),...
             p_trans(:,5),p_trans(:,6),p_trans(:,7),p_trans(:,8),p_trans(:,9),...
             p_trans(:,10),p_trans(:,11),p_trans(:,12)];
             title=strcat(['Historic.' char(stationname(i))]);
             disp (title)
             T = table(Stats,TP_months);
             prettyprint (T)
%            figure('Name',['Historic.' char(stationname(i))],'NumberTitle','off');
%            uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
%            'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
end







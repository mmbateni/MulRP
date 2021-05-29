%Challenge:
% problem arises in that the normally distributed random numbers generated are normally distributed
% whereas p00 and p11 require a uniform random number. The simplest solution is
% to transform the normally distributed random numbers into uniform ones using
% the cumulative distribution function. Alternatively,
% p00 and p11 could be transformed into normally distributed numbers. 
% A second problem arises in that the generated occurrences are less correlated
% than the observed ones. In other words, in order to generate occurrences 
% with the same correlation as the observed ones, it is necessary to use random numbers
% that are more correlated than the occurrences.
% What this code does: generate spatially correlated precip occurrence with two steps: 
% (1) automatic determination of correlation matrix of random number using
% algorithm.Diagonalization and repalcement of negative eigenvalues
% if matrix is non-positive definite 
% (2) generate precipitation occurrence using first-order Markov chain with
% Cholesky factorization

function [corr_occ_rand,corr_occ_gen,occurrences_gen,normal_rand]=multisite_occ_generation(...
      corr_occ,trans_prob,stationname,nstations,length_month,years_sim,months,Tolerance,...
      maxiter,seed)
warning('off','all');
n=length_month*years_sim;       % number of days in the iterative process
prob_occ_gen=zeros(2,12,nstations);
randn ("seed",seed)% for setting Octave random seed
% obtain correlation matrix of occurences
for q=1:12    % do for each month
    C=corr_occ(1,q).cor;
	  transitions=trans_prob(1,q).prob;    
  	% since random numbers generated have a normal distribution, each p00 and
  	% p10 has to be recalculated according to a normal number	
  	transitions_normal=zeros(size(transitions));
    for i=1:nstations
        transitions_normal(1,i)=norminv(transitions(1,i));%-sqrt(2)*erfcinv(2*transitions(1,i));
        transitions_normal(2,i)=norminv(transitions(2,i));%-sqrt(2)*erfcinv(2*transitions(2,i));
    end
    % produce independently normally ditributed random number
	  random=zeros(nstations,n(q));
	  for i=1:nstations
        random(i,:)=randn(1,n(q));
	  end
	% automatic determination of correlation matrix of random number using
  % algorithm of eq.7 of Brissette, Khalili, and Leconte, 2007.
  % tttt=cputime;
	  [M,K,norm_random]= spatial_iterate_occ(C,nstations,random,...
        transitions_normal,n(q),years_sim,Tolerance,maxiter);
  % time_in_second=cputime-tttt;
    occurrences=zeros(n(q),nstations);
    for i=1:nstations
        occurrences_i=occ_generate(transitions_normal(:,i),norm_random(:,i),length_month(q));	%'correlation of random numbers needed'
        occ1=reshape(occurrences_i,length_month(q),[]);
        occ2=occ1';
        [p00,p10,~] = transition(occ2);
        prob_occ_gen(:,q,i)=[p00;p10];
        occurrences(:,i)=occurrences_i;
    end
    eval(['corr_occ_rand_' char(months(q)) '=M;' ]);
	% resulting correlation of occurence
    eval(['corr_occ_gen_' char(months(q)) '=K;' ]);
    eval(['occurrences_' char(months(q)) '=occurrences;' ]);
    eval(['normal_rand_' char(months(q)) '=norm_random;' ]); 
	  corr_occ_rand(q).cor=M;
    corr_occ_rand(q).month=char(months(q));    
    corr_occ_gen(q).cor=K;
    corr_occ_gen(q).month=char(months(q));   
    occurrences_gen(q).occ=occurrences;
    occurrences_gen(q).month=char(months(q));
    normal_rand(q).rand=norm_random;
end
%
% produce graphics
%
LName = {'p00';'p10'};
for k=1:nstations
%            T = table(prob_occ_gen(:,1,k),prob_occ_gen(:,2,k),prob_occ_gen(:,3,k),prob_occ_gen(:,4,k),...
%                prob_occ_gen(:,5,k),prob_occ_gen(:,6,k),prob_occ_gen(:,7,k),prob_occ_gen(:,8,k),prob_occ_gen(:,9,k),...
%                prob_occ_gen(:,10,k),prob_occ_gen(:,11,k),prob_occ_gen(:,12,k),'VariableNames',cellstr(char(months)),...
%                'RowNames',LName);
             title_t=strcat(['Simulated.' char(stationname(k))]);
             disp (title_t)
             TP_months=[prob_occ_gen(:,1,k),prob_occ_gen(:,2,k),prob_occ_gen(:,3,k),prob_occ_gen(:,4,k),...
             prob_occ_gen(:,5,k),prob_occ_gen(:,6,k),prob_occ_gen(:,7,k),prob_occ_gen(:,8,k),prob_occ_gen(:,9,k),...
             prob_occ_gen(:,10,k),prob_occ_gen(:,11,k),prob_occ_gen(:,12,k)];
             T = table(LName,TP_months);
             prettyprint (T)
%            figure('Name',['Simulated.' char(stationname(k))],'NumberTitle','off');
%            uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
%            'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
end
% correlation of generated precip occurrence
figure('Name','occurrence correlation of historical and generated precip','NumberTitle','off');
for ii=1:12
        C=corr_occ(1,ii).cor;
        CC=reshape(triu(C,1),1,[]);
        i=find(CC~=0);
        CC=CC(i);
        K=corr_occ_gen(1,ii).cor;
        KK=reshape(triu(K,1),1,[]);
        i=find(KK~=0);
        KK=KK(i);
        subplot(3,4,ii)
        plot(CC,KK,'o',[0 1],[0 1]);
        title_p=strcat(['all stations, month ' char(months(ii))]);
        title(title_p);
        xlabel('observed correlation');
        ylabel('generated correlation');    
        axis ([0 1 0 1], "square");
        set(gca, "fontsize", 12)
        set(gcf,'Color',[1 1 1])      
%       set(gca,'FontSize',12)
end
%     % some other graphics commands 
%         xlabel('observed correlation');
%         ylabel('random numbers correlation');
%         axis([0 1 0 1])
%         axis square
%         set(gca,'FontSize',12)
%         set(gcf,'Color',[1 1 1])





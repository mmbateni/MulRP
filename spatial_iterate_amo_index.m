function [mat,K,prec] = spatial_iterate_amo_index(C,phat_multiexp,dist,nstations,...
    random,n,threshold,occurrences,Char_season,Tolerance,maxiter)
% this function iterates from a starting correlation matrix, and modifies it
% so that it becomes the correlation matrix of random numbers giving
% an amount correlation matrix equal to the original
%
for i=1:nstations
    stationname{i,1}=['S',num2str(i)]; 
end
val=5;
ii=0;
kiter=0.1;   % iteration parameter in calculation of new estimate of matrix 'mat'
mat=C;  % start with observed correlation matrix as best estimate of random numbers correlation matrix
count=0;
corr=zeros(nstations,nstations);

while val>Tolerance & ii<maxiter;
      ii=ii+1;         % 'live' counter fo follow progress on screen
      q=find(mat>1);  % avoid solutions with correlation coefficient greater than 1
      mat(q)=0.999;
      % insure that matrix is positive semi-definite (all eigen values >= 0).  if
      % not, cholesky deconposition will not work
      [V,D] = eig(mat);
      i=find(D<0);      % find eigenvalues smaller than 0
      if length(i)>=1 
            count=count+1;
            D(i)=0.000001;    % they should be only a little bit below zero
            Q=(V*D*inv(V));   % reconstruct matrix, which becomes a covariance matrix diagonal slighly larger than 1
                              % we keep on iterating with the covariance matrix and
                              % only noramlize at the end, allows for better accuracy
            %mat=Q./sqrt(diag(Q)*diag(Q)');   % transform covariance matrix in correlation matrix that is positive definite
            mat=real(Q);
      end
      U=chol(mat);
      L=U';% lower diagonal matrix  L*U=C
      values=L*random;
      values=values';    
      corr_random=values;
      mean_random=mean(corr_random);
      std_random=std(corr_random);
      unit_vector=ones(n,1);
      mean_matrix=unit_vector*mean_random;
      std_matrix=unit_vector*std_random;
      norm_random=(corr_random-mean_matrix)./(std_matrix);   % this matrix contains the N(0,1) serially independant but spatially correlated random numbers 
      %
      % at this point, for each station we will generate a series of length
      % 'n' of precipitation using a matrix of occurrrences already determined
      prec=zeros(n,nstations);
      randU01=0.5*erfc(-norm_random/sqrt(2)); % takes normal number and transform back into 0-1 number using CDF
      for i=1:nstations
            ParameterValues=eval(['phat_multiexp.',char(stationname(i)),'.phat',Char_season]);
            DistName=eval(['dist.',char(stationname(i)),'.Dist',Char_season]);
    %       p = -expm1(log1p(-randU01(:,i))/phat_mexp(1)); %if lower.tail=T p=randU01(:,i).^(1/phat_mexp(1));
            prec_t = distinv(randU01,DistName,ParameterValues,threshold);%threshold replaced by 0
    %       p = randU01(:,i).^(1/phat_mexp(1));
    %       prec_t = gpinv(p,phat_mexp(3),phat_mexp(2),0);%threshold replaced by 0
    %       prec_t = gpinv(p,phat_mexp(1),phat_mexp(2),phat_mexp(3));%
    %       res = log(phat_mexp(1)) + gpinv(randU01(:,i),phat_mexp(3),phat_mexp(2),threshold) + 
    %       (phat_mexp(1) - 1) .* gpcdf(randU01(:,i), phat_mexp(3),phat_mexp(2),threshold)
    %       exp_quantiles = expinv(p,1);%MU=1
    %       prec_t=phat_mexp(2) *((expm1(exp_quantiles * phat_mexp(3))./(exp_quantiles * phat_mexp(3))).* exp_quantiles)+ threshold;
            for j=1:n
                if  occurrences(j,i)==1 
                    prec(j,i)=prec_t(j);
                end
            end
      end
    % calculate correlation of occurences using only data when there is
    % precip at each stations pair simultaneously
      for i=1:nstations
            for j=1:nstations
                  ij=find(prec(:,i)>=0 & prec(:,j)>=0 );  % find all elements where there is data in both
                  T=corrcoef(prec(ij,i),prec(ij,j));
                  if numel(T)>1%for Matlab
                      T=T(1,2); 
                  end
                  corr(i,j)=T;
            end
      end  
      K=corr;
      val=max(max(abs(K-C))); 
      if ii~=maxiter & val>Tolerance % change value of random number matrix, only if computations are to continue
            for i=1:nstations
                for j=1:nstations
                   if i==j
                       mat(i,j)=1;
                   else
                       mat(i,j)=mat(i,j)+kiter*(C(i,j)-K(i,j));
                   end
                end
            end
      end
        val;
end
mat=mat./sqrt(diag(mat)*diag(mat)'); % have to really think whether or not this line is needed.

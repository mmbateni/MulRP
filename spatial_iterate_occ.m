function [mat,K,norm_random] = spatial_iterate_occ(C,nstations,random,...
    transitions_normal,nmonth,years_sim,Tolerance,maxiter)
%
% this function iterates from a starting correlation matrix, and modifies it
% so that it becomes the correlation matrix of random numbers giving
% an occurence correlation matrix equal to the original
%
%
val=5;
ii=0;
kiter=0.1;   % iteration parameter in calculation of new estimate of matrix 'mat'
mat=C;  % start with observed correlation matrix as best estimate of random numbers correlation matrix
while val>Tolerance && ii<maxiter
      ii=ii+1;         % 'live' counter fo follow progress on screen
	    q=find(mat>1);  % avoid solutions with correlation coefficient greater than 1
	    mat(q)=0.999;
	    % insure that matrix is positive semi-definite (all eigen values >= 0).  if
	    % not, cholesky deconposition will not work
	    [V,D] = eig(mat);
      D_vec = diag(D);
	    i=find(D_vec<=0);      % find eigenvalues smaller than 0
      if length(i)>=1 
          D_vec(i)=0.00001;    % they should be only a little bit above zero
          D= diag(D_vec);
          Q=(V*D*inv(V));   % reconstruct matrix, which becomes a covariance matrix diagonal slighly larger than 1
          % we keep on iterating with the covariance matrix and
          % only noramlize at the end, allows for better accuracy
          %mat=Q./sqrt(diag(Q)*diag(Q)');   % transform covariance matrix in correlation matrix that is positive definite
          mat=Q;
      end
      U=chol(mat);
      L=U';   % lower diagonal matrix  L*U=C
      values=L*random;
      values=values';     
      corr_random=values; 
      mean_random=mean(corr_random);
      std_random=std(corr_random);
      unit_vector=ones(nmonth,1);
      mean_matrix=unit_vector*mean_random;
      std_matrix=unit_vector*std_random;
      norm_random=(corr_random-mean_matrix)./(std_matrix);   % this matrix contains
      % N(0,1) serially independant but spatially correlated random numbers 
      % at this point, for each station,a series of length 'n' of precip
      % occurence will be generated
      occ=zeros((nmonth/years_sim),years_sim,nstations);
      for i=1:nstations
          for k=1:years_sim
          for j=2:(nmonth/years_sim)	
                j_k=j+((k-1)*(nmonth/years_sim));
                if occ(j-1,k,i)==0 && (norm_random(j_k,i)<= transitions_normal(1,i))
                    occ(j,k,i)=0;
                elseif occ(j-1,k,i)==1 && (norm_random(j_k,i)<= transitions_normal(2,i))
                    occ(j,k,i)=0;
                else
                    occ(j,k,i)=1;
                end
        end
        end
      end  
      occ=reshape(occ,nmonth,nstations);
	    K=corrcoef(occ);   
      val=max(max(abs(K-C))); 
      if ii~=maxiter & val>Tolerance    % change value of random number matrix, only if computations are to continue
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
end
mat=mat./sqrt(diag(mat)*diag(mat)');   % have to really think whether or not this line is needed.
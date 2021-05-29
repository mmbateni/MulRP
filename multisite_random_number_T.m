% produce normally distributed random number for generating temperature,
% the induced the correlation is expected to be equal to the observed
% correlation

function [nor_ran_tmax,nor_ran_tmin]=multisite_random_number_T(corr_tmax,...
    corr_tmin,nstations,years_sim,days)

 %
% insure that matrix is positive semi-definite (all eigen values >= 0).  if
% not, cholesky deconposition will not work

n=years_sim*days; % the days will be generated

% for Tmax
[V,D] = eig(corr_tmax);
i=find(D<0);      % find eigenvalues smaller than 0
if length(i)>=1 

    D(i)=0.000001;    % they should be only a little bit below zero
    Q=(V*D*inv(V));   % reconstruct matrix, which becomes a covariance matrix diagonal slighly larger than 1
                      % we keep on iterating with the covariance matrix and
                      % only noramlize at the end, allows for better
                      % accuracy
%    Q=Q./sqrt(diag(Q)*diag(Q)');   % transform covariance matrix in correlation matrix that is positive definite
    corr_tmax=real(Q);
end

% Cholesky factorization
LL=chol(corr_tmax);
LLL=ctranspose(LL);
L1=LLL*LL;
L = chol(L1,'lower');
        
% produce correlated nomally distributed random numbers for generating Tmax
nor_ran_tmax=L*(randn(nstations,n));
nor_ran_tmax=nor_ran_tmax';

% for Tmin
[V,D] = eig(corr_tmin);
i=find(D<0);      % find eigenvalues smaller than 0
if length(i)>=1 

    D(i)=0.000001;         % they should be only a little bit below zero
    Q=(V*D*inv(V));   % reconstruct matrix, which becomes a covariance matrix diagonal slighly larger than 1
                      % we keep on iterating with the covariance matrix and
                      % only noramlize at the end, allows for better
                      % accuracy
    Q=Q./sqrt(diag(Q)*diag(Q)');   % transform covariance matrix in correlation matrix that is positive definite
    corr_tmin=real(Q);
end

% Cholesky factorization
LL=chol(corr_tmin);
LLL=ctranspose(LL);
L1=LLL*LL;
L = chol(L1,'lower');
        
% produce correlated nomally distributed random numbers for generating Tmin
nor_ran_tmin=L*(randn(nstations,n));
nor_ran_tmin=nor_ran_tmin';


function cdf = empirical2cdf(X,i_n,j_n)
% EMPIRICALCDF Empirical cdf function for bivariate case
% Function Returned values are in the interval (0,1).
%
%   References:
%   Berg, D. Bakken, H. (2006) Copula Goodness-of-fit Tests: A Comparative Study
%
n = size(X, 1);
i = floor(n * i_n);
j = floor(n * j_n);
sortX1=sort(X(:,1));
sortX2=sort(X(:,2));
ith_order_u = sortX1(i);
jth_order_v = sortX2(j);
cdf= (1/n) * sum((X(:,1) <= ith_order_u ) .* (X(:,2) <= jth_order_v) );
end
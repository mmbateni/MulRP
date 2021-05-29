 function c = fclaytoncdf(u,theta)
% fclaytoncdf Copula cumulative distribution function.
%C = fclaytoncdf(U,THETA) returns the copula cumulative
%distribution function of the bivariate copula family FAMILY with
%parameter THETA at the values in U. U has the size [N, 2], where N is
%the number of samples. FAMILY can be one of
if theta == 0
   c = prod(u,2);
else
   c = u(:,1) + u(:,2) - 1 + max((1-u(:,1)).^(-theta) + (1-u(:,2)).^(-theta) - 1,0).^(-1./theta);
end
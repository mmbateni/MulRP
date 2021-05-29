
function Lam=fitLambda(u,p)
%Computing non-parametric estimators of the (matrix of) tail-dependence coefficients
%u 2col-matrix of (pseudo-)observations in [0,1]^d for estimating the (matrix of) tail-dependence coefficients.
%the method with which the tail-dependence coefficients are computed:
% method = "Schmid.Schmidt":
% nonparametric estimator of Schmid and Schmidt (2007) 
% (see also Jaworksi et al. (2009, p. 231)) computed for all pairs.
u_  =  1 - u;
d=size(u_,2);
Lam = diag(ones(d,1));
M = reshape(max(0, p - u_), [],d);
%Lam(2, 1) = mean(M(:,1).*M(:,2));
for i=2:d 
    for j =1:(i - 1) 
        Lam(i, j)=mean(M(:, i).*M(:, j));
    end
end
int_over_Pi = (p^d/2)^d;
int_over_M = p^(d+1)/(d+1);
Lam = (Lam - int_over_Pi)./(int_over_M - int_over_Pi);
Lam(Lam < 0) = 0;
Lam(Lam > 1) = 1;
t_Lam=Lam';
Lam(1,2)=t_Lam(1,2);


function PAR = joetau2par(tau)
%% Convert Kendall's tau to parameter value of Joe Copula
if (abs(tau) > 0.99999) 
    PAR=Inf;
elseif (tau < 0)
    PAR=1.000001;
else
     myfun = @(x) tau-tauF(x);
     fun = @(x) myfun(x);
     options = optimset('TolFun',1e-8);
     PAR = fzero(fun,[1,10],options);
%    PAR = uniroot(function(x) tau - tauF(x), lower = 1, upper = 5e+05,
%    tol = .Machine$double.eps^0.5)$root
end
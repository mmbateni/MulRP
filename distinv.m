
function [pt] = distinv(randU01,DistName,ParameterValues,trhold)
%
% Inverse cumulative distribution function of distributions(the quantile (the inverse of the CDF))
% For each element of randU01, computes the cdf at randU01 of the distribution with parameters "ParameterValues"
% for case of extended generalized Pareto:
% x = gpinv(p,k,sigma,theta) returns the inverse cdf for a extended generalized Pareto (GP) distribution 
% with tail index (shape) parameter k, scale parameter sigma,
% and threshold (location) parameter theta of zero, evaluated at the values in p.
%
names_dist={'exponential','gamma','extreme value' ,'lognormal','normal','weibull'};%'generalized extreme value' 'generalized pareto' 
flag_cencored=1;
DistName=char(DistName);
if   strcmp('extended generalized pareto',DistName)
     phat_mexp=ParameterValues;
     p = randU01.^(1/phat_mexp(1));
     pt = gpinv(p,phat_mexp(3),phat_mexp(2),trhold);%threshold replaced by 0
elseif strcmp(names_dist{1},DistName)
     pt =expinv(randU01,ParameterValues);
    if flag_cencored==1  
        pt((pt<trhold))=trhold;
    end
else
    a=ParameterValues(1);
    b=ParameterValues(2);
    switch (DistName)
      case names_dist{2}
            pt =gaminv(randU01,a,b);
      case names_dist{3}
            pt =gevinv(randU01,0,b,a);
      case names_dist{4}
            pt =logninv(randU01,a,b);
      case names_dist{5}
            pt =norminv(randU01,a,b);
      case  names_dist{6}
            pt =wblinv(randU01,a,b);
    endswitch
    if flag_cencored==1
        pt((pt<trhold))=trhold;
    end
end
if any(pt<0)
   pt((pt<0))=0;
end
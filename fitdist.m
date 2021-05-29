function [Params,NLogL] = fitdist(x,distname)
%FITDIST Fit probability distribution to data.
%   PD = FITDIST(X,DISTNAME) fits the probability distribution DISTNAME to
%   the data in the column vector X, and returns the fitted distribution 
%   (parameters & -log liklihood)
%   DISTNAME can be any of the following parametric distribution
%   names:
%
%         'exponential'                      Exponential
%         'extreme value',                   Extreme value
%         'gamma'                            Gamma
%         'lognormal'                        Lognormal
%         'normal'                           Normal
%         'weibull',                         Weibull
%
names_dist={'exponential','gamma','extreme value' ,'lognormal','normal','weibull'};%'generalized extreme value' 'generalized pareto' 
k=numel(x);
    switch (distname)
      case names_dist{1}
            ParameterValues= expfit(x);
            [nlogL] = explike(ParameterValues,x);
            % n=1;
      case names_dist{2}
            ParameterValues=  gamfit(x);
            [nlogL] = gamlike(ParameterValues,x);
      case names_dist{3}
            ParameterValues= evfit(x);
            [nlogL] = evlike(ParameterValues,x);
      case names_dist{4}
            ParameterValues= lognfit(x);
            [nlogL] = lognlike(ParameterValues,x);
      case  names_dist{5}
            [ParameterValues(1),ParameterValues(2)]= normfit(x);
            [nlogL] = normlike(ParameterValues,x);
      case  names_dist{6}
            ParameterValues= wblfit(x);
            [nlogL] = wbllike(ParameterValues,x);
      endswitch
Params=ParameterValues;
NLogL=nlogL;
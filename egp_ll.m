%
% This function provides the log liklihood for the third pdf family presented in 
% Papastathopoulos and Tawn (2013).
% parameters are and Kappa-transformation  > 0, sigma-scale > 0  and 
% xi-shape (any real) ,respectively.
% construct the distribution for 
%

function [ll_egp3]=egp_ll(observation,threshold,parinit)
% if ~exist('parinit','var')
% parinit=[0.5,0.5,0];
% end
%   observation = observation(observation > threshold);
    kappa = parinit(1);
    sigma = parinit(2);
    xi = parinit(3);
%     args = max(0, (1 + xi * (observation - threshold)./sigma));
%     argsp = args.^(-1/xi);
%     if (abs(xi) > 1e-08) 
%         ll_egp3 = (length(observation) * (log(kappa) - log(sigma))) + ((kappa - 1) * sum(log(1 - argsp + eps))) - ((1/xi + 1) * sum(log(args)));
%     else 
%         ll_egp3 = (length(observation) * (-log(sigma) + log(kappa))) + ((kappa - 1) * sum(log(1 - exp(-(observation - threshold)/sigma)))) -sum(observation - threshold)./sigma;
%     end
x_cens = observation(observation < threshold);
x_not_cens = observation(observation >= threshold );
  
    if (kappa <= 0 || sigma <= 0 || xi <= 10^(-6) || xi > 0.99) 
        ll_egp3=Inf;
    else
    censor=gpcdf(threshold,xi,sigma,0)^kappa;
    if (length(x_cens) > 0) %#ok<ISMT>
        contrib_cens=length(x_cens) * log(censor);
    else
        contrib_cens=0;
    end
    %    contrib_not_cens=sum((log(kappa) + ((kappa - 1) .* log(gpcdf(x_not_cens,xi,sigma,0))))+ log(gppdf(x_not_cens,xi,sigma,0)),'omitnan');
        contrib_not_cens=sumskipnan((log(kappa) + ((kappa - 1) .* log(gpcdf(x_not_cens,xi,sigma,0))))+ log(gppdf(x_not_cens,xi,sigma,0)));
    %    contrib_not_cens=sum(log((pgpd(x_not_cens+rounded, scale = sigma, shape = xi)^kappa) - (pgpd(q, scale = sigma, shape = xi)^kappa) ));
        ll_egp3=-(contrib_cens + contrib_not_cens);  
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
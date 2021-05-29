function [p_t] = egpinv(randU01,phat_mexp,trehold)
%
% extended Generalized Pareto inverse cumulative distribution function
% x = gpinv(p,k,sigma,theta) returns the inverse cdf for a extended generalized Pareto (GP) distribution 
% with tail index (shape) parameter k phat_mexp(3),scale parameter sigma phat_mexp(2),
% and threshold (location) parameter theta of zero, evaluated at the values in p.
% The size of x is the common size of the input arguments. A scalar input functions as a
% constant matrix of the same size as the other inputs.
        p = randU01.^(1/phat_mexp(1));
        p_t = gpinv(p,phat_mexp(3),phat_mexp(2),trehold);%threshold replaced by 0
%
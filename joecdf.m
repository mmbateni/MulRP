function P = joecdf(u_,PAR)
%% Joe Copula Family CDF Values
% Based on Wikipedia: https://en.wikipedia.org/wiki/Copula_(probability_theory)
% Condition 1 <= PAR should not be violated
u = u_(:,1);
v = u_(:,2);
P = 1 - ( (1-u).^PAR + (1-v).^PAR - (1-u).^PAR .* (1-v).^PAR ) .^ (1/PAR);
 
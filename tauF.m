    function tau = tauF(par)
%% for converting tau to parameter in Joe copula
% Condition 1 <= PAR_Joe should not be violated
      param1 = 2/par + 1;
      tem = psi(2) - psi(param1);
      tau = 1 + tem * 2/(2 - par);
      tau(par == 2) = 1 - psi(2);

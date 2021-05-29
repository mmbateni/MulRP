
function res=taildep(u,p_u)
%Computing non-parametric estimators of the (matrix of) tail-dependence coefficients
%n x 2-matrix of (pseudo-)observations in [0,1]^d for estimating the (matrix of) tail-dependence coefficients.
%the method with which the tail-dependence coefficients are computed:

% Beirlant, J., Goegebeur, Y., Teugels, J. and Segers, J. (2004) Statistics of Extremes: Theory and Applications. Chichester, West Sussex, England, UK: Wiley, ISBN 9780471976479, 522pp.
% Coles, S. (2001) An introduction to statistical modeling of extreme values, London: Springer-Verlag.
% Coles, S., Heffernan, J. E., and Tawn, J. A. (1999) Dependence measures for extreme value analyses. Extremes, 2, 339–365.
% de Haan, L. and Ferreira, A. (2006) Extreme Value Theory: An Introduction. New York, NY, USA: Springer, 288pp.
% Reiss, R.-D. and Thomas, M. (2007) Statistical Analysis of Extreme Values: with applications to insurance, finance, hydrology and other fields. Birkh\"auser, 530pp., 3rd edition.
% Sibuya, M. (1960) Bivariate extreme statistics. Ann. Inst. Math. Statist., 11, 195–210.
%     type = c("chi", "chibar"), 
x=u(:,1); 
y=u(:,2); 
good =(~isnan(x) & ~isnan(y));
x=x(good);
y=y(good);
n=length(x);
xun=sort(x);
yun=sort(y);
xun=xun(floor(n * p_u));
yun=yun(floor(n * p_u));
id = (x > xun) & (y > yun);
chi = sum(id)/(n * (1 - p_u));
chibar = 2 * (log(1 - p_u)/log(mean(id))) - 1;
res = [chi, chibar];

  

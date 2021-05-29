function tdc_value = TDC(X,i_n,region)
%% Lower/Upper Tail Dependence Coefficient for a given threshold i/n
if region=="lower"
tdc_value=empirical2cdf(X,i_n,i_n)/i_n;
end
if region=="upper"
tdc_value=(empirical2cdf(X,i_n,i_n)+(1-2*i_n))/(1-i_n);
end

function  term=taildepplot(u,region)
%chiplot : Chi-squared plot for independence
%chi. Returned values are in the interval (0,1).
%   References: Aghakouchak, A., Ciach, G., & Habib, E. (2010). Estimation of tail dependence coefficient
%   in rainfall accumulation fields. Advances in Water Resources, 33(9), 1142–1149.
%% 
u = u(sum(u,2)>0,:);
n = size(u, 1);
u_s = u(randi(n,n,1),:);
n1 = floor(0.60*n);
tdc_values(1:n1)=NaN;
tdc_values(n)=NaN;
i_n=(1:n)/n;
if region=="lower"
    for i =(n1+1):n-1
        tdc_values(i)=empirical2cdf(u_s,i_n(i),i_n(i))/i_n(i);
    end
end
if region=="upper"
        for i =(n1+1):n-1
            tdc_values(i)=(empirical2cdf(u_s,i_n(i),i_n(i))+(1-2*i_n(i)))/(1-i_n(i));
        end
end
term=[i_n((n1+1):n-1);tdc_values((n1+1):n-1)];
plot((1:n),tdc_values);
title('Tail dependence index Versus Threshold')
end
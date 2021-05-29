function [ind3] = pattern3(X)
% 3-day pattern recognition
% Reference:Park, H., & Chung, G. (2020). A Nonparametric Stochastic Approach for Disaggregation of Daily to Hourly Rainfall Using 3-Day Rainfall Patterns. Water, 12(8), 2306.
%  ...doi:10.3390/w12082306
ind_daily=zeros(length(X),1);
diff_daily=diff(X);
for j=1:(length(X)-2)
    if diff_daily(j)>0 && diff_daily(j+1)>0
        ind_daily(j+1)=2;
    end
    if diff_daily(j)>0 && diff_daily(j+1)<0
        if X(j+2)>X(j)
           ind_daily(j+1)=6;
        end
        if X(j+2)<X(j)
           ind_daily(j+1)=5;
        end
        if X(j+2)==X(j)
            ind_daily(j+1)=7;
        end
    end
    if diff_daily(j)<0 && diff_daily(j+1)<0
       ind_daily(j+1)=1;
    end
    if diff_daily(j)<0 && diff_daily(j+1)>0
       if X(j+2)>X(j)
          ind_daily(j+1)=4;
       end
       if X(j+2)<X(j)
          ind_daily(j+1)=3;
       end
    end
end
if  X(1)>0
    ind_daily(1)=7 ;
end
if  X(end)>0
    ind_daily(end)=7 ;
end
ind3=ind_daily;
end
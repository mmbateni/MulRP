function [tdc_value,kn,b] = optimal_tdc(X,region)
% Returns optimal Empirical Tail Dependence coefficient (TDC)
% Based on the heuristic plateau-finding algorithm from Frahm et al (2005) 
% "Estimating the tail-dependence coefficient: properties and pitfalls"
fprintf('initializing the optimal_tdc function\n');
ind_nz=(X(:,1)>0 &  X(:,2)>0);
data=X(ind_nz,:);
% data = X ; %#creates a copy 
q85=quantile(data,0.60,1);
ind_ext=(data(:,1)>q85(1) &  data(:,2)>q85(2));
data=data(ind_ext,:);
% data.reset_index(inplace=True,drop=True)
n = size(data,1);
% b is chosen such that ~1% of the data fall into the mooving average
b = floor(n/200); %#n/200: length to apply a moving average on 2b + 1 consecutive points
%%
TDC_col=zeros(n,1);
TDC_smoothed_col=zeros(n,1);
% data["TDC"] = np.NaN
% data["TDC_smoothed"] = np.NaN
if region == "upper"
%% Compute the Upper TDC for every possible threshold i/n
   for i=2:n-1 %Range(1,n-1):
       TDC_col(i,1) = TDC(data,(i/n),region);
   end
data=[data,TDC_col];
data =flip(data,1); %Reverse the order, the plateau finding algorithm starts with lower values
% data.reset_index(inplace=True,drop=True)
elseif  region =="lower"
%% Compute the Lower TDC for every possible threshold i/n
        for i =2:n-1 %range(1,n-1):
            TDC_col(i,1) = TDC(data,(i/n),region);
        end
        data=[data,TDC_col];
else
         print("Takes upper or lower argument as region only")
end
fprintf('TDC values are computed for all data\n');
%% Smooth the TDC with a mooving average of lenght 2*b+1
fprintf('Smooth the TDC with a mooving average of lenght 2*b+1\n');
for i=b+1:(n-b) %in range(b,n-b):
    TDC_smooth = 0 ;
    for j=(i-b):(i+b) %in range(i-b,i+b+1):
        TDC_smooth = TDC_smooth +data(j,3); %,"TDC"]
    end
    TDC_smooth = TDC_smooth/(2*b+1);
    TDC_smoothed_col(i) = TDC_smooth;
end
data=[data,TDC_smoothed_col];
m = floor(sqrt(n-2*b) ); %# lenght of the plateau
fprintf('computing the standard deviation of the smoothed series\n');
stdx = std(TDC_smoothed_col); %# the standard deviation of the smoothed series
% data["cond"] = np.NaN # column that will contains the condition: plateau - 2*sigma
len_k=n-2*b-m+1;
cond_col=ones(len_k,1);
for k=1:len_k %in range(0,n-2*b-m+1):
    plateau = 0 ;
    for i=(k+2):(k+m-1+1) %in range(k+1,k+m-1+1):
        plateau = plateau + abs(TDC_smoothed_col(i) - TDC_smoothed_col(k));
    end
    cond_col(k)= plateau - 2*stdx;
end
fprintf('Finding the first k such that: plateau - 2*sigma <= 0\n');
all_k = find(cond_col <= 0) ;
k = all_k(1) ;
%% The optimal TDC is defined as the average of the corresponding plateau
plateau_k = 0 ;
fprintf('optimal TDC is defined as the average of the corresponding plateau\n');
for i=2:m-1+1 %in range(1,m-1+1):
    plateau_k = plateau_k + data(k+i-1,4); 
end
tdc_value = plateau_k/m ;
kn=k/n;
end
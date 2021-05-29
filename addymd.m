% add year,month and day for a dataset, a year is supposed to have 365 days, no
% leap year is considered

function [X]=addymd(X)

y=zeros(size(X,1),3);
m=1;
n=1;
month=[0;31;59;90;120;151;181;212;243;273;304;334;365];
% add year
for i=365:365:length(X)
    y(m:i,1)=n;    
    m=i+1;
    n=n+1;
end

% add month
for i=1:length(X)/365
    for j=1:12
        y(((i-1)*365+month(j)+1):((i-1)*365+month(j+1)),2)=j;
    end
end

% add day
for i=1:length(X)/365
    for j=1:12
        k=find(y(:,1)==i&y(:,2)==j);
        y(k,3)=1:length(k);
    end
end

X=[y X];
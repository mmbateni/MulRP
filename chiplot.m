function  chiplot(u)
%chiplot : Chi-squared plot for independence
%chi. Returned values are in the interval (0,1).
%Translation of /CRAN/asbio/R/chi.plot.R 
%   References: Abberger, K. (2005). A simple graphical method to explore tail-dependence
%   in stock-return pairs. Applied Financial Economics, 15(1), 43-51. 
sz=18; %Plot Marker Size
u = u(prod(u,2)>0,:);
% to have illustrative figures a portion of random data points is plotted
u = u(floor(size(u,1)*0.70):floor(size(u,1)*0.775),:);
mid_y1=median(u(:,1));
mid_y2=median(u(:,2));
n = size(u, 1);
outer_1=zeros(n,n);
% outer_1 = sparse(n,n);
outer_2=zeros(n,n);
% outer_2 = sparse(n,n);
for i=1:n
    row_i=u(i,1)- u(:,1);
    outer_1(i,:)=row_i;
end
for j=1:n
    row_j=u(j,2)- u(:,2);
    outer_2(j,:)=row_j;
end 
outer_1(logical(eye(n)))=NaN;
outer_2(logical(eye(n)))=NaN;
t1 = outer_1 > 0;
t2 = outer_2 > 0;
t3 = t1 & t2;
%H =sum(t3,2,'omitnan')/(n - 1);
%Z =sum(t1,2,'omitnan')/(n - 1);
%G =sum(t2,2,'omitnan')/(n - 1);
H =nansum(t3,2)/(n - 1);
Z =nansum(t1,2)/(n - 1);
G =nansum(t2,2)/(n - 1);
S = sign((Z - 0.5).*(G - 0.5));
chi = (H - (Z .* G))./sqrt(Z .* (1 - Z) .* G .* (1 - G));
lambda = 4 * S .* max((Z - 0.5).^2, (G - 0.5).^2);
n = size(u, 1);
%% Top left: Scatter of the original data
subplot(2,2,1);%tiledlayout
%ax1 = nexttile;
xplot = u(:,1);
yplot = u(:,2);
scatter(eprob(xplot),eprob(yplot),sz);%'MarkerEdgeAlpha',.4
title('Scatter plot')
xlabel('Exceedance probability of Precipitation Monza') 
ylabel('Exceedance probability of Precipitation Lambro SUM') 
% Top right plot
subplot(2,2,2);%tiledlayout
hold on;
%ax2 = nexttile;
thresh = 4 * (1/(n - 1) - 0.5)^2;
refline_y=1.78/sqrt(n);
crit=[abs(lambda)<thresh,chi<0.642];
xplot = lambda(all(crit,2));
yplot = chi(all(crit,2));
scatter(xplot,yplot,sz);%'MarkerEdgeAlpha',1
line ("xdata",[-1,1], "ydata",[refline_y,refline_y], "linewidth", 1)
line ("xdata",[-1,1], "ydata",[-refline_y,-refline_y], "linewidth", 1)
%hline1 = refline([0 refline_y]);
%hline2 = refline([0 -refline_y]);
%hline1.Color = 'r';
%hline2.Color = 'r';
title('Chi Plot')
xlabel('λ') 
ylabel('χ') 
%% Left bottom plot
subplot(2,2,3);%tiledlayout
hold on;
%ax3 = nexttile;
u_ind_1q=and((u(:,1)<mid_y1),(u(:,2)<mid_y2));
u_1d=u(u_ind_1q,:);
n = size(u_1d, 1);
outer_1=zeros(n,n);
% outer_1 = sparse(n,n);
outer_2=zeros(n,n);
for i=1:n
    row_i=u_1d(i,1)- u_1d(:,1);
    outer_1(i,:)=row_i;
end
for j=1:n
    row_j=u_1d(j,2)- u_1d(:,2);
    outer_2(j,:)=row_j;
end 
outer_1(logical(eye(n)))=NaN;
outer_2(logical(eye(n)))=NaN;
t1 = outer_1 > 0;
t2 = outer_2 > 0;
t3 = t1 & t2;
%H =sum(t3,2,'omitnan')/(n - 1);
%Z =sum(t1,2,'omitnan')/(n - 1);
%G =sum(t2,2,'omitnan')/(n - 1);
H =nansum(t3,2)/(n - 1);
Z =nansum(t1,2)/(n - 1);
G =nansum(t2,2)/(n - 1);
S = sign((Z - 0.5).*(G - 0.5));
chi = (H - (Z .* G))./sqrt(Z .* (1 - Z) .* G .* (1 - G));
lambda = 4 * S .* max((Z - 0.5).^2, (G - 0.5).^2);
thresh = 4 * (1/(n - 1) - 0.5)^2;
refline_y=1.78/sqrt(n);
xplot = lambda(and(abs(lambda)<thresh,lambda>0));
yplot = chi(and(abs(lambda)<thresh,lambda>0));
scatter(xplot,yplot,sz);
line ("xdata",[0,1], "ydata",[refline_y,refline_y], "linewidth", 1)
line ("xdata",[0,1], "ydata",[-refline_y,-refline_y], "linewidth", 1)
%hline1 = refline([0 refline_y]);
%hline2 = refline([0 -refline_y]);
%hline1.Color = 'r';
%hline2.Color = 'r';
title('Chi Plot of 1st (lower-left) Quadrant')
xlabel('λ') 
ylabel('χ') 
%% Right bottom plot
subplot(2,2,4);%tiledlayout
hold on;
%ax4 = nexttile;
u_ind_4q=and((u(:,1)>mid_y1),(u(:,2)>mid_y2));
u_4d=u(u_ind_4q,:);
n = size(u_4d, 1);
outer_1=zeros(n,n);
% outer_1 = sparse(n,n);
outer_2=zeros(n,n);
for i=1:n
    row_i=u_4d(i,1)- u_4d(:,1);
    outer_1(i,:)=row_i;
end
for j=1:n
    row_j=u_4d(j,2)- u_4d(:,2);
    outer_2(j,:)=row_j;
end 
outer_1(logical(eye(n)))=NaN;
outer_2(logical(eye(n)))=NaN;
t1 = outer_1 > 0;
t2 = outer_2 > 0;
t3 = t1 & t2;
%H =sum(t3,2,'omitnan')/(n - 1);
%Z =sum(t1,2,'omitnan')/(n - 1);
%G =sum(t2,2,'omitnan')/(n - 1);
H =nansum(t3,2)/(n - 1);
Z =nansum(t1,2)/(n - 1);
G =nansum(t2,2)/(n - 1);
S = sign((Z - 0.5).*(G - 0.5));
chi = (H - (Z .* G))./sqrt(Z .* (1 - Z) .* G .* (1 - G));
lambda = 4 * S .* max((Z - 0.5).^2, (G - 0.5).^2);
thresh = 4 * (1/(n - 1) - 0.5)^2;
refline_y=1.78/sqrt(n);
crit=[abs(lambda)<thresh,lambda>0,chi>refline_y,chi<0.642];
xplot = lambda(all(crit,2));
yplot = chi(all(crit,2));
scatter(xplot,yplot,sz);
line ("xdata",[0,1], "ydata",[refline_y,refline_y], "linewidth", 1)
line ("xdata",[0,1], "ydata",[-refline_y,-refline_y], "linewidth", 1)
%hline1 = refline([0 refline_y]);
%hline2 = refline([0 -refline_y]);
%hline1.Color = 'r';
%hline2.Color = 'r';
title('Chi-Plot of 3rd (upper-right) Quadrant')
xlabel('λ') 
ylabel('χ') 
end

function [ x_s ] = disagg_knn( X,Zsimu,modindx )
%The sum of hourly rainfall Xd is daily rainfall Z. The sum of the total 24-h rainfall (X1, X2, : : : , Xd, d = 4) should
% be equal to the daily rainfall Z.To disaggregate daily rainfall data into hourly data, a set of n historical hourly rainfall data must
% be expressed in the form of a matrix X(d *n). the daily rainfall data vector, Z(1 *n), must be de?ned. We use the Gram–Schmidt
% process and its inverse process to disaggregate the Z value. The Gram–Schmidt process is an
% orthonormalizing column vector of X to the plane. The rotation matrix R is obtained using the method
% described by [Nowak, K., Prairie, J., Rajagopalan, B., & Lall, U. (2010). A nonparametric stochastic approach for multisite disaggregation of annual to daily streamflow. Water Resources Research, 46(8)].
% Using the rotation matrix R, Y(d *n) can be obtained as Y = R *X.
% Using matrix R, matrix X is rotated to matrix Y and the last row of the transformed matrix Y is
% changed to Yd = Z/sqrt(d).
if nargin<3
  modindx = 1;
end
n=size(X,1);
Ntot=size(X,2);
R= Rot_mat(n);
Y=R*X;
Zsim1=Zsimu./sqrt(n);
NN=round(sqrt(Ntot));
W3=(1:NN);
W4=1./W3;
Sum_K=sum(W4);
W5=Sum_K.*W3;
W=1./W5;
Y_end=Y(end,:);
for ij=1:length(Zsimu)
    Diffs=abs(Y_end-Zsim1(ij));
    [~,I] = sort(Diffs);
    I_k=I(1:NN);
    Y_NH=Y(1:(end-1),I_k);
    W_tot=repmat(W,(n-1),1);
    u=sum((W_tot.*Y_NH),2);
    y_s(:,ij)=[u;Zsim1(ij)];
end
x_s=R'*y_s;
% if modindx==1;
for ij=1:size(x_s,2)
    x_s_i=x_s(:,ij);
    k_n = find(x_s_i<0);
    k_p = find(x_s_i>0);
    add_fact=sum(x_s_i(k_n));
    contrib_tot=sum(x_s_i(k_p));
    add_e_fact=add_fact/contrib_tot;
    v_pos_add=add_e_fact.*x_s_i(k_p);
    x_s_i(k_p)=x_s_i(k_p)+v_pos_add;
    x_s_i(k_n)=0;
    x_s(:,ij)=x_s_i;
end
end


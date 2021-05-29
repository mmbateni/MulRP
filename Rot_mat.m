function [U3] = Rot_mat(n)
% Rotation Matrix required for disaggregation 
V=eye(n);
V(n,:)=sum(V,1);
U = zeros(n,n);
U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
for i = 2:n
    U(:,i) = V(:,i);
    for j = 1:i-1
        U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )*U(:,j);
    end
    U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
end
U2=U';
ones_V1=ones(n);
ones_V2=ones(n);
ones_V1(:,end)=(-1)*ones(n,1);
ones_V2(end,:)=(-1)*ones(1,n);
ones_V=ones_V1.*ones_V2;
U3=U2.*ones_V;
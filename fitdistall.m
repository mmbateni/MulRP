function [selected_dist,pd,AIC_cand_dist] = fitdistall(data,threshold,maxiter,Tolerance,data_season)
names_dist={'exponential','gamma','extreme value' ,'lognormal','normal','weibull','extended generalized pareto' };%'generalized extreme value' 'generalized pareto' 
candidate_dist={'NAN','NAN','NAN' ,'NAN','NAN','NAN','NAN' };
h=ones(1,6);
AIC_cand_dist=zeros(1,7);
[zeroa]=find(data<threshold);
XS_MW_nozero=data;XS_MW_nozero(zeroa)=[];
q=length(zeroa)/length(data);
      for ij=1:6
%         pd = fitdist(XS_MW_nozero,names_dist{ij});
          [Params,NLogL] = fitdist(XS_MW_nozero,names_dist{ij});
          NLL=NLogL;
          k=numel(Params);
          n=numel(XS_MW_nozero);
%         BIC(ij)=-2*(-NLL)+k*log(n);
          AIC1(ij)=-2*(-NLL)+2*k;
          AICc(ij)=AIC1(ij)+((2*k*(k+1))/(n-k-1));
          [pval,~]=kolmogorov_smirnov (XS_MW_nozero, names_dist{ij}, Params,"<>");
          if pval>0.90
             h(ij)=0
          end
%               h(ij)=kstest(XS_MW_nozero,'CDF',pd);
%               h_AD(ij)=adtest(XS_MW_nozero,'Distribution',pd);
          prob_dist_obj{ij}=Params;
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for extended generalized pareto
h(7)=0;
%h_AD(7)=0;
opt = optimset('Algorithm','nonlin_min','GradObj','off',...
   'Hessian','off',...
   'LargeScale','off',...
   'Display','off',...
   'MaxIter', maxiter,...
   'MaxFunEvals', 100000, ...
   'TolFun', Tolerance ,...
   'TolX', 1.000e-006,...
   'FunValCheck','off',...
   'DerivativeCheck','off',...
   'Diagnostics','off',...
   'Algorithm','interior-point');    
%parinit=[1,2.5,0.08];
parinit=[15,20,0.08];
[phat_value,~] = fmincon(@(parinit) egp_ll(data_season,threshold,parinit),parinit,[0,0,0;0,-1,0;0,0,-1],[0;-0.5;0],[],[],[],[],[],opt);
parinit_xi=phat_value(3);
parinit(3)=parinit_xi;
try
  [phat_value,NLL_egp] = fmincon(@(parinit) egp_ll(data,threshold,parinit),parinit,[0,0,0;0,-1,0;0,0,-1],[0;-0.5;0],[0,0,0;0,0,0;0,0,1],[0;0;parinit_xi],[],[],[],opt);
catch
  phat_value=NA;
  NLL_egp=NA;
  h(7)=1;
end_try_catch
k=3;
n=numel(data);
AIC1(7)=-2*(-NLL_egp)+2*k;
AICc(7)=AIC1(7)+((2*k*(k+1))/(n-k-1));
%prob_dist_obj{ij}=phat_value;
%pd9.DistName=names_dist{7};
prob_dist_obj{7}=phat_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ind_dist=find(~h&~h_AD);
ind_dist=find(~h);%if h=o, then the theoritical distribution is appropriate
candidate_dist(ind_dist)=names_dist(ind_dist);
AIC_cand_dist(ind_dist)=AICc(ind_dist);
if numel(ind_dist)>0
    [~,I_AIC]=max(AIC_cand_dist);
else
    [~,I_AIC]=max(AICc);
end
selected_dist=names_dist(I_AIC);
pd=prob_dist_obj{I_AIC};
AIC_cand_dist=AICc(I_AIC);
% Prob_xs=q+(1-q)*cdf(pd,XS_MW);
% Prob_xs(zeroa)=(length(zeroa)+1)/(2*(length(XS_MW)+1));
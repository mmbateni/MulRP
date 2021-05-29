clc;
clear all;
name_percip={'Abajaloo','Band','Dizaj'};
name_hyd={'Nazloo','Barandooz','Shahrchay'};
label={'min','max','average','lag-1AC','skewness','Stdeviation','Hexponent'};
name_sta=horzcat(name_percip,name_hyd);
X_percip_tot=zeros(40,1000,3);
X_hyd_tot=zeros(51,1000,3);
X_percip_Abajaloo=zeros(40,1000,12);
X_percip_Band=zeros(40,1000,12);
X_percip_Dizaj=zeros(40,1000,12);
X_hyd_Nazloo=zeros(51,1000,12);
X_hyd_Barandooz=zeros(51,1000,12);
X_hyd_Shahrchay=zeros(51,1000,12);
histo_data=Xlsreadall('Input_knn_rep.xlsx');
rep_knn_percip=histo_data(1);
rep_percip1=cell2mat(rep_knn_percip);
rep_percip2=rep_percip1(:,1:(end-1));
rep_percip=rep_percip2';
rep_knn_hyd=histo_data(2);
rep_hyd1=cell2mat(rep_knn_hyd); 
rep_hyd2=rep_hyd1(:,1:(end-1));
rep_hyd=rep_hyd2';
%%
histo_temp_data=Xlsreadall('Historical_input.xlsx');
rep_temp_Nazloo=histo_temp_data(4);
rep_temp_Barandooz=histo_temp_data(5);
rep_temp_Shahrchay=histo_temp_data(6);
rep_Nazloo1=cell2mat(rep_temp_Nazloo);
rep_Nazloo=rep_Nazloo1';
rep_Nazloo_t=[rep_Nazloo;sum(rep_Nazloo)];
rep_Barandooz1=cell2mat(rep_temp_Barandooz);
rep_Barandooz=rep_Barandooz1';
rep_Barandooz_t=[rep_Barandooz;sum(rep_Barandooz)];
rep_Shahrchay1=cell2mat(rep_temp_Shahrchay);
rep_Shahrchay=rep_Shahrchay1';
rep_Shahrchay_t=[rep_Shahrchay;sum(rep_Shahrchay)];
%
rep_temp_Abajaloo=histo_temp_data(1);
rep_temp_Band=histo_temp_data(2);
rep_temp_Dizaj=histo_temp_data(3);
rep_Abajaloo1=cell2mat(rep_temp_Abajaloo);
rep_Abajaloo=rep_Abajaloo1';
rep_Abajaloo_t=[rep_Abajaloo;sum(rep_Abajaloo)];
rep_Band1=cell2mat(rep_temp_Band);
rep_Band=rep_Band1';
rep_Band_t=[rep_Band;sum(rep_Band)];
rep_Dizaj1=cell2mat(rep_temp_Dizaj);
rep_Dizaj=rep_Dizaj1';
rep_Dizaj_t=[rep_Dizaj;sum(rep_Dizaj)];

%%
input_knn=Xlsreadall('Input_knn.xlsx');
input_knn_percip=input_knn(1);
data_percip=cell2mat(input_knn_percip);
[~,~,nz_data_percip] = find(data_percip);
mat_data_percip=reshape(nz_data_percip,40,[]);
input_knn_hyd=input_knn(2);
data_hyd=cell2mat(input_knn_hyd);
[~,~,nz_data_hyd] = find(data_hyd);
mat_data_hyd=reshape(nz_data_hyd,51,[]);
%%
for jp=1:length(nz_data_percip)
[ x_percip ] = disagg_knn( rep_percip,mat_data_percip(jp),1);
d3_ind=ceil(jp/40);
d2_ind=jp-((d3_ind-1)*40);
X_percip_tot(d2_ind,d3_ind,:)=x_percip;
end
for jh=1:length(nz_data_hyd)
[ x_hyd ] = disagg_knn( rep_hyd,mat_data_hyd(jh),1);
d3_ind_h=ceil(jh/51);
d2_ind_h=jh-((d3_ind_h-1)*51);
X_hyd_tot(d2_ind_h,d3_ind_h,:)=x_hyd;
end
%%
for kp=1:length(nz_data_percip)
d3_ind=ceil(kp/40);
d2_ind=kp-((d3_ind-1)*40);
[ x_percip1 ] = disagg_knn( rep_Abajaloo,X_percip_tot(d2_ind,d3_ind,1),1 );
[ x_percip2 ] = disagg_knn( rep_Band,X_percip_tot(d2_ind,d3_ind,2),1 );
[ x_percip3 ] = disagg_knn( rep_Dizaj,X_percip_tot(d2_ind,d3_ind,3),1 );
X_percip_Abajaloo(d2_ind,d3_ind,:)=x_percip1;
X_percip_Band(d2_ind,d3_ind,:)=x_percip2;
X_percip_Dizaj(d2_ind,d3_ind,:)=x_percip3;

end
for kh=1:length(nz_data_hyd)
d3_ind_h=ceil(kh/51);
d2_ind_h=kh-((d3_ind_h-1)*51);
[ x_hyd1 ] = disagg_knn( rep_Nazloo,X_hyd_tot(d2_ind_h,d3_ind_h,1),1 );
[ x_hyd2 ] = disagg_knn( rep_Barandooz,X_hyd_tot(d2_ind_h,d3_ind_h,2),1 );
[ x_hyd3 ] = disagg_knn( rep_Shahrchay,X_hyd_tot(d2_ind_h,d3_ind_h,3),1);
X_hyd_Nazloo(d2_ind_h,d3_ind_h,:)=x_hyd1;
X_hyd_Barandooz(d2_ind_h,d3_ind_h,:)=x_hyd2;
X_hyd_Shahrchay(d2_ind_h,d3_ind_h,:)=x_hyd3;
end
%%
X_tot_Nazloo = cat(3, X_hyd_Nazloo, X_hyd_tot(:,:,1));
X_tot_Barandooz = cat(3, X_hyd_Barandooz, X_hyd_tot(:,:,2));
X_tot_Shahrchay = cat(3, X_hyd_Shahrchay, X_hyd_tot(:,:,3));
X_tot_Abajaloo = cat(3, X_percip_Abajaloo, X_percip_tot(:,:,1));
X_tot_Band = cat(3, X_percip_Band, X_percip_tot(:,:,2));
X_tot_Dizaj = cat(3, X_percip_Dizaj, X_percip_tot(:,:,3));
%%
    naz_barandooz=zeros(1000,13);
    naz_shahr=zeros(1000,13);
    barandooz_shahr=zeros(1000,13);
    abajaloo_band=zeros(1000,13);
    abajaloo_dizaj=zeros(1000,13);
    band_dizaj=zeros(1000,13);
    n_nz=13*(13-1)/2;
    naz_naz=zeros(1000,n_nz);
    barandooz_barandooz=zeros(1000,n_nz);
    shahr_shahr=zeros(1000,n_nz);
    abajaloo_abajaloo=zeros(1000,n_nz);
    band_band=zeros(1000,n_nz);
    dizaj_dizaj=zeros(1000,n_nz);
 for jj=1:1000
    nazloo_corr1=X_tot_Nazloo(:,jj,:);
    nazloo_corr=reshape(nazloo_corr1,51,[]);
    barandooz_corr1=X_tot_Barandooz(:,jj,:);
    barandooz_corr=reshape(barandooz_corr1,51,[]);
    shahrchay_corr1=X_tot_Shahrchay(:,jj,:);
    shahrchay_corr=reshape(shahrchay_corr1,51,[]);
    abajaloo_corr1=X_tot_Abajaloo(:,jj,:);
    abajaloo_corr=reshape(abajaloo_corr1,40,[]);
    band_corr1=X_tot_Band(:,jj,:);
    band_corr=reshape(band_corr1,40,[]);
    dizaj_corr1=X_tot_Dizaj(:,jj,:);
    dizaj_corr=reshape(dizaj_corr1,40,[]);
    naz_barandooz(jj,:)=diag(corr(nazloo_corr,barandooz_corr))';
    naz_shahr(jj,:)=diag(corr(nazloo_corr,shahrchay_corr))';
    barandooz_shahr(jj,:)=diag(corr(barandooz_corr,shahrchay_corr))';
    abajaloo_band(jj,:)=diag(corr(abajaloo_corr,band_corr))';
    abajaloo_dizaj(jj,:)=diag(corr(abajaloo_corr,dizaj_corr))';
    band_dizaj(jj,:)=diag(corr(band_corr,dizaj_corr))';
    [~,~,naz_naz(jj,:)]=find(tril(corr(nazloo_corr,nazloo_corr),-1));
    [~,~,barandooz_barandooz(jj,:)]=find(tril(corr(barandooz_corr,barandooz_corr),-1));
    [~,~,shahr_shahr(jj,:)]=find(tril(corr(shahrchay_corr,shahrchay_corr),-1));
    [~,~,abajaloo_abajaloo(jj,:)]=find(tril(corr(abajaloo_corr,abajaloo_corr),-1));
    [~,~,band_band(jj,:)]=find(tril(corr(band_corr,band_corr),-1));
    [~,~,dizaj_dizaj(jj,:)]=find(tril(corr(dizaj_corr,dizaj_corr),-1));   
 end
    N_Nazloo = cell(13,1);
    N_Barandooz = cell(13,1);
    N_Shahrchay = cell(13,1);
    N_Abajaloo = cell(13,1);
    N_Band = cell(13,1);
    N_Dizaj = cell(13,1);
%     nbins=7;
%     edge_Nazloo_tot=zeros(13,nbins+1);
%     edge_Barandooz_tot=zeros(13,nbins+1);
%     edge_Shahrchay_tot=zeros(13,nbins+1);
%     edge_Abajaloo_tot=zeros(13,nbins+1);
%     edge_Band_tot=zeros(13,nbins+1);
%     edge_Dizaj_tot=zeros(13,nbins+1);
    edge_Nazloo_tot={};edge_Barandooz_tot={};edge_Shahrchay_tot={};edge_Abajaloo_tot={};edge_Band_tot={};edge_Dizaj_tot={};
for ji=1:13
    X_Nazloo_temp1=X_tot_Nazloo(:,:,ji);
    [~,edge_Nazloo] = histcounts(rep_Nazloo_t(ji,:));
    edge_Nazloo=[edge_Nazloo,inf];
    if edge_Nazloo(1)>0
        edge_Nazloo=[0,edge_Nazloo];
    end
    edge_Nazloo_tot{ji}=edge_Nazloo;
%     max_Nazloo_temp=max(max(X_Nazloo_temp1));
%     min_Nazloo_temp=min(min(X_Nazloo_temp1));
%     diff_Nazloo=((max_Nazloo_temp-min_Nazloo_temp)/7);
%     edge_Nazloo=min_Nazloo_temp:diff_Nazloo:max_Nazloo_temp;
    X_Barandooz_temp1=X_tot_Barandooz(:,:,ji);
    [~,edge_Barandooz] = histcounts(rep_Barandooz_t(ji,:));
    edge_Barandooz=[edge_Barandooz,inf];
    if edge_Barandooz(1)>0
        edge_Barandooz=[0,edge_Barandooz];
    end
    edge_Barandooz_tot{ji}=edge_Barandooz;
%     max_Barandooz_temp=max(max(X_Barandooz_temp1));
%     min_Barandooz_temp=min(min(X_Barandooz_temp1));
%     diff_Barandooz=((max_Barandooz_temp-min_Barandooz_temp)/7);
%     edge_Barandooz=min_Barandooz_temp:diff_Barandooz:max_Barandooz_temp;
    X_Shahrchay_temp1=X_tot_Shahrchay(:,:,ji);
    [~,edge_Shahrchay] = histcounts(rep_Shahrchay_t(ji,:));
    edge_Shahrchay=[edge_Shahrchay,inf];
    if edge_Shahrchay(1)>0
        edge_Shahrchay=[0,edge_Shahrchay];
    end
    edge_Shahrchay_tot{ji}=edge_Shahrchay;
%     max_Shahrchay_temp=max(max(X_Shahrchay_temp1));
%     min_Shahrchay_temp=min(min(X_Shahrchay_temp1));
%     diff_Shahrchay=((max_Shahrchay_temp-min_Shahrchay_temp)/7);
%     edge_Shahrchay=min_Shahrchay_temp:diff_Shahrchay:max_Shahrchay_temp;
    X_Abajaloo_temp1=X_tot_Abajaloo(:,:,ji);
    [~,edge_Abajaloo] = histcounts(rep_Abajaloo_t(ji,:));  
    edge_Abajaloo=[edge_Abajaloo,inf];
    if edge_Abajaloo(1)>0
        edge_Abajaloo=[0,edge_Abajaloo];
    end
    edge_Abajaloo_tot{ji}=edge_Abajaloo;
%   max_Abajaloo_temp=max(max(X_Abajaloo_temp1));
%   min_Abajaloo_temp=min(min(X_Abajaloo_temp1));
%   diff_Abajaloo=((max_Abajaloo_temp-min_Abajaloo_temp)/7);
%   edge_Abajaloo=min_Abajaloo_temp:diff_Abajaloo:max_Abajaloo_temp;
    X_Band_temp1=X_tot_Band(:,:,ji);
    [~,edge_Band] = histcounts(rep_Band_t(ji,:)); 
    edge_Band=[edge_Band,inf];
    if edge_Band(1)>0
        edge_Band=[0,edge_Band];
    end
    edge_Band_tot{ji}=edge_Band;
%     max_Band_temp=max(max(X_Band_temp1));
%     min_Band_temp=min(min(X_Band_temp1));
%     diff_Band=((max_Band_temp-min_Band_temp)/7);
%     edge_Band=min_Band_temp:diff_Band:max_Band_temp;
    X_Dizaj_temp1=X_tot_Dizaj(:,:,ji);
    [~,edge_Dizaj] = histcounts(rep_Dizaj_t(ji,:)); 
    edge_Dizaj=[edge_Dizaj,inf];
    if edge_Dizaj(1)>0
        edge_Dizaj=[0,edge_Dizaj];
    end
    edge_Dizaj_tot{ji}=edge_Dizaj;
%     max_Dizaj_temp=max(max(X_Dizaj_temp1));
%     min_Dizaj_temp=min(min(X_Dizaj_temp1));
%     diff_Dizaj=((max_Dizaj_temp-min_Dizaj_temp)/7);
%     edge_Dizaj=min_Dizaj_temp:diff_Dizaj:max_Dizaj_temp;
    [N_Nazloo1] = zeros(length(edge_Nazloo)-1,1000);
    [N_Barandooz1] = zeros(length(edge_Barandooz)-1,1000);
    [N_Shahrchay1] = zeros(length(edge_Shahrchay)-1,1000);
    [N_Abajaloo1] = zeros(length(edge_Abajaloo)-1,1000);
    [N_Band1] = zeros(length(edge_Band)-1,1000);
    [N_Dizaj1] = zeros(length(edge_Dizaj)-1,1000); 
% %     N_Nazloo=zeros(7,1000,13); 
% %     N_Barandooz=zeros(7,1000,13);
% %     N_Shahrchay=zeros(7,1000,13);
% %     N_Abajaloo=zeros(7,1000,13);
% %     N_Band=zeros(7,1000,13);
% %     N_Dizaj=zeros(7,1000,13);
    for jk=1:1000
    [N_Nazloo1(:,jk)] = (histcounts(X_Nazloo_temp1(:,jk),edge_Nazloo))./51;
    [N_Barandooz1(:,jk)] = (histcounts(X_Barandooz_temp1(:,jk),edge_Barandooz))./51;
    [N_Shahrchay1(:,jk)] = (histcounts(X_Shahrchay_temp1(:,jk),edge_Shahrchay))./51;
    [N_Abajaloo1(:,jk)] = (histcounts(X_Abajaloo_temp1(:,jk),edge_Abajaloo))./40;
    [N_Band1(:,jk)] = (histcounts(X_Band_temp1(:,jk),edge_Band))./40;
    [N_Dizaj1(:,jk)] = (histcounts(X_Dizaj_temp1(:,jk),edge_Dizaj))./40; 
    end
    N_Nazloo{ji} = ( N_Nazloo1);
    N_Barandooz {ji}= ( N_Barandooz1);
    N_Shahrchay {ji}= ( N_Shahrchay1);
    N_Abajaloo {ji}= ( N_Abajaloo1);
    N_Band {ji}= ( N_Band1);
    N_Dizaj{ji} = ( N_Dizaj1);
%     N_Nazloo(:,:,ji)=N_Nazloo1;
%     N_Barandooz(:,:,ji)=N_Barandooz1;
%     N_Shahrchay(:,:,ji)=N_Shahrchay1;
%     N_Abajaloo(:,:,ji)=N_Abajaloo1;
%     N_Band(:,:,ji)=N_Band1;
%     N_Dizaj(:,:,ji)=N_Dizaj1;
end
 %%
[ statdata_Abajaloo ] = boxplot_stat( X_tot_Abajaloo );
[ statdata_Band ] = boxplot_stat( X_tot_Band );
[ statdata_Dizaj ] = boxplot_stat( X_tot_Dizaj );
[ statdata_Nazloo ] = boxplot_stat( X_tot_Nazloo );
[ statdata_Barandooz ] = boxplot_stat( X_tot_Barandooz );
[ statdata_Shahrchay ] = boxplot_stat( X_tot_Shahrchay );




    

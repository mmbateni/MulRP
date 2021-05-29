function [amounts_gen,amounts_gen2]=month_separate(data,years_sim,nstations,length_month)
% separate the generated precipitation and temperature by station
% data is a structure of 4*1, which each elemet is a matrix of
% (length_season*n_simyear)*n_stations
%'amounts_gen' has  long elements for each station
%
%% combine seasonal precip
%
warning('off','all');
months={'jan' ;'feb'; 'mar'; 'apr';'may';'jun';'jul';'aug';'sep';'oct';'nov';'dec'};
 % dec
multiexp_preDEC(:,(1:nstations))=data(1,1).pre(1:years_sim*length_month(1,12),:);
 % jan
multiexp_preJAN(:,(1:nstations))=data(1,1).pre(years_sim*length_month(1,12)+1 ...
    :years_sim*length_month(1,12)+years_sim*length_month(1,1),:);
 % feb
multiexp_preFEB(:,(1:nstations))=data(1,1).pre(years_sim*length_month(1,12) ...
    +years_sim*length_month(1,1)+1:end,:);
% mar
multiexp_preMAR(:,(1:nstations))=data(1,2).pre(1:years_sim*length_month(1,3),:);
 % apr
multiexp_preAPR(:,(1:nstations))=data(1,2).pre(years_sim*length_month(1,3)+1 ...
    :years_sim*length_month(1,3)+years_sim*length_month(1,4),:);
% may
multiexp_preMAY(:,(1:nstations))=data(1,2).pre(years_sim*length_month(1,3) ...
    +years_sim*length_month(1,4)+1:end,:);
 % jun
multiexp_preJUN(:,(1:nstations))=data(1,3).pre(1:years_sim*length_month(1,6),:);
 % jul
multiexp_preJUL(:,(1:nstations))=data(1,3).pre(years_sim*length_month(1,6)+1 ...
    :years_sim*length_month(1,6)+years_sim*length_month(1,7),:);
 % aug
multiexp_preAUG(:,(1:nstations))=data(1,3).pre(years_sim*length_month(1,6) ...
    +years_sim*length_month(1,7)+1:end,:);
% sep
multiexp_preSEP(:,(1:nstations))=data(1,4).pre(1:years_sim*length_month(1,9),:);
 % oct
multiexp_preOCT(:,(1:nstations))=data(1,4).pre(years_sim*length_month(1,9)+1 ...
    :years_sim*length_month(1,9)+years_sim*length_month(1,10),:);
 %nor
multiexp_preNOV(:,(1:nstations))=data(1,4).pre(years_sim*length_month(1,9) ...
    +years_sim*length_month(1,10)+1:end,:);
amounts_gen.month=months;
amounts_gen.val={multiexp_preJAN ;multiexp_preFEB; multiexp_preMAR; multiexp_preAPR;multiexp_preMAY;multiexp_preJUN;multiexp_preJUL;multiexp_preAUG;multiexp_preSEP;multiexp_preOCT;multiexp_preNOV;multiexp_preDEC};
amounts_gen2{1}=multiexp_preJAN;amounts_gen2{2}=multiexp_preFEB;amounts_gen2{3}=multiexp_preMAR;
amounts_gen2{4}=multiexp_preAPR;amounts_gen2{5}=multiexp_preMAY;amounts_gen2{6}=multiexp_preJUN;
amounts_gen2{7}=multiexp_preJUL;amounts_gen2{8}=multiexp_preAUG;amounts_gen2{9}=multiexp_preSEP;
amounts_gen2{10}=multiexp_preOCT;amounts_gen2{11}=multiexp_preNOV;amounts_gen2{12}=multiexp_preDEC;

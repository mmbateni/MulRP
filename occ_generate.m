%Challenge:
% This function make the synthtic sequence of occurences

function [occurrences]=occ_generate(trans_prob,normal_rand,length_month)
n=size(normal_rand,1);
%%%
% occ=zeros(n,1);
% for i=2:n
%                 if occ(i-1)==0 && (normal_rand(i)<= trans_prob(1))
%                     occ(i)=0;
%                 elseif occ(i-1)==1 && (normal_rand(i)<= trans_prob(2))
%                     occ(i)=0;
%                 else
%                     occ(i)=1;
%                 end
% end
%  occurrences=occ;
%%%
	occ=zeros(length_month,(n/length_month));
        for k=1:(n/length_month)
        for j=2:length_month	
            j_k=j+((k-1)*(length_month));
                if occ(j-1,k)==0 && (normal_rand(j_k)<= trans_prob(1))
                    occ(j,k)=0;
                elseif occ(j-1,k)==1 && (normal_rand(j_k)<= trans_prob(2))
                    occ(j,k)=0;
                else
                    occ(j,k)=1;
                end
        end
        end
  occurrences=reshape(occ,n,1);
  %%%
% occ=reshape(occ,length_month,[]);
% occ=occ';
% [p00_y,p10_y,total_y] = transition(occ)

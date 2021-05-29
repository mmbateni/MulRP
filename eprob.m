function [X_eprob] = eprob(X)
% eprob: calculates the exceedance probability for n column vectors in the
%        array [m n] X, where m are the observations. The probability is 
%        output in percent. eX is output.
%
% Usage: X_eprob = eprob(X);
%
% Input Arguments:
%
%   X - [m n] vector where m are the observations and n are the number of
%   datasets for which the exceedance probability is to be calculated. 
%   The size of m must be the same for all datasets.
%
% Output Arguments:
%
%   x - input data X [m 1]
%   x.r - the number of rows, m
%   x.sort - X input data sorted in descending order
%   x.rank - single column matrix of the sorted data rank
%   x.eprob - calculated exceedance probability (rank/m+1)
%
% Example:
%
%   X = randn(1000,1) % create randomly distributed dataset with 1000 values
%   X_eprob = eprob(X);
%

%Scap = 10; % active operational energy storage capacity
% Scap = StorCapPercent eX average annual generation
%eX = struct;
%X_data = X;
X_r = size(X,1); % no. rows
[~,I] = sort(X,'ascend'); % sorts data in descending order
[~,II] =sort(I,'ascend');
X_rank = (1:X_r)';
I_rank=X_rank(II);
%X_eprob = zeros(X_r,1);
X_eprob = I_rank./(X_r+1);

% % plotting eeXceedance probability curve (in %)
% plot(X_eprob * 100,X_sort,'r-','LineWidth',2);
% xlabel('Exceedance Probability (%)','FontWeight','Bold');
% ylabel('Value','FontWeight','Bold');
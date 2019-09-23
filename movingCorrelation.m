% Author: Michael J. Bommarito II
% Contact: michael.bommarito@gmail.com
% Date: Oct 3, 2010, Edited by Willem Huiskamp, 2015
% Provided as-is, informational purposes only, public domain.
%
% Inputs:
%   1. dataMatrix: variables (X1,X2,...,X_M) in the columns
%   2. windowSize: number of samples to include in the moving window
%   3. indexColumn: the variable X_i against which correlations should be 
% returned
%
% Output:
%   1. correlationTS: correlation between X_{indexColumn} and X_j for j !=
% indexColumn from windowSize+1 to the number of observations.  The first
% windowSize rows are NaN.
%   2. PvalTS: P values for each correlation
%   3. correlationTS2: Elements of correlationTS where P < 0.05

function [correlationTS correlationTS2 PvalTS] = movingCorrelation(dataMatrix, windowSize, indexColumn)

[N,M] = size(dataMatrix);
correlationTS = nan(N, M-1);
PvalTS = nan(N, M-1);

for t = windowSize:N
    [C P]= corrcoef(dataMatrix(t-windowSize+1:t, :));
    idx = setdiff(1:M, [indexColumn]);
    correlationTS(t-floor(windowSize/2), :) = C(indexColumn, idx);
    PvalTS(t-floor(windowSize/2), :) = P(indexColumn, idx);
end


[row,col] = find(PvalTS>0.05);
correlationTS2 = correlationTS;

for i = 1:length(col)
    correlationTS2(row(i),col(i)) = NaN;
end

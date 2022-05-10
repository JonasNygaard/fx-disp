function csum = nancumsum(data,dim)

%% nan_sum.m
% ########################################################################### %
% function  csum = nancumsum(data,dim)
% Purpose:  Compute cumulative sum while ignoring NaN entries 
%
% Input:    data    = TxN matrix of data
%           dim     = Scalar inditing whether to sum over columns (1) or rows (2)
%
% Output:   csum    = Cumulative Sum of input matrix
%               
% Author:
% Jonas N. Eriksen
% Department of Economics and Business Economics
% Aarhus University and CREATES
%
% Encoding: UTF8
% Last modified: April, 2018
% ########################################################################### %

% Error checking
if isempty(data)
    error('nancumsum.m: Input matrix of portfolio returns is empty');
end

% Setting default dimension
if nargin < 2
    dim = 1;
end

%% Computing the cumulative sum while treating NaNs as zeros
% ########################################################################### %
%{
    We first identify any NaN entries, then we set them to zero so as not
    to influence the sum and then add back NaN for all NaN columns/rows.  
%}
% ########################################################################### %

csum            = data;
nan_indx        = find(csum~=csum);
csum(nan_indx)  = 0;
csum            = cumsum(csum,dim);
csum(nan_indx)  = NaN;

end

% ########################################################################### %
% [EOF]
% ########################################################################### %
function sumStats = portfolioStatistics(returns,annualize)

%% linearGMM.m
% ########################################################################### %
% function  sumStats = portfolioStatistics(returns,annualize)
% Purpose:  Compute standard statistics for a series of returns
%
% Input:    returns     = TxN matrix of portfolio excess returns
%           annualize   = Flag for annualizing (true) or not (false)
%
% Output:   sumStats    = Matrix of output results
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
if isempty(returns)
    error('portfolioStatistics.m: Portfolio return matrix is empty');
end

if (nargin > 2)
    error('portfolioStatistics.m: To many input arguments');
end

if (nargin == 2) && ~ismember(annualize,{'true','false'})
    error('portfolioStatistics.m: Wrong argument for annualizing returns');
end

% Setting default
if (nargin == 1)
    annualize = 'true';
end

%% Compute standard summary statistics
% ########################################################################### %
%{
    We compute standard summary statistics including means, standard 
    deviations, Sharpe ratio, skewness, kurtosis, and first-order
    autocorrelation coefficients.  
%}
% ########################################################################### %

% Getting dimension of portfolio return data
nAsset = size(returns,2);

% Compute descriptive statistics for portfolio returns
if strcmp(annualize,'true')

    meanReturns     = nanmean(returns).*12;
    stdReturns      = nanstd(returns)*sqrt(12);
    skewReturns     = skewness(returns);
    kurtReturns     = kurtosis(returns);
    sharpeReturns   = meanReturns./stdReturns;

elseif strcmp(annualize,'false')

    meanReturns     = nanmean(returns);
    stdReturns      = nanstd(returns);
    skewReturns     = skewness(returns);
    kurtReturns     = kurtosis(returns);
    sharpeReturns   = meanReturns./stdReturns;

end

% Computing t-statistic for mean return for each asset/portfolio
tstatReturns    = NaN(nAsset,1);

for iAsset = 1:nAsset

    % Remove row with NaN entries for regression
    pfret_loop = returns(~any(isnan(returns(:,iAsset)),2),iAsset);

    % Compute t-statistic using robust standard errors
    res = linRegNWA(pfret_loop,ones(size(pfret_loop,1),1),0);
    tstatReturns(iAsset,1) = res.t;

end

% Collect summary statistics in output matrix
sumStats(1,:)   = meanReturns;
sumStats(2,:)   = tstatReturns;
sumStats(3,:)   = stdReturns;
sumStats(4,:)   = sharpeReturns;
sumStats(5,:)   = skewReturns;
sumStats(6,:)   = kurtReturns;

end

% ########################################################################### %
% [EOF]
% ########################################################################### %
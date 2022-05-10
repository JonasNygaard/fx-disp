function regResults = linRegNWA(y,x,const)

%% linRegNWA.m
% ########################################################################### %
% function  function regResults = linRegNWA(y,x,const)
% Purpose:  Estimate linear regression model and compute HAC standard errors
%           using Newey and West (1987) with Andrews (1991) lag selection.
%
% Input:    y       = T x N matrix of dependent variables
%           x       = T x K matrix of independent variables
%           const   = Scalar for adding/omitting intercept in the regression
%               - Default is no cross-sectional intercept (const = 1)
%
% Output:   regResults  = struct with model characteristics
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
if size(x,1) ~= size(y,1)
  error('linRegNWA.m: Unequal number of time series observations in y and x'); 
end

if (nargin < 2)
    error('linRegNWA.m: Not enough input parameters')
end

if (nargin > 3)
    error('linRegNWA.m: Too many input parameters');
end

%Setting constant as default
if (~isequal(ones(size(x,1),1),x(:,1))) && (nargin == 2)
    const = 1;
end

%% Setting preliminaries
% ########################################################################### %
%{
    Adding constant to design matrix per default and obtaining data
    dimensions. 
%}
% ########################################################################### %

% Adding constant to design matrix
if const == 1
    x = [ones(size(x,1),1) x];
end

% Getting data dimensions
[nObs,nReg] = size(y);
nVar        = size(x,2);

%% Estimating model coefficients and goodness of fit measures
% ########################################################################### %
%{
    We estimate the model parameters and compute goodness of fit measures.
%}
% ########################################################################### %

% Estimating X'X matrix, coefficients, and residuals
Exx         = x'*x/nObs;
coef        = x\y;
errMat      = y - x*coef;

% Computing r-squares and setting preliminaries for HAC constructions
s2          = mean(errMat.^2);
vary        = mean((y - ones(nObs,1) * mean(y)).^2);
rSquared    = 100.*(1 - s2./vary)';
adjRsquared = 100.*(1 - (s2./vary)*(nObs-1)/(nObs-nVar))';

%% Computing Newey-West (1987) errors with Andrews (1991) lag selection
% ########################################################################### %
%{
    We compute HAC standard errors following Newey and West (1987) using a 
    Bartlett kernel and automatic bandwidth selection using the method
    outlines in Andrews (1991). 
%}
% ########################################################################### %

% Preallocations
stdError    = zeros(nVar,nReg);
tStat       = zeros(nVar,nReg);
fTest       = zeros(nVar,3);

% Computing Newey-West standard errors with Andrews bandwidth selection
for iReg = 1:nReg

    % Picking out residuals for model i  
    errVec  = errMat(:,iReg);

    % Computing moment conditions errors
    u       = x.*(errVec*ones(1,nVar));

    % Computing Newey-West-Andrews covariance matrix 
    S       = hacAndrews(u);
    varCoef = Exx\S/Exx/nObs;

    % F-test on coefficient being zero (except constant)
    chi2val         = coef(2:end,iReg)'*((varCoef(2:end,2:end))\coef(2:end,iReg));
    df              = size(coef(2:end,1),1);
    pval            = 1 - cdf('chi2',chi2val,df);
    fTest(iReg,1:3) = [chi2val df pval];
      
    % Computing standard errors and t-statistics
    stdError(:,iReg)    = sqrt(diag(varCoef));
    tStat(:,iReg)       = coef(:,iReg)./stdError(:,iReg);

end

% Building struct with output results
regResults.b    = coef;
regResults.se   = stdError;
regResults.t    = tStat;
regResults.r2   = rSquared;
regResults.ar2  = adjRsquared;
regResults.f    = fTest;

end

% ########################################################################### %
% [EOF]
% ########################################################################### %
function [firstStage,secondStage] = famaMacBeth(pfReturns,riskFactors,xscons)

%% famaMacBeth.m
% ########################################################################### %
% function  [firstStage,secondStage] = famaMacBeth(pfReturns,riskFactors,xscons)
% Purpose:  Estimate linear asset pricing models using the Fama and MacBeth
%           (1973) two-pass cross-sectional regression methodology.
%
% Input:    pfReturns   = TxN maxtrix of portfolio excess returns
%           riskFactors = TxK matrix of common risk factors
%           xscons      = Flag for adding/omitting intercept in CS regression
%               - Default is no cross-sectional intercept (xscons = "false")
%
% Output:   Two structures including results from the two steps
%
% Ref:      The function is based on Cochrane (2005), and Burnside (2011)
%               
% Author:
% Jonas N. Eriksen
% Department of Economics and Business Economics
% Aarhus University and CREATES
%
% Encoding: UTF8
% Last modified: April, 2018
% ########################################################################### %

% Error checking on input parameters
if (nargin < 2)
    error('famaMacBeth.m: Not enough input parameters');
end

if (nargin > 3)
    error('famaMacBeth.m: Too many input parameters');
end

if (size(pfReturns,1) ~= size(riskFactors,1))
    error('famaMacBeth.m: Unequal number of time series observations');
end

if (nargin == 3) && ~ismember(xscons,{'false','true'})
    error('famaMacBeth.m: xscons should either be "false" or "true"');
end

% Setting defaults
if (nargin == 2)
    xscons   = 'false';
end

%% Setting preliminaries
% ########################################################################### %
%{
    Obtain data dimensions for returns and risk factors, and construct
    vector of ones for use in first-step regression
%}
% ########################################################################### %

% Getting data dimensions
[nObs,nSeries]  = size(pfReturns);
nFactors        = size(riskFactors,2);
iota            = ones(nObs,1);

%% First-pass time series regressions
% ########################################################################### %
%{
    We rub time series regressions for each portfolio/asset of the kind
            r_{i,t}-r_{f,t} = 𝛼_{i} + β_{i}F_{t} + ε_{i,t}
    where we always include a constant in the regression.
%}
% ########################################################################### %

% Estimating factor betas (with constant per default)
tsReg               = linRegNWA(pfReturns,riskFactors,1);

% Estimating residual matrix for later use
errMat              = pfReturns - [iota riskFactors]*tsReg.b;

% Constructing struct with results
firstStage.beta     = tsReg.b;
firstStage.tstat    = tsReg.t;
firstStage.r2       = tsReg.r2;

%% Second-pass cross-sectional regressions
% ########################################################################### %
%{
    We run cross-sectional regressions for each time t = 1,...,T of the kind
            r_{i,t}-r_{f,t} = λ'β_{i} + u_{i,t}
    with or without an intercept as specified by the user. 
%}
% ########################################################################### %

if strcmp(xscons,'false')

    % Preallocations
    tsLambda    = NaN(nObs,nFactors);

    % Estimating lambda for each time period (assuming constant factor betas)
    for iObs = 1:nObs

        tsLambda(iObs,:)  = tsReg.b(2:end,:)'\pfReturns(iObs,:)';

    end

    % Estimating risk prices and Fama-MacBeth t-statistics
    riskPrices      = mean(tsLambda);
    seLambda        = sqrt(sum((tsLambda-repmat(riskPrices,nObs,1)).^2)/nObs^2);
    tLambdaFM       = riskPrices./seLambda;

    % Computing fitted values and R2 
    meanReturns     = mean(pfReturns);
    fittedValues    = tsReg.b(2:end,:)' * riskPrices';
    errResid        = mean(pfReturns)' - fittedValues;
    s2              = mean(errResid.^2);
    vary            = mean((meanReturns - ones(1,nSeries) * mean(meanReturns)).^2);
    rSquared        = 100*(1-s2./vary)';

elseif strcmp(xscons,'true')

    % Preallocations
    tsLambda    = NaN(nObs,nFactors+1);

    % Estimating lambda for each time period (assuming constant factor betas)
    for iObs = 1:nObs

        tsLambda(iObs,:)  = [ones(nSeries,1) tsReg.b(2:end,:)']\pfReturns(iObs,:)';

    end

    % Estimating risk prices and Fama-MacBeth t-statistics
    riskPrices      = mean(tsLambda);
    seLambda        = sqrt(sum((tsLambda-repmat(riskPrices,nObs,1)).^2)/nObs^2);
    tLambdaFM       = riskPrices./seLambda;

    % Computing fitted values and R2 
    meanReturns     = mean(pfReturns);
    fittedValues    = [ones(nSeries,1) tsReg.b(2:end,:)']*riskPrices';
    errResid        = mean(pfReturns)' - fittedValues;
    s2              = mean(errResid.^2);
    vary            = mean((meanReturns - ones(1,nSeries) * mean(meanReturns)).^2);
    rSquared        = 100*(1-s2./vary)';

end

%% Estimating GMM standard errors using Newey-West + Andrews method
% ########################################################################### %
%{
    We consider the following set of moment conditions as in Cochrane (2005)
                ε           0
        g = E[  εf  ]  = [  0  ]
                u           0      
%}
% ########################################################################### %

if strcmp(xscons,'false')

    % Setting up first set of moment conditions (time series regressions)
    u1          = (errMat(:,1)*ones(1,nFactors+1)).*[iota riskFactors];
    for iFactor = 2:nSeries
        u1  = [u1 (errMat(:,iFactor)*ones(1,nFactors+1)).*[iota riskFactors]];
    end

    % Setting up first set of moment conditions (cross-sectional regression)
    u2          = pfReturns - ones(nObs,1)*(tsReg.b(2:end,:)'*riskPrices')';

    % Estimating the spectral density matrix
    S           = hacAndrews([u1 u2]);

    % Estimating the covariance matrix for risk prices
    D11         = -kron(eye(nSeries),([iota riskFactors]'*[iota riskFactors])/nObs);
    D12         = zeros(nSeries*(nFactors+1),nFactors);
    D21         = -kron(eye(nSeries),[0 riskPrices]);
    D22         = -tsReg.b(2:end,:)';
    D           = [D11 D12 ; D21 D22 ];
    A           = [eye(nSeries*(nFactors+1)) zeros(nSeries*(nFactors+1),nSeries) ; ...
                        zeros(nFactors,nSeries*(nFactors+1)) (-D22')];
    gmmCov      = (A*D)\(A*S*(A'))/(D'*(A'))/nObs; %inv(A*D)*(A*S*(A'))*inv(D'*(A'))/nObs
    idx         = nSeries*(nFactors+1)+1:nSeries*(nFactors+1)+nFactors;
    lambdaCov   = gmmCov(idx,idx);

    % Computing GMM standard errors and t-statistics
    seLambda    = sqrt(diag(lambdaCov))';
    tLambda     = riskPrices./seLambda;

    % Testing for overidentifying restrictions
    nErrors     = size([u1 u2],2);
    mMat        = eye(nErrors) - D*((A*D)\A); % eye(nErrors) - D*inv(A*D)*A
    Vg          = mMat*S*mMat'/nObs;
    g           = mean([u1 u2])';
    chiValue    = g'*pinv(Vg)*g;
    pvalChi     = 1-chi2cdf(chiValue,nSeries-nFactors);

elseif strcmp(xscons,'true')

    % Setting up first set of moment conditions (time series regressions)
    u1          = (errMat(:,1)*ones(1,nFactors+1)).*[iota riskFactors];
    for iFactor = 2:nSeries
        u1  = [u1 (errMat(:,iFactor)*ones(1,nFactors+1)).*[iota riskFactors]];
    end

    u1alt       = kron(errMat,ones(1,nFactors+1)).*kron(ones(1,nSeries),[iota riskFactors]);

    % Setting up first set of moment conditions (cross-sectional regression)
    u2          = pfReturns - ones(nObs,1)*([ones(nSeries,1) tsReg.b(2:end,:)']*riskPrices')';

    % Estimating the spectral density matrix
    S           = hacAndrews([u1 u2]);

    % Estimating the covariance matrix for risk prices
    D11         = -kron(eye(nSeries),([iota riskFactors]'*[iota riskFactors])/nObs);
    D12         = zeros(nSeries*(nFactors+1),nFactors+1);
    D21         = -kron(eye(nSeries),[0 riskPrices(2:end)]);
    D22         = -[ones(nSeries,1) tsReg.b(2:end,:)'];
    D           = [ D11 D12 ; D21 D22 ];
    A           = [eye(nSeries*(nFactors+1)) zeros(nSeries*(nFactors+1),nSeries) ; ...
                        zeros(nFactors+1,nSeries*(nFactors+1)) (-D22')];
    gmmCov      = (A*D)\(A*S*(A'))/(D'*(A'))/nObs; %inv(A*D)*(A*S*(A'))*inv(D'*(A'))/nObs
    idx         = nSeries*(nFactors+1)+1:nSeries*(nFactors+1)+nFactors+1;
    lambdaCov   = gmmCov(idx,idx);

    % Computing GMM standard errors and t-statistics
    seLambda    = sqrt(diag(lambdaCov))';
    tLambda     = riskPrices./seLambda;

    % Testing for overidentifying restrictions
    nErrors     = size([u1 u2],2);
    mMat        = eye(nErrors) - D*((A*D)\A); % eye(nErrors) - D*inv(A*D)*A
    Vg          = mMat*S*mMat'/nObs;
    g           = mean([u1 u2])';
    chiValue    = g'*pinv(Vg)*g;
    pvalChi     = 1-chi2cdf(chiValue,nSeries-nFactors);

end

% Constructing structure with results
secondStage.lambda  = riskPrices;
secondStage.tstat   = tLambda;
secondStage.tstatFM = tLambdaFM;
secondStage.r2      = rSquared;
secondStage.chi2    = [chiValue pvalChi];
secondStage.fit     = fittedValues;
secondStage.mean    = meanReturns';

end

% ########################################################################### %
% [EOF]
% ########################################################################### %
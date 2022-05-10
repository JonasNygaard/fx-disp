function S = hacAndrews(data)

%% hacAndrew.m
% ########################################################################### %
% function  S = hacAndrews(errMat)
% Purpose:  Estimate the SDF parameters of linear asset pricing models using 
%           Hansen's (1982) Generalized Method of Moments (GMM). 
%
% Input:    pfReturns   = TxN maxtrix of portfolio excess returns
%           riskFactors = TxK matrix of common risk factors
%
% Output:   Structures including estimation results and implied risk prices
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

% Error checking on input
if (nargin > 1)
    error('hacAndrews.m: Too many input arguments');
end

if (nnz(isnan(data)) ~= 0)
    error('hacAndrews.m: Data contains NaN entries');
end

%% Computing bandwidth using Andrews (1991) automatic selection
% ########################################################################### %
%{
    Optimal bandwidth = 1.1447 * ( 𝛼(1) * T ) ^ (1/3)
%}
% ########################################################################### %

% Getting data dimensions
[nObs,nElem]    = size(data);
iota            = ones(nObs-1,1); 

% Estimate 𝛼(1) elements for each column in the data matrix
for iElem = 1:nElem

    % Estimating AR(1) coefficient for each components 
    rho             = [iota data(1:end-1,iElem)]\data(2:end,iElem);
    sigma2          = mean(([iota data(1:end-1,iElem)]*rho - data(2:end,iElem)).^2);
    alphaN(iElem)   = ( 4*rho(2)^2 * sigma2^2 / ((1-rho(2))^6 * (1+rho(2))^2) );
    alphaD(iElem)   = ( sigma2^2 / (1-rho(2))^4 );

end

% Estimate 𝛼(1) as equal-weighted average as in Andrews (1991)
alphaParm   = sum(alphaN)/sum(alphaD);

% Estimating the bandwidth parameter (ceil to get even number, conservative)
bandWidth   = ceil( 1.1447*(alphaParm*nObs)^(1/3) );

%% Estimate covariance matrix using Bartlett kernel (Newey and West, 1987)
% ########################################################################### %
%{
    We estimate the HAC covariance matrix using the Bartlett kernel following
    Newey and West (1987) and determine the optimal bandwidth using the 
    data-driven method of Andrews (1991). 
%}
% ########################################################################### %

% Specify Bartlett kernel weights
bartw   = ( bandWidth + 1- (0:bandWidth) ) ./ ( bandWidth+1 );

% Initialize covariance matrix
V       = data'*data/nObs;

% Estimate addition for each lag
for iLag = 1:bandWidth
    Gammai      = (data((iLag+1):nObs,:)'*data(1:nObs-iLag,:))/nObs;
    GplusGprime = Gammai+Gammai';
    V           = V + bartw(iLag+1)*GplusGprime;
end

% Final estimate
S = V;

% ########################################################################### %
% [EOF]
% ########################################################################### %
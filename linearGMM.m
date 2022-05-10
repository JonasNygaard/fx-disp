function gmmResults = linearGMM(pfReturns,riskFactors)

%% linearGMM.m
% ########################################################################### %
% function  gmmResults = linearGMM(pfReturns,riskFactors)
% Purpose:  Estimate the SDF parameters of linear asset pricing models using 
%           Hansen's (1982) Generalized Method of Moments (GMM). 
%
% Input:    pfReturns   = TxN matrix of portfolio excess returns
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

% Error checking on input parameters
if (nargin < 2)
    error('linearGMM.m: Not enough input parameters');
end

if (nargin > 2)
    error('linearGMM.m: Too many input parameters');
end

if (size(pfReturns,1) ~= size(riskFactors,1))
    error('linearGMM.m: Unequal number of time series observations');
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

% Computing number of elements, nuisance parameters, and degrees of freedom
nElem           = nFactors*(nFactors+1)/2;
nNuisance       = nFactors+nElem;
freedomDegrees  = nSeries-nFactors;

%% Estimating mean and covariances
% ########################################################################### %
%{
    We stack risk factors and excess returns in a matrix and compute
    means and a joint covariance matrix. This is similar to Kan, Robotti, 
    and Shanken (2013): "Pricing Model Performance and the Two-Pass 
    Cross-Sectional Regression Methodology". Note that it is the same that
    Burnside (2011) and Cochrane (2005) are using as d = V12. 
%}
% ########################################################################### %

% Computing means and covariance matrix of returns and factors
Y               = [riskFactors pfReturns];
mu1             = mean(riskFactors)';
mu2             = mean(pfReturns)';
mumu2           = mean(mu2'); 
covarY          = cov(Y,1);

% Picking out blocks of the covariance matrix
V11             = covarY(1:nFactors,1:nFactors);
V12             = covarY(1:nFactors,nFactors+1:end);
V21             = V12';
V22             = covarY(nFactors+1:end,nFactors+1:end);

%% Estimating SDF parameters and GMM Newey-West + Andrews standard errors
% ########################################################################### %
%{
    We consider a simple stochastic discount factor (SDF) given as
        m = 1 - bf
    and estimate the parameters of the model using GMM of Hansen (1982) as
        b = (d'd)^-1 * d'R
    where d=V21 is the second moment matrix of returns and factors. We refer
    to Burnside (2011) and Cochrane (2005) for further details. 
%}
% ########################################################################### %

% Estimating SDF parameters
b       = (V12*V21)\(V12*mu2);

% Building moment matrix for estimating error covariance matrix
u2      = (riskFactors-ones(nObs,1)*mu1');
sdf     = ones(nObs,1)-u2*b;
u1      = pfReturns.*(sdf*ones(1,nSeries));

% Building moments for factor covariance matrix (needed later)
idx     = tril(true(nFactors));
for iObs = 1:nObs

    tmp         = u2(iObs,:)'*u2(iObs,:);
    u3(iObs,:)  = tmp(idx) - V11(idx);

end

% Estimating the spectral density matrix
S       = hacAndrews([u1 u2]);

% Estimating the covariance matrix for b (see eq. 19, page page 13 in Burnside)
P       = [eye(nSeries) V21*b*(b')];
Vb      = inv(V21'*V21)*(V21'*P*S*(P')*V21)*inv(V21'*V21)./nObs;

% Getting standard errors and t-statistics for SDF parameters
seb     = sqrt(diag(Vb));
btstat  = b./seb;

%% Estimating risk prices and GMM Newey-West + Andrews standard errors
% ########################################################################### %
%{
    We estimate implied risk prices using the relation (Cochrane, 2005)
        ∑(f'f)*b = λ
    and estimate standard errors using the Delta method. 
%}
% ########################################################################### %

% Computing risk price estimates
lambda  = V11*b;

% Estimating the spectral density matrix
S       = hacAndrews([u1 u2 u3]);

% % Building matrices determining moments set to zero and determining weights
% a       = [V12*eye(nSeries) zeros(nFactors,nFactors) zeros(nFactors,nElem) ;...
%         zeros(nFactors,nSeries) eye(nFactors) zeros(nFactors,nElem) ;...
%         zeros(nElem,nSeries) zeros(nElem,nFactors) eye(nElem)];
% d       = [-V21 mu2*b' zeros(nSeries,nElem) ; ...
%         zeros(nFactors,nFactors) -eye(nFactors) zeros(nFactors,nElem) ; ...
%         zeros(nElem,nFactors) zeros(nElem,nFactors) -eye(nElem)];
% thetaCov    = inv(a*d)*a*S*a'*inv(d'*a')./nObs;

% Estimating the variance of theta (parameter vector)
Q           = [inv(V21'*V21)*(V21'*P) zeros(nFactors,nElem) ; ...
            zeros(nNuisance,nSeries) eye(nNuisance)]; % Q = inv(a*d)*a
thetaCov    = Q*S*(Q')./nObs;

% Estimating the standard errors of lambda via Delta method
[idx1,idx2] = find(tril(true(nFactors)) == 1);
derivSigmaf = zeros(nFactors,nElem);
for iElem=1:nElem
    derivSigmaf(idx1(iElem),iElem) = b(idx2(iElem),1);
    derivSigmaf(idx2(iElem),iElem) = b(idx1(iElem),1);
end
derivLambda = [V11 zeros(nFactors,nFactors) derivSigmaf];
lambdaCov   = derivLambda*thetaCov*derivLambda';

% Estimating the covariance matrix for lambda
seLambda    = sqrt(diag(lambdaCov));
tstatLambda = lambda./seLambda;

% Computing the r-square
fitVal      = V21*b;
r2          = 1-(mu2-fitVal)'*(mu2-fitVal)/((mu2-mumu2)'*(mu2-mumu2));

%% Estimating and testing the HJ distance
% ########################################################################### %
%{
    We estimate the Hansen-Jagannathan (1997) (HJ) distance for the asset 
    pricing model and simulate p-values following Parker and Julliad (2005). 
%}
% ########################################################################### %

% Computing HJ distance
uMat    = V22+mu2*mu2'; % Second moment matrix of returns
W       = inv(uMat);
hjDist  = sqrt(mu2'*W*mu2 - mu2'*W*V21*inv(V12*W*V21)*V12*W*mu2);

% Computing the A matrix and its N-K non-zero eigenvalues
S       = hacAndrews([u1 u2]);
S22     = S(1:nSeries,1:nSeries);
S1      = chol(S22,'upper');
W1      = chol(W,'upper');
W1      = W1';
A       = S1*W1*(eye(nSeries)-W1'*V21*inv(V12*W*V21)*V12*W1)*W1'*S1';
eigVal  = eig(A);

% Simulating the asymptotic distribution of the HJ distance
nSim    = 5000;
pVal    = zeros(nSim,1);
rChi    = chi2rnd(1,freedomDegrees,nSim);

for iSim = 1:nSim

    % Computing the i'th draw of the asymptotic distribution
    iDistDraw       = sum(eigVal(1:freedomDegrees).*rChi(:,iSim));

    % Determining if we reject or not. 
    pVal(iSim,1)    = ( nObs*hjDist^2 > iDistDraw );

end

% Computing p-value for HJ distance
pval = 1 - mean(pVal);

%% Computing Shanken's CSRT test
% ########################################################################### %
%{
    Following Kan, Robotti, and Shanken (2013), we consider a generalized
    version of Shanken's (1985) CSRT test.  
%}
% ########################################################################### %

% Shanken CSRT test
beta    = V21*inv(V11);
e       = mu2-beta*lambda;
Wh      = sqrtm(eye(nSeries));
WP      = Wh*null(beta'*Wh);
g       = (pfReturns*WP).*(sdf*ones(1,nSeries-nFactors));
S       = hacAndrews(g);
e1      = WP'*e;
csrt    = e1'*inv(S)*e1;

% Asymptotic p-value of CSRT statistic 
pAsym   = 1 - chi2cdf(nObs*csrt,nSeries-nFactors);

% Approximate finite sample p-value of CSRT statistic
pApprox = 1 - fcdf(csrt*(nObs-nSeries+1)/(nSeries-nFactors),nSeries-nFactors,nObs-nSeries+1);

% ########################################################################### %
%% Building the output structure
% ########################################################################### %

gmmResults.b        = b;
gmmResults.seb      = seb;
gmmResults.btstat   = btstat;
gmmResults.lambda   = lambda;
gmmResults.selambda = seLambda;
gmmResults.tlambda  = tstatLambda;
gmmResults.r2       = r2.*100;
gmmResults.hj       = hjDist;
gmmResults.hjp      = pval;
gmmResults.csrt     = csrt;
gmmResults.crstp    = [pAsym pApprox];

end
    
% ########################################################################### %
% [EOF]
% ########################################################################### %
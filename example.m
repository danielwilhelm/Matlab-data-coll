%% Example/Application of covSelection with OGA using Simulated Data
%  Daniel Wilhelm & Nicolas Cerkez
%  University College London, 2017

%% Clear Workspace
clear
clc
close

%% Simulating Data
% In order to be able to use the OGA algorithm, we need (at least) an outcome
% variable (Y), a matrix (X) of regressors, a cost function (c), and a budget
% constraint (B). 

% Simluate data for Y (Nx1 vector: N observations of the outcome)
Y = randperm(1000000, 100)';

% Simluate data for X (Nxp matrix: N observations of M covariates)
X = randperm(1000000, 100*10);
X = reshape(X, 100, 10);

% We assume that the cost function takes the following linear form: 
% c = n*S, where n = # of observations and S = a vector of indices indicating 
% what covariates of X should be included (price is normalized to 1). 
% This could easily be done by writing c = @(S, N) N*sum(S), however, since
% we wrote a function called linear_cfun with exactly this function, we here
% call this functino instead.
c = @linear_cfun;

% Assume that the budget constraint is equal to 200
B = 200;

%% Using OGA to get coefficients of covariates to include
%  Given that we are only looking at the OGA method here, one way to solve
%  this problem is to directly use the OGA function. The OGA application
%  uses the same input arguments as the covSelection function and gives the same
%  output arguments (plus Ihat, which is a vector indexing which covariates
%  to use):
[Ihat, coeff, intercept] = OGA(Y, X, c, B);

%% Using covSelection to get coefficients of covariates to include
%  We now use the covSelection function with four input arguments (Y, X, c, B)
%  and expect two output arguments: 
%  1) coeff = coefficient of regressors to be included (without intercept)
%  2) intercept = intercept of regression 
[coeff2, intercept2] = covSelection(Y, X, c, B);

% The outputs of both, the OGA and the covSelection functions, are, of course, the same!


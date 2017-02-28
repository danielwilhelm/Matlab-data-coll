%% Example/Application of covSelection with OGA using Simulated Data
%  Daniel Wilhelm & Nicolas Cerkez
%  University College London, 2017

%% Clear Workspace
clear
clc
close

%% Simulating Data
% In order to be able to use the OGA algorithm, we need (at least) an outcome
% variable (Y), a matrix (X) of regressors, a cost function (c), a budget
% constraint (B), and a vector of possible sample sizes (n). 

% Simluate data for Y (Nx1 vector: N observations of the outcome)
Y = randperm(1000000, 100)';

% Simluate data for X (Nxp matrix: N observations of M covariates)
X = randperm(1000000, 100*10);
X = reshape(X, 100, 10);

% We assume that the cost function takes the following linear form: 
c = @linear_cfun;

% We can also assume that the cost function takes "the general form" 
d = @general_cfun;
% for this cost function, we need to specify the following: 
%M = 10;         % this cannot exceed length of S?!?!
%tau = [1:M];
%tau0 = 100;
%phi = 1473;
%alpha = 0.4;
%ETA = 430;

% Assume that the budget constraint is equal to 200
B = 200;

% Suppose we have three different sample sizes (they must be in ascending order)
n = [40, 80, 100];

%% Using covSelection to get coefficients of covariates to include
%  We now use the covSelection function with four input arguments (Y, X, c, B, n)
%  and expect three output arguments: 
%  1) coeff = coefficient of regressors to be included (without intercept)
%  2) intercept = intercept of regression 
%  3) ssize = the optimal sample size
tic
[coeff, intercept, ssize] = covSelection(Y, X, c, B, n);
toc

disp('The coefficients of the covariates to be included are:') 
celldisp(coeff)
disp('The intercept of the regression is:')
celldisp(intercept)
disp('The desired sample size is:')
celldisp(ssize)


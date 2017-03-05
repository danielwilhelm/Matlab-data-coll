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

% Assume that the budget constraint is equal to the following
B = 450;

% Suppose we have the following sample sizes
n = [80, 123, 200, 197, 603];

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
disp(coeff)
disp('The intercept of the regression is:')
disp(intercept)
disp('The desired sample size is:')
disp(ssize)

% Notes/Comments:  
% covSelection also prints a table with the different sample sizes and the 
% corresponding MSE, ratio of cost over budget, and the number of covariates 
% to be inlcuded. The table only prints the sample sizes that select covariates. 
% If zero covariates are selected, the sample size, MSE etc. are not included in the table. 
% If the budget is specified in an unrealistic manner, it is possible that
% covariates have a coefficient==Inf. covSelection was written in a way not
% to include these cases (i.e. the table also does not print these). 

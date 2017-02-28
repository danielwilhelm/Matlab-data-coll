%% An Example of a Linear Cost Function
% Daniel Wilhelm & Nicolas Cerkez
% University College London, 2017

function [c] = linear_cfun(S, ss)
    % Create a function with the following input and output arguments: 
    % Inputs:  S = indices indexing which covariates to include
    %          ss = number of observations (sample size)
    % Output:  c = the cost function to be used
    
    % This cost function is of the following linear form: 
    %               c(S, ss) = ss * sum(S) 
    % where ss is the number of observations and S is the
    % a vector of indices indicating what covariates of X should be
    % included (the price of each covariate has been normalized to
    % one). 
    
    c = ss * sum(S);
    
end

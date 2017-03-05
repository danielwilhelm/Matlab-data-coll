%% Peforming the Covariate Selection and Returning the Regression coefficients
%  Daniel Wilhelm & Nicolas Cerkez
%  University College London, 2017

% Note: This is a version of the covSelection function that only uses the
% OGA algorithm

function [coeff, intercept, ssize] = covSelection(Y, X, c, B, n, varargin)
    % Create a function with the following input and output arguments: 
    % Inputs:   Y = outcome vector (Nx1)
    %           X = covariates matrix (Nxp)
    %           c = a cost function 
    %           B = the available budget
    %           n = a vector (1xn) of possible sample sizes (in ascending order!!)
    %           varargin = other (optional) variables (see below)
    % Outputs:  coeff = coefficient of regression (without intercept)
    %           intercept = intercept of regression 
    %           ssize = the optimal sample size
	
    % parse the input arguments
	iP = inputParser;
	iP.addRequired('Y', @isnumeric);                        % Nx1 vector: N observations of the outcome
	iP.addRequired('X', @isnumeric);                        % Nxp matrix: N observations of p covariates
    iP.addRequired('c', @(f) isa(f, 'function_handle'));    % cost function for a given sample size n, i.e. c: {0,1}^M -> R
	iP.addRequired('B', @isnumeric);                        % the available budget
	iP.addRequired('n', @isnumeric);                        % 1xn vector: all possible sample sizes
    iP.addParameter('method', 'OGA', @(x) any(strcmpi(x, {'OGA'})));	% method to be used for covariate selection (only OGA for now)
	iP.addParameter('includeReg', NaN, @isnumeric);        % Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be forced
	iP.addParameter('excludeReg', NaN, @isnumeric);        % Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be prevented

    iP.parse(Y, X, c, B, n, varargin{:});
	Y = iP.Results.Y; X = iP.Results.X; c = iP.Results.c; B = iP.Results.B; n=iP.Results.n; method = iP.Results.method; includeReg = iP.Results.includeReg; excludeReg = iP.Results.excludeReg;
    
    if (~isnan(includeReg) & (length(includeReg)~=size(X,2)))
		error('covSelection: includeReg of incorrect dimension!')
	end;
	if (~isnan(excludeReg) & (length(excludeReg)~=size(X,2)))
		error('covSelection: excludeReg of incorrect dimension!')
	end;
	if (~isnan(includeReg) & ~isnan(excludeReg) & includeReg'*excludeReg~=0)
		error('covSelection: cannot simultaneously include and exclude a covariate!')
	end;
	
	coeff = zeros(size(X,2),1);                 % define coeff = vector of zeros with rows = # of columns of X
    
    MSE =[];                                    % placeholder
    coeff2 =[];                                 % placeholder
    sampsize =[];                               % placeholder
    intercept1 =[];                             % placeholder
    cost =[];                                   % placeholder
    number =[];                                 % placeholder
    
    for i=1:length(n)                           % loop through the number of sample sizes suggested (n)
        ss = n(i);                              % for clarity, label the i'th sample size ss
        switch method
            case {'OGA'}                        % use orthogonal greedy algorithm (OGA) using general budget constraint
                [Ihat, coeffS, intercept] = OGA(Y, X, c, B, ss, 'includeReg', includeReg, 'excludeReg', excludeReg);    % see OGA function for a precise description
                coeff = zeros(size(X, 2), 1);   
                coeff(Ihat) = coeffS;           % define coeff (at indices/covariates to be included) = coefficients from OGA algorithm
        end;
        XS = [ones(size(X,1),1), X];            % include a constant into X matrix
        betahat = [intercept; coeff];           % create a vector betahat that combines intercept and coefficients
        yhat = XS*betahat;                      % calculate predicted outcome values using betahat
        MSE1 = (1/ss)*mean((yhat - Y).^2);      % calculate mean squared error MSE
        MSE(:,i) = MSE1;                        % store the MSE in a row vector 
        coeff2(:,i) = coeff;                    % store the coefficients in the i'th column of this matrix
        sampsize(:,i) = ss;                     % store the sample sizes in the i'th column of this matrix
        intercept1(:,i) = intercept;            % store the intercepts in the i'th column of this matrix
        cost(:,i) = c(sum(Ihat~=0), ss);        % store the associated cost in a matrix
        number(:,i) = sum(Ihat~=0);             % store the associated number of covariates selected in a matrix
    end
    
    % take care of case where zero covariates are selected
    columnsWithAllZeros = all(coeff2 == 0);             % find which column of coeff2 has no covariates selected
    coeff2 = coeff2(:, ~columnsWithAllZeros);           % create a new coeff2 matrix without the column with all zeros (without the sample sizes that pick zero covariates)
    intercept1 = intercept1(:, ~columnsWithAllZeros);   % adjust intercept matrix accordingly
    sampsize = sampsize(:, ~columnsWithAllZeros);       % adjust sample size matrix accordingly
    cost = cost(:, ~columnsWithAllZeros);               % adjust cost matrix accordingly
    number = number(:, ~columnsWithAllZeros);           % adjust number of covariates matrix accordingly
    MSE = MSE(:, ~columnsWithAllZeros);                 % adjust MSE vector accordingly
    n = n(:, ~columnsWithAllZeros);                     % adjust sample size vector accordingly
    
    % take care of case where covariates with coefficients==Inf are selected
    columnsWithAllZeros1 = all(coeff2 == Inf);           % find which column of coeff2 has coefficients==Inf
    coeff2 = coeff2(:, ~columnsWithAllZeros1);           % create a new coeff2 matrix without the column with all zeros (without the sample sizes that pick zero covariates)
    intercept1 = intercept1(:, ~columnsWithAllZeros1);   % adjust intercept matrix accordingly
    sampsize = sampsize(:, ~columnsWithAllZeros1);       % adjust sample size matrix accordingly
    cost = cost(:, ~columnsWithAllZeros1);               % adjust cost matrix accordingly
    number = number(:, ~columnsWithAllZeros1);           % adjust number of covariates matrix accordingly
    MSE = MSE(:, ~columnsWithAllZeros1);                 % adjust MSE vector accordingly
    n = n(:, ~columnsWithAllZeros1);                     % adjust sample size vector accordingly

    [~, ind] = min(MSE);                                % find the min MSE
    coeff=coeff2(:,ind);                                % define coeff as the coeff that are the max number of covariates for the given sample sizes
    intercept = intercept1(:,ind);                      % define the corresponding intercept
    ssize = sampsize(:,ind);                            % define the corresponding sample size
    
    % plot the various sample sizes, the corresponding MSE, cost, and number of covariates selected in a table
    table(n', MSE', cost'/B, number','VariableNames', {'SampleSize'; 'MSE'; 'CostoverBudget'; 'NumberCovariates'})
end
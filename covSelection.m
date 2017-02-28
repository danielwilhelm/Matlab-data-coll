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
    
    coeff1 = {};                                % placeholder for cell array below
    coeff3 = {};                                % placeholder for cell array below
    intercept1 = {};                            % placeholder for cell array below
    sampsize = {};                              % placeholder for cell array below
    
    for i=1:length(n)                           % loop through the number of sample sizes suggested (n)
        ss = n(i);                              % for clarity, label the i'th sample size ss
        switch method
            case {'OGA'}                        % use orthogonal greedy algorithm (OGA) using general budget constraint
                [Ihat, coeffS, intercept] = OGA(Y, X, c, B, ss, 'includeReg', includeReg, 'excludeReg', excludeReg);    % see OGA function for a precise description
                coeff(Ihat) = coeffS;           % define coeff (at indices/covariates to be included) = coefficients from OGA algorithm
        end;
        coeff2(:, i) = coeff;                   % store the resulting coefficients in the i'th cell of this matrix
        coeff1{i} = coeff;                      % store the resulting coefficients in the i'th cell of this cell array
        intercept1{i} = intercept;              % store the resulting intercepts in the i'th cell of this cell array
        sampsize{i} = ss;                       % store the i'th sample size (ss) in the i'th column of this matrix
    end
    for i = 1:size(coeff2, 2)                   % for each column of the coeff2 matrix
        [~, ~, v] = find(coeff2(:, i));         % create a vector v with the non-zero elements of coeff2
        coeff3{i} = {v};                        % store the resulting coefficients in the i'th cell of this cell array
    end
    x = cellfun('length',coeff3);               % create a (1xp) vector that counts the number of elements in the coeff1 array (so count which n gave how many coefficients)
    % the row vector x contains a bunch of numbers each counting how many
    % covariates each sample size chooses. The goal is to pick the largest,
    % however, we must take into account that two sample sizes may pick the
    % same number of covariates. In that case, we pick the larger of the
    % two sample sizes (this is why the possible sample sizes have to be 
    % in ascending order). 
    maxval = max(x);                            % find the max value (number of coefficients) of x 
    idx = find(x == maxval);                    % find all indices in x that correspond to this max value
    ind = max(idx);                             % define ind as the index of x with the largest sample size and the most coefficients (covariates)
    coeff=coeff1(ind);                          % define coeff as the coeff that are the max number of covariates for the given sample sizes
    intercept = intercept1(ind);                % define the corresponding intercept
    ssize = sampsize(ind);                      % define the corresponding sample size
end
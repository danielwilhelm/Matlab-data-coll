%% Peforming the Covariate Selection and Returning the Regression coefficients
%  Daniel Wilhelm & Nicolas Cerkez
%  University College London, 2017

% Note: This is a version of the covSelection function that only uses the
% OGA algorithm

function [coeff, intercept] = covSelection(Y, X, c, B, varargin)
    % Create a function with the following input and output arguments: 
    % Inputs:   Y = outcome vector (Nx1)
    %           X = covariates matrix (Nxp)
    %           c = a cost function 
    %           B = the available budget
    %           varargin = other (optional) variables (see below)
    % Outputs:  coeff = coefficient of regression (without intercept)
    %           intercept = intercept of regression 
	
    % parse the input arguments
	iP = inputParser;
	iP.addRequired('Y', @isnumeric);                        % Nx1 vector: N observations of the outcome
	iP.addRequired('X', @isnumeric);                        % Nxp matrix: N observations of p covariates
    iP.addRequired('c', @(f) isa(f, 'function_handle'));    % cost function for a given sample size n, i.e. c: {0,1}^M -> R
	iP.addRequired('B', @isnumeric);                        % the available budget
	iP.addParameter('method', 'OGA', @(x) any(strcmpi(x, {'OGA'})));	% method to be used for covariate selection (we only have OGA for now)
	iP.addParameter('includeReg', NaN, @isnumeric);        % Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be forced
	iP.addParameter('excludeReg', NaN, @isnumeric);        % Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be prevented

    iP.parse(Y, X, c, B, varargin{:});
	Y = iP.Results.Y; X = iP.Results.X; c = iP.Results.c; B = iP.Results.B; method = iP.Results.method; includeReg = iP.Results.includeReg; excludeReg = iP.Results.excludeReg;
    
    if (~isnan(includeReg) & (length(includeReg)~=size(X,2)))
		error('covSelection: includeReg of incorrect dimension!')
	end;
	if (~isnan(excludeReg) & (length(excludeReg)~=size(X,2)))
		error('covSelection: excludeReg of incorrect dimension!')
	end;
	if (~isnan(includeReg) & ~isnan(excludeReg) & includeReg'*excludeReg~=0)
		error('covSelection: cannot simultaneously include and exclude a covariate!')
	end;
	
	coeff = zeros(size(X,2),1);         % define coeff = vector of zeros with rows = # of columns of X
	
	switch method
		
		% use orthogonal greedy algorithm (OGA) using general budget constraint
		case {'OGA'}
			[Ihat, coeffS, intercept] = OGA(Y, X, c, B, 'includeReg', includeReg, 'excludeReg', excludeReg);    % see OGA function for a precise description
			coeff(Ihat) = coeffS;       % define coeff (at indices/covariates to be included) = coefficients from OGA algorithm
	end;
end
%% Orthogonal Greedy Algorithm (OGA) with a General Budged Constraint
%  Daniel Wilhelm & Nicolas Cerkez
%  University College London, 2017

function [Ihat, coeff, intercept] = OGA(Y, X, c, B, varargin)
    % Create a function with the following input and output arguments: 
    % Inputs:   Y = outcome vector (Nx1)
    %           X = covariates matrix (Nxp)
    %           c = a cost function 
    %           B = the available budget
    %           varargin = other (optional) variables (see below)
    % Outputs:  Ihat = indicates indices (of covariates) to be included in reg
    %           coeff = coefficient of regression (without intercept)
    %           intercept = intercept of regression 
	
    % parse input arguments
	iP = inputParser;
	iP.addRequired('Y', @isnumeric);                    % Nx1 vector: N observations of the outcome
	iP.addRequired('X', @isnumeric);                    % Nxp matrix: N observations of M covariates
	iP.addRequired('c', @(f) isa(f, 'function_handle'));% cost function for a given sample size n, i.e. c: {0,1}^M -> R
	iP.addRequired('B', @isnumeric);                    % the available budget
	iP.addParameter('maxIter', 1000, @isnumeric);       % max number of regressors to be selected
	iP.addParameter('includeReg', NaN, @isnumeric);     % Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be forced
	iP.addParameter('excludeReg', NaN, @isnumeric);     % Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be prevented

	iP.parse(Y, X, c, B, varargin{:});
	Y = iP.Results.Y; X = iP.Results.X; c = iP.Results.c; B = iP.Results.B; maxIter = iP.Results.maxIter; includeReg = iP.Results.includeReg; excludeReg = iP.Results.excludeReg;

	if (~isnan(includeReg) & (length(includeReg)~=size(X,2)))
		error('OGA: includeReg of incorrect dimension!')
	end;
	if (~isnan(excludeReg) & (length(excludeReg)~=size(X,2)))
		error('OGA: excludeReg of incorrect dimension!')
	end;
	if (~isnan(includeReg) & ~isnan(excludeReg) & includeReg'*excludeReg~=0)
		error('OGA: cannot simultaneously include and exclude a covariate!')
	end;

	% initializations
	M = size(X,2); N = size(X,1);       % M = number of columns of X; N = number of rows of X
	if (maxIter>M)                      
		maxIter=M;                      % make sure the max number of regressors to be selected is not greater than the number of columns of X
	end;
	if (isnan(includeReg))              % if no inclusion of covariates should be forced
		S = zeros(M,1);                 % define S as a (Mx1) vector of zeros
		r = Y;                          % define residuals r = Y 
	else
		S = includeReg;                 % define S = includeReg as a (Mx1) vector
		XS = [ones(N,1) X(:,find(S))];  % define XS as a matrix: first column = ones (for intercept); other columns are the columns of X where S==1
        coeffi = pinv(XS'*XS)*XS'*Y;    % coeffi = OLS estimator of reg Y on XS
        r = Y-XS*coeffi;                % compute residuals
	end;
	
    Ihat = 1:M;                         % define Ihat as a (1xM) vector, where M = # of columns of X
    coeff = Ihat*Inf;                   % define coeff
    intercept = 0;                      % define intercept

    % check whether budget constraint is satisfied when no covariates are selected
    if (c(S, N)<=B)                     % evaluate cost function with all covariates indicated by S (in this case either none or with covariates where inlcudeReg==1), if less than B
    	coeff = zeros(M,1);             % define coeff as a (Mx1) vector
    	intercept = mean(Y);            % define intercept = mean of Y (outcome) (since no covariates included in this case)
    else
    	return                          
    end;
	
    % add covariates as long as budget is satisfied and number of max iterations is not exceeded
	k = 0;                              % define k == 0
    while(k<=maxIter)                   % as long as k is smaller than the max number of regressors to be selected

		% compute correlations between residual r and remaining covariates
		if (isnan(excludeReg))                  % if excludeReg==NaN
			notS = find((1-S));                 % define notS = vector of indices if S==0
        else                                    % else
			notS = find((1-S).*(1-excludeReg)); % define notS = vector of indices when S==0 and excludeReg==0
		end;
		cortmp = corrcoef([r X(:,notS)]);       % compute correlation between residual r and all columns (variables) of X if they are indicated by notS
        cor = cortmp(2:end,1);                  % define cor (a vector) = rows 2:end & first column of cortmp

		% find covariate with maximum correlation
		[m I] = max(cor);               % find maximum correlation of cor and index of the max
		I = notS(I);                    % define I as the covariate with the highest correlation

		% check budget constraint
		Snew = S;                           % define Snew = S
		Snew(I) = 1;                        % indicate the position of this new covariate to be included
		if (c(Snew, N)<=B)                  % evaluate cost function with all covariates indicated by Snew, if less than B
			S = Snew;                       % define S = Snew
			Ihat = find(S);                 % define Ihat = indices if S==1
            k = k + 1;                      % add one iteration/count to k
            XS = [ones(N,1) X(:,Ihat)];     % define XS as above
            coeffi = pinv(XS'*XS)*XS'*Y;    % coeffi = OLS estimator of reg Y on XS
            r = Y-XS*coeffi;                % compute residuals
            coeff = coeffi(2:(size(XS,2))); % define coeff = coefficients of included covariates
            intercept = coeffi(1);          % define intercept = intercept of regression
        else
        	break
		end;
	end;
end
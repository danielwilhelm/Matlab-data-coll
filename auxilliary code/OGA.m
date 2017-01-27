% OGA algorithm with general budget constraint c(S,n) <= B
function [Ihat coeff intercept] = OGA(Y, X, c, B, varargin)

	% parse the input arguments
	iP = inputParser;
	iP.addRequired('Y', @isnumeric);		% Nx1 vector: N observations of the outcome
	iP.addRequired('X', @isnumeric);		% Nxp matrix: N observations of M covariates
	iP.addRequired('c', @(f) isa(f, 'function_handle'));	% cost function for a given sample size n, i.e. c: {0,1}^M -> R
	iP.addRequired('B', @isnumeric);	% budget
	iP.addParamValue('maxIter', 1000, @isnumeric);	% max number of regressors to be selected
	iP.addParamValue('includeReg', NaN, @isnumeric);	% Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be forced
	iP.addParamValue('excludeReg', NaN, @isnumeric);	% Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be prevented

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
	M = size(X,2); N = size(X,1);
	if (maxIter>M) 
		maxIter=M;
	end;
	if (isnan(includeReg))
		S = zeros(M,1); 
		r = Y;
	else
		S = includeReg;
		XS = [ones(N,1) X(:,find(S))];
        coeffi = pinv(XS'*XS)*XS'*Y;
        r = Y-XS*coeffi;
	end;
	k = 0;
	
    Ihat = 1:M; coeff = Ihat*Inf;
    intercept = 0;

    % check whether budget constraint is satisfied when no covariates are selected
    if (c(S)<=B)
    	coeff = zeros(M,1);
    	intercept = mean(Y);
    else
    	return
    end;

	% add covariates as long as budget is satisfied and number of max iterations is not exceeded
	while(k<=maxIter)

		% compute correlations between residual r and remaining covariates
		if (isnan(excludeReg))
			notS = find((1-S));
		else
			notS = find((1-S).*(1-excludeReg));
		end;
		cortmp = corrcoef([r X(:,notS)]); cor = cortmp(2:end,1);

		% find covariate with maximum correlation
		[m I] = max(cor);
		I = notS(I);

		% check budget constraint
		Snew = S;
		Snew(I) = 1;
		if (c(Snew)<=B)
			S = Snew;
			Ihat = find(S);
            k = k + 1;
            XS = [ones(N,1) X(:,Ihat)];
            coeffi = pinv(XS'*XS)*XS'*Y;
            r = Y-XS*coeffi;
            coeff = coeffi(2:(size(XS,2)));
            intercept = coeffi(1);
        else
        	break
		end;
	end;
end
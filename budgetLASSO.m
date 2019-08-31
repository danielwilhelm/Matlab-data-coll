function [coeff fitInfo lambda0] = budgetLASSO(Y, X, c, B, varargin)
	
	% parse the input arguments
	iP = inputParser;
	iP.addRequired('Y', @isnumeric);		% Nx1 vector: N observations of the outcome
	iP.addRequired('X', @isnumeric);		% Nxp matrix: N observations of M covariates
	iP.addRequired('c', @(f) isa(f, 'function_handle'));	% cost function for a given sample size n, i.e. c: {0,1}^M -> R
	iP.addRequired('B', @isnumeric);	% budget
	iP.addParamValue('maxIter', 20, @isnumeric);	% max number of iterations to find optimal lambda
	iP.addParamValue('lbounds', [0 10], @isnumeric);	% lower and upper bounds on lambda
	iP.addParamValue('includeReg', NaN, @isnumeric);	% Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be forced
	iP.addParamValue('excludeReg', NaN, @isnumeric);	% Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be prevented

	iP.parse(Y, X, c, B, varargin{:});
	Y = iP.Results.Y; X = iP.Results.X; c = iP.Results.c; B = iP.Results.B; maxIter = iP.Results.maxIter; lbounds = iP.Results.lbounds; includeReg = iP.Results.includeReg; excludeReg = iP.Results.excludeReg;

	if (~isnan(includeReg) & (length(includeReg)~=size(X,2)))
		error('budgetLASSO: includeReg of incorrect dimension!')
	end;
	if (~isnan(excludeReg) & (length(excludeReg)~=size(X,2)))
		error('budgetLASSO: excludeReg of incorrect dimension!')
	end;
	if (~isnan(includeReg(1)) & ~isnan(excludeReg(1)) & includeReg'*excludeReg~=0)
		error('budgetLASSO: cannot simultaneously include and exclude a covariate!')
	end;

	% Step 1: partial out forced regressors in includeReg and exclude covariates in excludeReg
	if (isnan(includeReg))
		Yres = Y;
		if (isnan(excludeReg))
			Xres = X;
		else
			Xres = X(:,find(1-excludeReg));
		end;
	else
		Xincl = X(:,find(includeReg));
		bhat = inv(Xincl'*Xincl)*Xincl'*Y;
		Yres = Y - Xincl*bhat;
		if (isnan(excludeReg))
			Xres = X(:,find(1-includeReg));
		else
			Xres = X(:,find((1-includeReg).*(1-excludeReg)));
		end;
	end;

	function bg = budgetgap(lambda)
		[coefflasso fitInfo] = lasso(Xres, Yres, 'Lambda', lambda);

		coeff = zeros(size(X,2),1);
		if (isnan(includeReg))
			if (isnan(excludeReg))
				coeff = coefflasso;
			else
				coeff(find(1-excludeReg)) = coefflasso;
			end;
		else
			coeff(find(includeReg)) = bhat;
			if (isnan(excludeReg))
				coeff(find(1-includeReg)) = coefflasso;
			else
				coeff(find((1-includeReg).*(1-excludeReg))) = coefflasso;
			end;
		end;
		S = coeff*0;

		if ~isempty(find(coeff))
			S(find(coeff)) = 1;
		end;
		bg = c(S) - B;
	end;

	% Step 2: perform lasso on residual from above regression
	maxc = c(ones(size(X,2),1));	% max cost when all variables are selected
	if maxc<B
		lambda0 = 0.0000001;
	else
		% run bisection algorithm to find the lambda value that minimizes the budget gap
		if (budgetgap(lbounds(2))<0)
			if (budgetgap(lbounds(1))>0)
				lambda0 = bisect(lbounds(1), lbounds(2), @budgetgap, maxIter);
			else
				lambda0 = 0.0000001;
			end;
		else
			coeff = Inf*ones(size(X,2),1);
			lambda0 = -1;
			fitInfo.Intercept = Inf;
		end;
	end;

	% run LASSO with optimal lambda
	if lambda0>=0
		[coefflasso fitInfo] = lasso(Xres, Yres, 'Lambda', lambda0);
		
		% combine regression coefficients from Steps 1 and 2
		coeff = zeros(size(X,2),1);
		if (isnan(includeReg))
			if (isnan(excludeReg))
				coeff = coefflasso;
			else
				coeff(find(1-excludeReg)) = coefflasso;
			end;
		else
			coeff(find(includeReg)) = bhat;
			if (isnan(excludeReg))
				coeff(find(1-includeReg)) = coefflasso;
			else
				coeff(find((1-includeReg).*(1-excludeReg))) = coefflasso;
			end;
		end;
		Shat = coeff*0;

		if ~isempty(find(coeff))
			Shat(find(coeff)) = 1;
		end;

		if c(Shat) > B
			disp('budget LASSO: selection exceed budget!')
		end;
	else
		disp('budget LASSO: bisection found negative lambda!')
	end;
	
end
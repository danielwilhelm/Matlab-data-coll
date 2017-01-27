% perform covariate selection and return the regression coefficients
function [coeff, intercept] = covSelection(Y, X, varargin)

	% parse the input arguments
	iP = inputParser;
	iP.addRequired('Y', @isnumeric);		% Nx1 vector: N observations of the outcome
	iP.addRequired('X', @isnumeric);		% Nxp matrix: N observations of p covariates
	iP.addParamValue('method', 'OGA', @(x) any(strcmpi(x, {'OGA','POST-OGA','wmpalg','largest','POST-LASSO','LASSO', 'random'})));	% method to be used for covariate selection
	iP.addParamValue('m', 10, @isnumeric);	% number of regressors to be selected
	iP.addParamValue('c', @(f) isa(f, 'function_handle'));	% cost function for a given sample size n, i.e. c: {0,1}^M -> R
	iP.addParamValue('B', @isnumeric);	% budget
	iP.addParamValue('lambda', -1, @isnumeric);	% penalization parameter for lasso (negative value: do cross-validation to select lambda)
	iP.addParamValue('trueCoeff', @isnumeric);	% number of regressors to be selected
	iP.addParamValue('post', 0, @isnumeric);	% do post-lasso?
	iP.addParamValue('cv', 0, @isnumeric);	% LASSO: do cross-validation to select lambda? otherwise choose lambda so as to satisfy the budget constraint
	iP.addParamValue('perm', -1, @isnumeric);	% random choice: permutations to use; if <0, then draw randomly here
	iP.addParamValue('includeReg', NaN, @isnumeric);	% Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be forced
	iP.addParamValue('excludeReg', NaN, @isnumeric);	% Mx1 vector of zeros and ones, indicating (=1) regressors whose inclusion should be prevented

	iP.parse(Y, X, varargin{:});
	Y = iP.Results.Y; X = iP.Results.X; m = iP.Results.m; lambda = iP.Results.lambda; trueCoeff = iP.Results.trueCoeff; method = iP.Results.method; c = iP.Results.c; B = iP.Results.B; cv = iP.Results.cv; perm = iP.Results.perm; includeReg = iP.Results.includeReg; excludeReg = iP.Results.excludeReg;
	
	if (~isnan(includeReg) & (length(includeReg)~=size(X,2)))
		error('covSelection: includeReg of incorrect dimension!')
	end;
	if (~isnan(excludeReg) & (length(excludeReg)~=size(X,2)))
		error('covSelection: excludeReg of incorrect dimension!')
	end;
	if (~isnan(includeReg) & ~isnan(excludeReg) & includeReg'*excludeReg~=0)
		error('covSelection: cannot simultaneously include and exclude a covariate!')
	end;
	post = 0;
	if sum(strcmp(method, { 'POST-LASSO', 'POST-OGA' }))>0
		post = 1;
	end;
	coeff = zeros(size(X,2),1);
	
	switch method
		
		% use own orthogonal greedy algorithm using general budget constraint
		case { 'OGA', 'POST-OGA' }
			[Ihat, coeffS, intercept] = OGA(Y, X, c, B, 'includeReg', includeReg, 'excludeReg', excludeReg);
			coeff(Ihat) = coeffS;

		% use the orthogonal greedy algorithm implemented in Matlab
		% case 'wmpalg'
		% 	[f, r, g, Ihat] = wmpalg('OMP', Y, X, 'itermax', m);
		% 	coeff(Ihat) = pinv(X(:,Ihat)'*X(:,Ihat))*X(:,Ihat)'*Y;
		
		% return the m largest population coefficients
		case 'largest'
			if (~isnan(includeReg))
				error('includeReg is not NaN: forcing the inclusion of certain covariates is not available for largest covariate selection algorithm!')
			end;
			if (~isnan(excludeReg))
				error('excludeReg is not NaN: preventing the inclusion of certain covariates is not available for largest covariate selection algorithm!')
			end;
			[Bhat Im] = sort(trueCoeff,'descend');
			coeff(Im(1:m)) = Bhat(1:m);
			intercept = 0;

		% use the LASSO algorithm implemented in Matlab
		case { 'LASSO', 'POST-LASSO' }

			if lambda < 0	% find optimal lambda
				
				if cv 	% do cross-validation
					% [Bhat fitInfo] = lasso(X, Y, 'CV', length(Y));
					% coeff = Bhat(:,fitInfo.IndexMinMSE);
					[Bhat fitInfo] = lasso(X, Y, 'NumLambda', 10);
					[m ind] = min(fitInfo.MSE);
					coeff = Bhat(:,ind);

				else 	% find smallest lambda that satisfies the budget constraint
					[coeff fitInfo lambda0] = budgetLASSO(Y, X, c, B, 'includeReg', includeReg, 'excludeReg', excludeReg);
				end;
			else 
				[coeff fitInfo] = lasso(X, Y, 'Lambda', lambda);
% 				if lambda0<0
% 					coeff = zeros(size(X,2),1);
% 				end;
			end;
			intercept = fitInfo.Intercept;

		% randomly select m covariates
		case 'random'
			if (~isnan(includeReg))
				error('includeReg is not NaN: forcing the inclusion of certain covariates is not available for random covariate selection algorithm!')
			end;
			if (~isnan(excludeReg))
				error('excludeReg is not NaN: preventing the inclusion of certain covariates is not available for random covariate selection algorithm!')
			end;
			if perm<0
				Ihat = randsample(size(X,2), m);
			else
				Ihat = perm(1:m);
			end;
			XS = [ ones(size(X,1),1) X(:,Ihat) ];
			coeffi = pinv(XS'*XS)*XS'*Y;
			coeff(Ihat) = coeffi(2:size(X,2));
			intercept = coeffi(1);
	end;

	% in the case of POST-LASSO or POST-OGA, fit the regression one more time
	if (post)
		if (sum(abs(coeff)>0)>0) & (sum(coeff==Inf)==0)
			Ihat = find(coeff);
			XS = [ ones(size(X,1),1) X(:,Ihat) ];
			coeffi = pinv(XS'*XS)*XS'*Y;
			coeff(Ihat) = coeffi(2:size(XS,2));
			intercept = coeffi(1);
		else
			intercept = mean(Y);
		end;
	end;
end
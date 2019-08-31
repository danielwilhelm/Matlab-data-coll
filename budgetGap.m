function gap = budgetGap(MSEtarget, Ytrain, Xtrain, Yeval, Xeval, scriptN, price, varargin)

	iP = inputParser;
	iP.addRequired('MSEtarget', @isnumeric);		% scalar: target MSE
	iP.addRequired('Ytrain', @isnumeric);		% observations of outcome in training sample
	iP.addRequired('Xtrain', @isnumeric);		% observations of p covariates in training sample
	iP.addRequired('Yeval', @isnumeric);		% observations of outcome in evaluation sample
	iP.addRequired('Xeval', @isnumeric);		% observations of p covariates in evaluation sample
	iP.addRequired('scriptN', @isnumeric);		% nNx1 vector: grid of sample sizes
	iP.addRequired('price', @isnumeric);		% Mx1 vector: price vector
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

	iP.parse(MSEtarget, Ytrain, Xtrain, Yeval, Xeval, scriptN, price, varargin{:});
	price = iP.Results.price; scriptN = iP.Results.scriptN; MSEtarget = iP.Results.MSEtarget; Ytrain = iP.Results.Ytrain; Xtrain = iP.Results.Xtrain; Yeval = iP.Results.Yeval; Xeval = iP.Results.Xeval; m = iP.Results.m; lambda = iP.Results.lambda; trueCoeff = iP.Results.trueCoeff; method = iP.Results.method; c = iP.Results.c; B = iP.Results.B; cv = iP.Results.cv; perm = iP.Results.perm; includeReg = iP.Results.includeReg; excludeReg = iP.Results.excludeReg;
	post = 0;
	if sum(strcmp(method, { 'POST-LASSO', 'POST-OGA' }))>0
		post = 1;
	end;

	Ntrain = size(Xtrain,1); Neval = size(Xeval,1); M = size(Xtrain,2); nN = length(scriptN); intercept = zeros(nN,1); gammahat = zeros(nN, M);

	% loop through sample sizes in scriptN
    for n=1:nN

    	% define cost function given sample size
		cofS = @(S) c(S,scriptN(n),price);

        % loop through different covariate selection methods
        if cofS(zeros(M,1))<= B
	        	        	
        	% estimate regression coefficients by covariate selection method
        	switch method
    			case 'LASSO'
		        	[gammahat(n,:) intercept(n)] = covSelection(Ytrain, Xtrain, 'method', method, 'lambda', lambda, 'c', cofS, 'B', B, 'trueCoeff', trueCoeff, 'includeReg', includeReg, 'excludeReg', excludeReg);
		        case 'POST-LASSO'
		        	[gammahat(n,:) intercept(n)] = covSelection(Ytrain, Xtrain, 'method', method, 'lambda', lambda, 'c', cofS, 'B', B, 'trueCoeff', trueCoeff, 'includeReg', includeReg, 'excludeReg', excludeReg);
        		otherwise
        			[gammahat(n,:) intercept(n)] = covSelection(Ytrain, Xtrain, 'method', method, 'm', m(n), 'c', cofS, 'B', B, 'trueCoeff', trueCoeff, 'includeReg', includeReg, 'excludeReg', excludeReg);
    		end;

        	% compute predicted values
        	XSeval = [ ones(Neval,1) Xeval ];
        	fhatN = XSeval*[intercept(n); gammahat(n,:)'];

			% compute prediction error
    		RN(n) = mean((Yeval-fhatN).^2);

        	% compute MSE (i.e. prediction error divided by n)
        	MSEN(n) = RN(n) / scriptN(n);
	    else
	    	MSEN(n) = Inf; 
	    end;
    end;

    % pick the sample size yielding the smallest MSE 
    gap = MSEtarget - min(MSEN);


end
function [Y, X, D] = generateData(N, M, model, TREATMENT, varargin)

	% parse the input arguments
	iP = inputParser;
	iP.addRequired('N', @isnumeric);		
	iP.addRequired('M', @isnumeric);		
	iP.addRequired('model', @isstruct);
	iP.addRequired('TREATMENT', @isnumeric);
	iP.addParamValue('covariateData', @isnumeric);
	iP.addParamValue('outcomeData', @isnumeric);

	iP.parse(N, M, model, TREATMENT, varargin{:});
	N = iP.Results.N; M = iP.Results.M; model = iP.Results.model; TREATMENT = iP.Results.TREATMENT; covariateData = iP.Results.covariateData; outcomeData = iP.Results.outcomeData;


	% generate covariates
	if strcmp(model.correlations, 'empirical')	% generate covariates from data
		NEmp = size(covariateData, 1);

		% draw individuals
		rind = randsample(NEmp, N, true);

		% collect the covariates for each individual
		X = covariateData(rind,:);

		% estimate coefficients in data
		coeff0hat = (covariateData'*covariateData)^(-1)*covariateData'*outcomeData;

		% estimate residual standard deviation in the data
		epsilon0hat = outcomeData - covariateData*coeff0hat;
		sigmaEpsilon = sqrt(var(epsilon0hat));

		% generate the regression coefficient sequence
		coeff = coeff0hat + sign(coeff0hat).*getCoefficients(M, model.coefficients, model.scale);

	else 
		% construnct the covariance matrix of the covariates
		switch model.correlations
			
			case 'orthogonal'		% covariates that are all mutually uncorrelated
				Sigma = eye(M);
				
			case 'constant'			% covariates with equal correlation between all of them
				correlation = 0.6;
				Sigma = correlation*ones(M) + (1 - correlation)*eye(M);

			case 'decaying'			% covariates with correlation decrasing with distance between their indices	
				rho = 0.6;
				A = ones(M,1)*linspace(1,M,M);
				Sigma = (ones(M)*rho) .^(abs(A-A'));			
		end

		% generate covariates X
		mu = zeros(M,1);
		X = mvnrnd(mu, Sigma, N);	

		% set residual standard deviation
		sigmaEpsilon = 1;

		% generate the regression coefficient sequence
		coeff = getCoefficients(M, model.coefficients, model.scale);
	end;

	

	% generate random noise, treatment variable and outcome
	epsilon = sigmaEpsilon*randn(N,1);

	if TREATMENT
		D = randi([0 1], N, 1);
	else
		D = zeros(N, 1);
	end;

	Y = model.beta*D + X*coeff + epsilon;
end
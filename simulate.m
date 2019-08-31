function [gammahat, nSelCov, fhat, nhatInd, R, Roracle, MSE, betahat, sebetahat, m, err, methodLabelsLambda, cost, eqb] = simulate(nMC, N, M, model, methods, lambda, scriptN, price, c, B, varargin)

	% parse the input arguments
	iP = inputParser;
	iP.addRequired('nMC', @isnumeric);
	iP.addRequired('N', @isnumeric);
	iP.addRequired('M', @isnumeric);
	iP.addRequired('model', @isstruct);
	iP.addRequired('methods');
	iP.addRequired('lambda', @isnumeric);
	iP.addRequired('scriptN', @isnumeric);
	iP.addRequired('price', @isnumeric);
	iP.addRequired('c', @(f) isa(f, 'function_handle'));
	iP.addRequired('B', @isnumeric);
	iP.addParamValue('nodisplay', 1, @isnumeric);
	iP.addParamValue('covariateData', @isnumeric);
	iP.addParamValue('outcomeData', @isnumeric);

	iP.parse(nMC, N, M, model, methods, lambda, scriptN, price, c, B, varargin{:});
	nMC = iP.Results.nMC; N = iP.Results.N; M = iP.Results.M; model = iP.Results.model; methods = iP.Results.methods; lambda = iP.Results.lambda; scriptN = iP.Results.scriptN; price = iP.Results.price; c = iP.Results.c; B = iP.Results.B; nodisplay = iP.Results.nodisplay; covariateData = iP.Results.covariateData; outcomeData = iP.Results.outcomeData;

	% get population regression coefficients
	gamma = getCoefficients(M, model.coefficients, model.scale);

	% construct labels for methods
	if (strcmp(methods{1}, 'experiment'))
		methodIndexStart = 2;
	else
		methodIndexStart = 1;
	end;
	methodLabels = { };
	for i=1:length(methods)
		if strcmp(methods{i}, 'LASSO') & (length(lambda)>1)
			for j=1:length(lambda)
				methodLabels = [ methodLabels, 'LASSO' ];
			end;
		else
			methodLabels = [ methodLabels, methods{i} ];
		end;
	end;

	% other initializations
	nB = length(B); nN = length(scriptN); nMethods = length(methodLabels);
	fhatN=zeros(nMC,N,nMethods,nN); gammahatN=zeros(nMC,M,nMethods,nN); interceptN=zeros(nMC,nMethods,nN); intercept=zeros(nMC,nMethods); err=zeros(nMC,1); nhatInd=zeros(nMC,nMethods); nSelCov=nhatInd;
	fhat=zeros(nMC,N,nMethods); gammahat=zeros(nMC,M,nMethods); R=zeros(nMC,nMethods,nB); Roracle=zeros(nMC,nB); RN=zeros(nMC,nMethods,nN,nB); MSEN=RN; MSE=R;
	betahat=zeros(nMC,nMethods); sebetahat=betahat; cost=zeros(nMC,nMethods,nB); eqb = zeros(nMC, nMethods);

	% ------------------------------------------------------------------------
	% perform Monte Carlo simulations
	% ------------------------------------------------------------------------

	% loop through all budgets
	if (~nodisplay)
		h = waitbar(0,'Please wait ...'); step=0;
		steps = nB*nMC*nN*nMethods;
	end;
	tic;
	for b=1:length(B)

		% compute the number of covariates to be collected
		m = ones(1,length(scriptN));

		% loop through MC samples
		for mc=1:nMC
		        
		        % generate data
		        [Y, X, D] = generateData(N, M, model, 0, 'outcomeData', outcomeData, 'covariateData', covariateData, 'outcomeData', outcomeData);

		        % estimate model as in experiment
		        if (strcmp(methods{1}, 'experiment'))
			        XS = [ ones(size(X,1),1) X ];
					coeffi = pinv(XS'*XS)*XS'*Y;
					gammahat(mc,:,1) = coeffi(2:end);
					fhat(mc,:,1) = XS*coeffi;
					R(mc,1,b) = mean((Y-fhat(mc,:,1)').^2);
					MSE(mc,1,b) = R(mc,1,b) / size(X,1);
					cost(mc,1,b) = c(ones(M,1), N, price);
				end;

		        % loop through sample sizes in scriptN
		        for n=1:nN

		        	% define cost function given sample size
					cofS = @(S) c(S,scriptN(n),price);

			        % loop through different covariate selection methods
			        if cofS(zeros(M,1))<= B
				        
				        LASSOind=1; PLASSOind=1;
				        for i=methodIndexStart:nMethods
				        	
				        	% estimate regression coefficients by covariate selection method
				        	switch methodLabels{i}
			        			case 'LASSO'
						        	[gammahatN(mc,:,i,n) interceptN(mc,i,n)] = covSelection(Y, X, 'method', methodLabels{i}, 'lambda', lambda(LASSOind), 'c', cofS, 'B', B(b), 'trueCoeff', gamma);
						        	LASSOind = LASSOind + 1;
						        case 'POST-LASSO'
						        	[gammahatN(mc,:,i,n) interceptN(mc,i,n)] = covSelection(Y, X, 'method', methodLabels{i}, 'lambda', lambda(PLASSOind), 'c', cofS, 'B', B(b), 'trueCoeff', gamma);
						        	PLASSOind = PLASSOind + 1;
				        		otherwise
				        			[gammahatN(mc,:,i,n) interceptN(mc,i,n)] = covSelection(Y, X, 'method', methodLabels{i}, 'm', m(n), 'c', cofS, 'B', B(b), 'trueCoeff', gamma);
			        		end;

				        	% compute predicted values
				        	XS = [ ones(size(X,1),1) X ];
				        	fhatN(mc,:,i,n) = XS*[interceptN(mc,i,n); gammahatN(mc,:,i,n)'];

							% compute prediction error
			        		RN(mc,i,n,b) = mean((Y-fhatN(mc,:,i,n)').^2);

				        	% compute MSE (i.e. prediction error divided by n)
				        	MSEN(mc,i,n,b) = RN(mc,i,n,b) / scriptN(n);
				        	
				        	% update waitbar
				        	if (~nodisplay)
					        	step = step + 1;
					        	waitbar(step / steps)
					        end;
				        end;
				    else
				    	MSEN(mc,:,n,b) = Inf; 
				    end;
			    end;

			    % pick the sample size yielding the smallest MSE 
			    [C nhatInd(mc,:)] = min(squeeze(MSEN(mc,:,:,b)),[],2); 

			    % store the gammahat, fhat, and R for the best sample size
			    for i=methodIndexStart:nMethods
			    	fhat(mc,:,i) = fhatN(mc,:,i,nhatInd(mc,i));
			    	gammahat(mc,:,i) = gammahatN(mc,:,i,nhatInd(mc,i));
			    	R(mc,i,b) = RN(mc,i,nhatInd(mc,i),b);
			    	MSE(mc,i,b) = MSEN(mc,i,nhatInd(mc,i),b);
			    	Shat = (1-(gammahat(mc,:,i)==0))';
			    	cost(mc,i,b) = c(Shat, scriptN(nhatInd(mc,i)), price);
			    end;

			    % compute the oracle risk
		    	Roracle(mc,b) = mean((Y-X*gamma).^2);

		    	% determine budgets necessary to attain experiment's MSE
			    if (strcmp(methods{1}, 'experiment'))
			    	LASSOind=1; PLASSOind=1;
			    	for i=2:nMethods
			    		switch methodLabels{i}
		        			case 'LASSO'
					        	bGap = @(Bval) budgetGap(MSE(mc,1,b), Y, X, Y, X, scriptN, price, 'method', methodLabels{i}, 'lambda', lambda(LASSOind), 'c', c, 'B', Bval, 'trueCoeff', gamma);
					        	LASSOind = LASSOind + 1;
					        case 'POST-LASSO'
					        	bGap = @(Bval) budgetGap(MSE(mc,1,b), Y, X, Y, X, scriptN, price, 'method', methodLabels{i}, 'lambda', lambda(PLASSOind), 'c', c, 'B', Bval, 'trueCoeff', gamma);
					        	PLASSOind = PLASSOind + 1;
			        		otherwise
			        			bGap = @(Bval) budgetGap(MSE(mc,1,b), Y, X, Y, X, scriptN, price, 'method', methodLabels{i}, 'm', m, 'c', c, 'B', Bval, 'trueCoeff', gamma);
		        		end;			    	

				    	if bGap(B) > 0
				    		Bupper = B; Blower = B/10;
				    	else
				    		Bupper = 10*B; Blower = B;
				    	end;
				    	eqb(mc,i) = bisect(Blower, Bupper, bGap, 5);

			    	end

			    end

		        % simulate the experiment
		        for i=1:nMethods

		        	% generate new data with treatment
		        	if i==1
		        		nSel = size(X,1);
		        	else
		        		nSel = scriptN(nhatInd(mc,i));
		        	end;
			        [Y, X, D] = generateData(nSel, M, model, 1, 'outcomeData', outcomeData, 'covariateData', covariateData);

			    	% estimate treatment effect    	
			        Ytilde = Y - X*gammahat(mc,:,i)';
			        Dt = [ones(length(D),1) D];
			        betahatc = pinv(Dt'*Dt)*Dt'*Ytilde;
			        betahat(mc,i) = betahatc(2);
			        sebetahat(mc,i) = sqrt(var(Ytilde-betahat(mc,i)*D) / (D'*D) / N);
			    end;
		end;
	end;
	toc;
	if (~nodisplay)
		close(h)
	end;


	% construct labels for methods
	methodLabelsLambda = { };
	for i=1:length(methods)
		if strcmp(methods{i}, 'LASSO') & (length(lambda)>1)
			for j=1:length(lambda)
				methodLabelsLambda = [ methodLabelsLambda, strcat('LASSO', num2str(lambda(j))) ];
			end;
		else
			methodLabelsLambda = [ methodLabelsLambda, methods{i} ];
		end;
	end;

end







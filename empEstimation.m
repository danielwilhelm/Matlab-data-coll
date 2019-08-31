TEMPLATE = 1;


% ------------------------------------------------------------------------
% parameter choices
% ------------------------------------------------------------------------


	% choose empirical application ("schoolgrants" or "childcare")
	if ~TEMPLATE
		application = 'childcare';
	end;

	% ============ methods for selection of covariates
 	methods = { 'OGA', 'LASSO', 'POST-LASSO' };

	% ============ k-fold out-of-sample evaluation
	% K = 1;	% K=1 means training and evaluation sample are the same (i.e. full sample)
	
	switch application

		case 'childcare'

			% ============ cost structure
			% S 		is Mx1 vector of ones and zeros, indicating selected regressors
			% price 	is Mx1 vector of prices extracted from childcare_description.xlsx
			% n 		is a scalar indicating the sample size
			% c = @(S,n,price) n*(price'*S);

			% cost features of the experiment
			TS_exp = 120;
			n_exp = 1466;
			sumS_exp = 40;
			cAdmin_exp = 10000;
			cTrain_exp = 25000;
			cInterv_exp = 630000;

			% assumed parameter values
			alpha = 0.4;
			eta = 200;

			% calibrated cost function
			phi = cAdmin_exp/TS_exp^alpha;
			kappa = @(n) 150*(n>0 & n<=1400) + (cTrain_exp/TS_exp)*(n>1400 & n<=3000) + 250*(n>3000 & n<=4500) + 300*(n>4500 & n<=6000) + 350*(n>6000);
			T = @(S,price) 1+price'*S;
			p = (cInterv_exp/n_exp - eta)/TS_exp;

			cadmin = @(S,price) phi*T(S,price)^alpha;
			ctrain = @(S,n,price) kappa(n)*T(S,price);
			cinterv = @(S,n,price) n*( eta + p*T(S,price) );
			c = @(S,n,price) cadmin(S,price) + ctrain(S,n,price) + cinterv(S,n,price);


			% ============ outcome variable (1=cognitive test, 2=health assessment)
			if ~TEMPLATE
				outcomeInd = 2;
			end;

			% ============ budget: 
			%	- if B is a scalar, then perform covariate selection only for this budget and report performance statistics instead of plot
			%	- set to NaN if budget should be set to cost of the experiment
			B = NaN;						
			% B = [1 1.2 1.5 1.7 2 5 10];		%	- if B is a vector, then perform covariate selection for each of those budgets and plot the prediction error as a function of B

			% ============ data collection parameters
			if ~TEMPLATE
				minN = 500;						% min sample size to be considered in the grid scriptN
				maxN = 3300;					% max sample size to be considered in the grid scriptN
				nN = 100;						% number of sample sizes in the grid scriptN
			end;

			% ============ LASSO penalization parameter (negative value: do LASSO that satisfies the budget)
			lambda = -1;

		case 'schoolgrants'

			if ~TEMPLATE
				% ============ outcome variable (1=French, 2=math, 3=oral)
					outcomeInd = 2;

				% ============ use baseline (=1) or follow-up (=0) outcomes?
					baselineOutcomes = 0;

				% ============ exclude high cost covariates (i.e. baseline outcomes)?
					dropAllHighCostCovariates = 0;

				% ============ force inclusion of lagged outcome?
					forceOwnBaselineOutcome = 0;

				% ============ factor that multiplies price of high cost covariates (i.e. 1 means no scaling)
					highCostPriceScalingFactor = 1;

				% ============ factor that multiplies the correlation of baseline and follow-up outcomes (i.e. 1 means no change, >1 means baseline outcome becomes more predictive)
					outcomeCorrScalingFactor = 1.5;
			end;

			% ============ for follow-up outcomes, run variable selection using only control group?
				controlGroupOnly = 1;


			% ============ budget: 
			%	- if B is a scalar, then perform covariate selection only for this budget and report performance statistics instead of plot
			%	- set to NaN if budget should be set to cost of the experiment
			B = NaN;
			% B = [1 1.2 1.5 1.7 2 5 10];		%	- if B is a vector, then perform covariate selection for each of those budgets and plot the prediction error as a function of B

			% ============ cost function
			% S 		is Mx1 vector of ones and zeros, indicating selected regressors
			% price 	is Mx1 vector of prices extracted from childcare_description.xlsx
			% n 		is a scalar indicating the sample size

			% cost features of the experiment
			TS_low_exp = 60;  TS_high_exp = 45;
            c_exp = 350; nc_exp = 24;
			sumS_low_exp = 254; sumS_high_exp = 3;
			cAdmin_low_exp = 5000; cAdmin_high_exp = 24000;
			cTrain_low_exp = 1600; cTrain_high_exp = 1600;
			cInterv_low_exp = 10000; cInterv_high_exp = 150000;
			cInterv_fixed_low_exp = 500; cInterv_fixed_high_exp = 1000;

			% assumed parameter values
			alpha_low = 0.7; alpha_high = 0.7;
			lambda_low = 1/7; lambda_high = 1/52.5;

			% calibrated cost function
			nc = nc_exp;
			phi_low = cAdmin_low_exp/TS_low_exp^alpha_low; 
			phi_high = cAdmin_high_exp/TS_high_exp^alpha_high;
			
			mu_low = @(n) floor(lambda_low*(n/nc));
			nlower = 0:10:60; nupper = 10:10:70; nseq = 1:7;
			mu_n_high = nseq(nlower < nc & nc <= nupper);
			mu_high = @(n) floor(lambda_high*(n/nc)*mu_n_high);

			nlower = 0:20:360; nupper = 20:20:380; nseq = 1:19;
			kappa_low = @(n) 20* nseq(nlower < mu_low(n) & mu_low(n) <= nupper);
			nlower = 0:4:64; nupper = 4:4:68; nseq = 1:17;
			kappa_high = @(n) 12* nseq(nlower < mu_high(n) & mu_high(n) <= nupper);
			
			psi_low = @(n) mu_low(n); psi_high = @(n) mu_high(n); 

			eta_low = cInterv_fixed_low_exp/mu_low(c_exp*nc_exp);
			eta_high = cInterv_fixed_high_exp/mu_high(c_exp*nc_exp);
			p_low = (cInterv_low_exp-cInterv_fixed_low_exp)/(c_exp*TS_low_exp);
			p_high = (cInterv_high_exp-cInterv_fixed_high_exp)/(c_exp*nc_exp*TS_high_exp);

			tau_high = TS_high_exp / (sumS_high_exp);

			if baselineOutcomes
				T_low = @(S,price) price'*S;
				T_high = @(S,price) tau_high*1;
			else
				highCostCovariates = 3:5;
				T_low = @(S,price) price(setdiff(1:length(S),highCostCovariates))'*S(setdiff(1:length(S),highCostCovariates));
				T_high = @(S,price) price(highCostCovariates(1))*1 + price(highCostCovariates)'*S(highCostCovariates);
			end;
			
			cadmin = @(S,price) 1*(phi_low*T_low(S,price)^alpha_low + phi_high*T_high(S,price)^alpha_high);
			ctrain = @(S,n,price) 1*(kappa_low(n)*T_low(S,price) + kappa_high(n)*T_high(S,price));
			cinterv = @(S,n,price) 1*((sum(S)>0)*(psi_low(n)*eta_low + n/nc*p_low*T_low(S,price)) + psi_high(n)*eta_high + n*p_high*T_high(S,price) );
			c = @(S,n,price) cadmin(S,price) + ctrain(S,n,price) + cinterv(S,n,price);

			% ============ data collection parameters
			if ~TEMPLATE
				minN = 500;						% minimum sample size to be considered in the grid scriptN
				maxN = 3200;					% max sample size to be considered in the grid scriptN
				nN = 10;						% number of sample sizes in the grid scriptN
			end;

			% ============ LASSO penalization parameter (negative value: do LASSO that satisfies the budget)
			% lambda = [ 0.01 0.05 0.1 ];
			lambda = -1;

			% ============ select household variables in addition to the others?
			selectHH0 = 0;
	end;

	% report RMSE as ratio?
	reportRatios = 0;

	% perform budget comparisons?
	if ~TEMPLATE
		budgetComparison = 1;
	end;


% ------------------------------------------------------------------------
% Initializations
% ------------------------------------------------------------------------

	switch application
		case 'childcare'
			outcomes = { 'test1', 'test4' };
			outcomeLabels = { 'cognitive test', 'health assessment' };

			dataFilename = 'data/childcare.txt';
			descriptionFilename = 'data/childcare_description.xlsx';
			descriptionCells = 'A2:A44';
			priceCells = 'B5:B44';

			outcomeCells = 1:2;
			covariateCells = 4:43;
			treatCell = 3;

			selectHH0 = 0;
			clusterDummies = 0;
			controlGroupOnly = 0;
			excludeReg = NaN;
			includeReg = NaN;

		case 'schoolgrants'
			outcomeLabels = { 'French test', 'math test', 'oral test' };

			dataFilename = 'data/schoolgrants.txt';
			descriptionFilename = 'data/schoolgrants_description.xlsx';
			descriptionCells = 'A2:A263';
			priceCells = 'B9:B263';
			pricesIDCell = 'B2:B3';

			if baselineOutcomes
				outcomes = { 'w0_testf', 'w0_testm', 'o0_total' };
				outcomeCells = 10:12;
				covariateCells = [ 8:9 13:262 ];
				excludeReg = NaN;
				includeReg = NaN;
			else
				outcomes = { 'w1_testf', 'w1_testm', 'o1_total' };
				outcomeCells = 3:5;
				covariateCells = [ 8:262 ];	
			end;
			treatCell = 6;
	end

	if (budgetComparison)
		methods = [ 'experiment' methods ];
	end;
	
	% construct outcome labels
	outcome = outcomes{outcomeInd};
	outcomeLabel = outcomeLabels{outcomeInd};

	% load data and variable names
	dat=dlmread(dataFilename);
	[~,~,varnames]=xlsread(descriptionFilename, 1, descriptionCells);
	Yorig = dat(:,outcomeCells); Yvarnames = varnames(outcomeCells);
	Xorig = dat(:,covariateCells); Xorigvarnames = varnames(covariateCells);
	treatorig = dat(:,treatCell);

	% load price data
	[~,~,priceOrig]=xlsread(descriptionFilename, 1, priceCells);
	priceOrig = cell2mat(priceOrig);

	% additional initializations for schoolgrants application:
	if (strcmp(application,'schoolgrants'))

		% load household variables starting with "hh0"?
		if (selectHH0)
			selXorig = Xorig;
			Xvarnames = Xorigvarnames;
		else
			selIndicators = ~strncmpi('hh0', Xorigvarnames, 3);
			selXorig = Xorig(:,selIndicators);
			priceOrig = priceOrig(selIndicators);
			Xvarnames = Xorigvarnames(selIndicators);
		end;
	else
		selXorig = Xorig;
		Xvarnames = Xorigvarnames;
	end;

	% find missing values in X and Y
	missingX = any(isnan(selXorig),2);
	outcomeIndicator = strcmpi(outcome, Yvarnames);
	selYorig = Yorig(:,outcomeIndicator);
	missingY = isnan(selYorig);
	notmissing = (~missingY) & (~missingX);

	% select observations without any missing values in Y or X
	X = selXorig(notmissing,:);
	Y = selYorig(notmissing);
	treat = treatorig(notmissing);

	% in schoolgrants application, use only control group?
	if (strcmp(application,'schoolgrants') & (~baselineOutcomes) & controlGroupOnly)
		X = X(find(treat==0),:);
		Y = Y(find(treat==0),:);
	end;

	% standardize covariates by their standard deviation and drop those covariates that have variance=0
	var0Indicators = (std(X)==0); N = size(X,1);
	selX = X(:,~var0Indicators);
	stdmat=ones(N,1)*std(selX);
	X = selX ./ stdmat;
	price = priceOrig(~var0Indicators);
	Xvarnames = Xvarnames(~var0Indicators);

	% construct labels for methods
	methodLabels = { };
	for i=1:length(methods)
		if sum(strcmp(methods{i}, { 'LASSO', 'POST-LASSO' })) & (length(lambda)>1)
			for j=1:length(lambda)
				methodLabels = [ methodLabels, methods{i} ];
			end;
		else
			methodLabels = [ methodLabels, methods{i} ];
		end;
	end;
	nMethods = length(methodLabels);
	if (strcmp(methods{1}, 'experiment'))
		methodIndexStart = 2;
	else
		methodIndexStart = 1;
	end;

	% set the budget equal to the cost in the experiment
	M = size(X,2);
	NK = floor(N/K); 
	if (K==1)
		Ntrain = N; Neval = N;
	else
		Ntrain = (K-1)*NK; Neval = NK;
	end;
	if isnan(B)
		B = c(ones(M,1), Ntrain, price);
	end;

	% adjustments for sensitivity analysis	
	if (strcmp(application,'schoolgrants'))
		if dropAllHighCostCovariates & forceOwnBaselineOutcome
			error('cannot drop and simultaneously force inclusion of baseline outcomes!')
		else
			if (~baselineOutcomes)
				highCostCovariatesNames = strrep(outcomes, '1', '0');

				% regressors to be forced for inclusion
				if dropAllHighCostCovariates
					excludeReg = zeros(M,1);
					for oname=1:length(highCostCovariatesNames)
						ind = strmatch(highCostCovariatesNames{oname}, Xvarnames, 'exact');
						excludeReg(ind) = 1;
					end;
				else
					excludeReg = NaN;
				end;

				if forceOwnBaselineOutcome
					includeReg = zeros(M,1);
					ownBaselineName = strrep(outcome, '1', '0');
					ind = strmatch(ownBaselineName, Xvarnames, 'exact');
					includeReg(ind) = 1;
				else
					includeReg = NaN;
				end;

				% price adjustments for high cost covariates (if any)
				price(highCostCovariates) = highCostPriceScalingFactor*price(highCostCovariates);

				% adjust correlation of baseline and follow-up outcomes
				if (outcomeCorrScalingFactor ~= 1)

					% regress follow up outcome on baseline outcome
					ownBaselineInd = strmatch(highCostCovariatesNames(outcomeInd), Xvarnames, 'exact');
					Ybase = [ ones(N,1) X(:,ownBaselineInd) ];
					bhat = inv(Ybase'*Ybase)*Ybase'*Y;
					ehat = Y-Ybase*bhat;
					shat2 = var(ehat);
					R2 = var(Ybase*bhat);
					cFactor = outcomeCorrScalingFactor;
					dFactor = sqrt((1-cFactor^2)*R2/shat2 + 1);
					Y = bhat(1) + bhat(2)*cFactor*X(:,ownBaselineInd) + dFactor*ehat;
				end;
			else
				excludeReg = NaN;
				includeReg = NaN;
			end;
		end;

		
	end;

	

	% other initializations
	switch application
		case 'schoolgrants'
			if (baselineOutcomes)
				bl = 'BASE';
			else
				bl = 'FOLLOW';
				if (dropAllHighCostCovariates)
					bl = 'FOLLOW_NOHIGH';
				end;
				if (forceOwnBaselineOutcome)
					bl = 'FOLLOW_FORCEBASE';
				end;
				if (highCostPriceScalingFactor ~= 1)
					bl = [ bl '_PRICE_ADJ' ];
				end;
				if (outcomeCorrScalingFactor ~= 1)
					bl = [ bl '_CORR_ADJ' ];
				end;	
			end;
			suffix = [ '_' bl ];
		case 'childcare'
			suffix = [ ];
	end;
	scriptN = floor(linspace(minN, maxN, nN));
	nB = length(B); 
	fhatN=zeros(K,Neval,nMethods,nN); gammahatN=zeros(M,K,nMethods,nN); interceptN=zeros(nMethods,nN); nhatInd=zeros(nMethods,1);
	fhat=zeros(K,Neval,nMethods); gammahat=zeros(M,K,nMethods); R=zeros(K,nMethods,nB); RN=zeros(K,nMethods,nN,nB); MSE=RN; eqb=zeros(K,nMethods,1);
	nhatInd=zeros(nMethods,K);


% ------------------------------------------------------------------------
% perform covariate selection
% ------------------------------------------------------------------------

	% loop through all budgets
	tic;
	for iK=1:K
		for b=1:length(B)

			% compute the number of covariates to be collected
			m = ones(1,length(scriptN));

			% permutation for random benchmark choice of covariates
			perm = randsample(M, max(m));

			% if K>1, select training and evaluation samples
			if (K==1)
				evalInd = 1:N; trainInd = 1:N;
			else
				evalInd = ((iK-1)*NK+1):(iK*NK); train = 1:(K*NK); trainInd = train(train<((iK-1)*NK+1) | train>(iK*NK));
			end;
			Xtrain = X(trainInd,:);
			Ytrain = Y(trainInd);
			Xeval = X(evalInd,:);
			Yeval = Y(evalInd);

			% estimate model as in experiment
	        if (strcmp(methods{1}, 'experiment'))
		        XStrain = [ ones(Ntrain,1) Xtrain ];
				coeffi = pinv(XStrain'*XStrain)*XStrain'*Ytrain;
				gammahat(:,iK,1) = coeffi(2:end);
				XSeval = [ ones(Neval,1) Xeval ];
				fhat(iK,:,1) = XSeval*coeffi;
				R(iK,1,b) = mean((Yeval-fhat(iK,:,1)').^2);
				MSEmin(iK,1,b) = R(iK,1,b)/Ntrain;
				cost(1,b) = c(ones(M,1), Ntrain, price);
				eqb(iK,1) = cost(1,b);
			end;

		    % loop through sample sizes in scriptN
		    parfor n=1:nN

		    	% define cost function given sample size
				cofS = @(S) c(S,scriptN(n),price);

		        % loop through different covariate selection methods
		        LASSOind=1; PLASSOind=1;
		        for i=methodIndexStart:nMethods
		        	
		        	% estimate regression coefficients by covariate selection method
		        	if m(n)>0
		        		switch methodLabels{i}
		        			case 'LASSO'
					        	[ gammahatN(:,iK,i,n) interceptN(i,n) ] = covSelection(Ytrain, Xtrain, 'method', methodLabels{i}, 'lambda', lambda(LASSOind), 'c', cofS, 'B', B(b), 'includeReg', includeReg, 'excludeReg', excludeReg);
					        	LASSOind = LASSOind + 1;
					        case 'POST-LASSO'
					        	[ gammahatN(:,iK,i,n) interceptN(i,n) ] = covSelection(Ytrain, Xtrain, 'method', methodLabels{i}, 'lambda', lambda(PLASSOind), 'c', cofS, 'B', B(b), 'includeReg', includeReg, 'excludeReg', excludeReg);
					        	PLASSOind = PLASSOind + 1;
			        		otherwise
			        			[ gammahatN(:,iK,i,n) interceptN(i,n) ] = covSelection(Ytrain, Xtrain, 'method', methodLabels{i}, 'm', m(n), 'c', cofS, 'B', B(b), 'perm', perm, 'includeReg', includeReg, 'excludeReg', excludeReg);
		        		end;
		        	end;

		        	% compute predicted values
		        	XSeval = [ ones(Neval,1) Xeval];
		        	fhatN(iK,:,i,n) = XSeval*[interceptN(i,n); gammahatN(:,iK,i,n)];

					% compute prediction error
		    		RN(iK,i,n,b) = mean((Yeval-fhatN(iK,:,i,n)').^2);

		        	% compute MSE (i.e. prediction error divided by n)
		        	if (cofS(gammahatN(:,iK,i,n)>0)<= B)
			        	MSE(iK,i,n,b) = RN(iK,i,n,b) / scriptN(n);
			        else
			        	MSE(iK,i,n,b) = Inf;
			        end;

		        end;

		    end;

		    % pick the sample size yielding the smallest MSE 
		    [C nhatInd(:,iK)] = min(squeeze(MSE(iK,:,:,b)),[],2); 

		    % store the gammahat, fhat, and R for the best sample size
		    for i=methodIndexStart:nMethods
		    	fhat(iK,:,i) = fhatN(iK,:,i,nhatInd(i,iK));
		    	gammahat(:,iK,i) = gammahatN(:,iK,i,nhatInd(i,iK));
		    	MSEmin(iK,i,b) = MSE(iK,i,nhatInd(i,iK),b);
		    end;

		    % determine budgets necessary to attain experiment's MSE
	        if (strcmp(methods{1}, 'experiment'))
		    	LASSOind=1; PLASSOind=1;
		    	for i=2:nMethods
		    		switch methodLabels{i}
	        			case 'LASSO'
				        	bGap = @(Bval) budgetGap(MSEmin(iK,1,b), Ytrain, Xtrain, Yeval, Xeval, scriptN, price, 'method', methodLabels{i}, 'lambda', lambda(LASSOind), 'c', c, 'B', Bval, 'includeReg', includeReg, 'excludeReg', excludeReg);
				        	LASSOind = LASSOind + 1;
				        case 'POST-LASSO'
				        	bGap = @(Bval) budgetGap(MSEmin(iK,1,b), Ytrain, Xtrain, Yeval, Xeval, scriptN, price, 'method', methodLabels{i}, 'lambda', lambda(PLASSOind), 'c', c, 'B', Bval, 'includeReg', includeReg, 'excludeReg', excludeReg);
				        	PLASSOind = PLASSOind + 1;
		        		otherwise
		        			bGap = @(Bval) budgetGap(MSEmin(iK,1,b), Ytrain, Xtrain, Yeval, Xeval, scriptN, price, 'method', methodLabels{i}, 'm', m, 'c', c, 'B', Bval, 'includeReg', includeReg, 'excludeReg', excludeReg);
	        		end;			    	

			    	if bGap(B) > 0
			    		Bupper = B; Blower = B/500;
			    	else
			    		Bupper = 500*B; Blower = B;
			    	end;
			    	eqb(iK,i) = bisect(Blower, Bupper, bGap, 15);

		    	end

		    end
			
		end;
	end;
	toc;



% ------------------------------------------------------------------------
% show results
% ------------------------------------------------------------------------
	
	% find missing values in treatment indicator
	treatmissing = isnan(treat);

	% treatment effect in the experiment w/o covariates
	if (~controlGroupOnly)
		Xexp = [ ones(N-sum(treatmissing),1) treat(~treatmissing) ];
		betahatc = pinv(Xexp'*Xexp)*Xexp'*Y(~treatmissing);
		betahat = betahatc(2);
		disp(' '); disp(' ');
		disp(strcat('treatment effect w/o covariates:  ', num2str(betahat)));
	end;
	

	% create directory for results
	rfolder = datestr(now,'yyyymmdd');
	mkdir(['results-' application], rfolder)


	% construct labels for methods including multiple lambdas in case of LASSO
	methodLabelsLambda = { };
	for i=1:length(methods)
		if sum(strcmp(methods{i}, { 'LASSO', 'POST-LASSO' })) & (length(lambda)>1)
			for j=1:length(lambda)
				methodLabelsLambda = [ methodLabelsLambda, strcat(methods{i}, num2str(lambda(j))) ];
			end;
		else
			methodLabelsLambda = [ methodLabelsLambda, methods{i} ];
		end;
	end;

	if nB==1

		% print the set scriptN and the corresponding max no. of covariates
		disp(' '); disp('price of collecting covariate per person:'); disp(dataset({price 'price'}, 'obsnames', Xvarnames));


		% display summary statistics
		nhat = zeros(1,nMethods);
		if K==1
			nSelCov = squeeze(sum(1-(gammahat==0),1));
		else
			nSelCov = mean(squeeze(sum(1-(gammahat==0),1)),1);
		end;
		
		if (budgetComparison)
			nhat(1) = Ntrain;
			if K==1
				nhat(2:end) = squeeze(scriptN(nhatInd(2:end,:)));
			else
				nhat(2:end) = mean(scriptN(nhatInd(2:end,:)), 2);
			end;
		else
			if K==1
				nhat = squeeze(scriptN(nhatInd));
			else
				nhat = mean(scriptN(nhatInd), 2);
			end;
		end;
		costToBK = zeros(K,nMethods);
		for iK=1:K
			for i=1:nMethods
				if i==1
					costToBK(iK,i) = 1;
				else
					costToBK(iK,i) = c((1-(gammahat(:,iK,i)==0)), scriptN(nhatInd(i,iK)), price)/B;
				end;
			end;	
		end;
		if K==1
			summaryTable = [ nhat' nSelCov squeeze(costToBK)' sqrt(squeeze(MSEmin))' squeeze(eqb)' ones(nMethods,1) ];
		else
			summaryTable = [ nhat' nSelCov' mean(costToBK,1)' sqrt(mean(MSEmin,1))' mean(eqb,1)' ones(nMethods,1) ];
		end;

		for i=2:nMethods
			summaryTable(i,end) = summaryTable(i,end-1) / summaryTable(1,end-1);
		end;
		

		disp(' '); disp('--------------- Results ---------------')
		disp(' '); disp(dataset({summaryTable 'n','selection','costratio','RMSE','EQB','relEQB'}, 'obsnames', methodLabelsLambda)); disp(' ');
		disp('where:');
		disp(' ');
		disp('     n:           selected sample size');
		disp('     selection:   no. of selected covariates #(S)');
		disp('     costratio:   the ratio of cost, i.e. c(S,n,price), divided by the budget B');
		disp('     RMSE:        objective function');
		disp('     EQB:         equivalent budget');
		disp('     relEQB:      equivalent budget relative to experiment');
		disp(' ');

	else
		figure('Position',[1,1,800,800]);
		marker = {'+-','o-','*-','x-','s-','d-','^-','v-','>-','<-','p-','h-'};
		plot(B, R(1,:), '-', 'Linewidth', 2, 'Color', [0 0 1])
		hold on
		for i=2:nMethods
			plot(B, R(i,:), marker{i}, 'Color', [0.6 0.6 0.6]);
		end;
		legend(methodLabelsLambda);
		axis 'auto y';
		xlabel('B');
		title('prediction error');
		xlim([min(B) max(B)]);
		ylim('auto');
	end;

	% display selected covariates
	if K==1
		disp(' '); disp(' '); disp('selected covariates:'); disp(' ');
		for i=methodIndexStart:nMethods
			disp([ '----------- ' methodLabelsLambda{i} ' -----------']);
			disp(Xvarnames(find(squeeze(gammahat(:,1,i))))); disp(' ');
		end;

		% save list of selected covariates to latex file
		maxNSelCov = max(nSelCov);
		varList = cell(maxNSelCov,nMethods);
		for i=1:nMethods
			newCol = Xvarnames(find(squeeze(gammahat(:,1,i))));
			for k=1:length(newCol)
				varList(k,i) = newCol(k);
			end;
		end;
		varList(cellfun('isempty',varList))={' '};
		
		rfilename = [ 'results-' application '/' rfolder '/Table' tableno '_varlist.tex'];
		addText = @(s) [ '\text{' strrep(s, '_', '\_') '}' ];
		varList = cellfun(addText, varList, 'uniformoutput', 0);
		matrix2latex(varList, rfilename, 'headerRow', methodLabelsLambda, 'alignment', 'l', 'format', '%10s', 'caption', [ application ' (outcome: ' outcomeLabel ')'], 'size', 'tiny' );
	end;

	% save summary table to latex file
	if (reportRatios)
		reportVars = { 'Method', '$\hat{n}$','$|\hat{I}|$','Cost/B', 'RMSE ratio', 'EQB' , 'Relative EQB'};
	else
		reportVars = { 'Method', '$\hat{n}$','$|\hat{I}|$','Cost/B', 'RMSE', 'EQB', 'Relative EQB' };
	end;
	rfilename = [ 'results-' application '/' rfolder '/Table' tableno '.tex'];
	if K==1
		fc = { '$%d$', '$%d$', '$%.5f$', '$%.5f$', '$%.2f$', '$%.2f$' };
	else
		fc = { '$%.1f$', '$%.1f$', '$%.5f$', '$%.5f$', '$%.2f$', '$%.2f$' };
	end;
	matrix2latex(summaryTable, rfilename, 'headerRow', reportVars, 'headerColumn', methodLabelsLambda, 'formatColumn', fc, 'alignment', 'r', 'caption', [ application ' (outcome: ' outcomeLabel ')']);
	disp(['results are saved in ' rfilename]);


	% report summary statistics about objective function across different n
	if K==1
		costv = zeros(nN,1);
		costvAdmin = zeros(nN,1); costvTrain = zeros(nN,1); costvInterv = zeros(nN,1);
		for k=1:nN
			cval = c((squeeze(gammahatN(:,1,2,k))>0), scriptN(k), price);
			cvalAdmin = cadmin((squeeze(gammahatN(:,1,2,k))>0), price);
			cvalTrain = ctrain((squeeze(gammahatN(:,1,2,k))>0), scriptN(k), price);
			cvalInterv = cinterv((squeeze(gammahatN(:,1,2,k))>0), scriptN(k), price);
			if isempty(cval)
				costv(k) = NaN; costvAdmin(k) = NaN; costvTrain(k) = NaN; costvInterv(k) = NaN;
			else
				costv(k) = cval; costvAdmin(k) = cvalAdmin; costvTrain(k) = cvalTrain; costvInterv(k) = cvalInterv;
			end;
		end;
		gamat = squeeze(gammahatN(:,1,2,:));
		objTable = [ squeeze(MSE(1,2,:)), scriptN', (sum(~(gamat==0 | gamat==Inf),1))' costv/B costvAdmin./costv costvTrain./costv costvInterv./costv];
		objTable = objTable(objTable(:,1)<Inf,:);
		disp(' '); disp(dataset({objTable 'MSE','ngrid','selcov','costratio','costratioAdmin','costratioTrain','costratioInterv'})); disp(' ');
	end;
	

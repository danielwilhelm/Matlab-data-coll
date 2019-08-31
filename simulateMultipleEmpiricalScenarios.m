
	rng('default');
	rng(110);


% ------------------------------------------------------------------------
% simulation parameters
% ------------------------------------------------------------------------

	% number of MC samples
	% nMC = 10;

	% choose empirical application ("schoolgrants" or "childcare")
	% application = 'childcare';

	% ============ methods for selection of covariates ('OGA', 'largest', 'random', 'LASSO', 'POST-LASSO')
	methods = { 'OGA', 'LASSO', 'POST-LASSO' };
	
	switch application

		case 'childcare'

			% ============ cost structure
			% S 		is Mx1 vector of ones and zeros, indicating selected regressors
			% price 	is Mx1 vector of prices extracted from childcare_description.xlsx
			% n 		is a scalar indicating the sample size
			% c = @(S,n,price) n*(price'*S);

			alpha = 0.4;
			phi = 10000/120^alpha;
			kappa = @(n) 150*(n>0 & n<=1400) + (25000/120)*(n>1400 & n<=3000) + 250*(n>3000 & n<=4500) + 300*(n>4500 & n<=6000) + 350*(n>6000);
			T = @(S,price) 1+price'*S;
			eta = 200;
			p = (630000/1466 - eta)/120;

			cadmin = @(S,price) phi*T(S,price)^alpha;
			ctrain = @(S,n,price) kappa(n)*T(S,price);
			cinterv = @(S,n,price) n*( eta + p*T(S,price) );
			c = @(S,n,price) cadmin(S,price) + ctrain(S,n,price) + cinterv(S,n,price);


			% ============ outcome variable (1=cognitive test, 2=health assessment)
			outcomeInd = 1;

			% ============ budget: 
			B = 569074.0637880735;								%	- if B is a scalar, then perform covariate selection only for this budget and report performance statistics instead of plot
			% B = [1 1.2 1.5 1.7 2 5 10];		%	- if B is a vector, then perform covariate selection for each of those budgets and plot the prediction error as a function of B

			% ============ data collection parameters
			minN = 1000;						% min sample size to be considered in the grid scriptN
			maxN = 4000;					% max sample size to be considered in the grid scriptN
			nN = 50;						% number of sample sizes in the grid scriptN

			% ============ LASSO penalization parameter (negative value: do LASSO that satisfies the budget)
			lambda = -1;


		case 'schoolgrants'

			% ============ cost function
			% S 		is Mx1 vector of ones and zeros, indicating selected regressors
			% price 	is Mx1 vector of prices extracted from childcare_description.xlsx
			% n 		is a scalar indicating the sample size
			c = @(S,n,price) n*(price'*S);

			% ============ budget: 
			B = 1;								%	- if B is a scalar, then perform covariate selection only for this budget and report performance statistics instead of plot
			% B = [1 1.2 1.5 1.7 2 5 10];		%	- if B is a vector, then perform covariate selection for each of those budgets and plot the prediction error as a function of B

			% ============ data collection parameters
			minN = 50;						% minimum sample size to be considered in the grid scriptN
			maxN = 1000;					% max sample size to be considered in the grid scriptN
			nN = 10;						% number of sample sizes in the grid scriptN

			% ============ outcome variable (1=French, 2=math, 3=oral)
			outcomeInd = 1;

			% ============ LASSO penalization parameter (negative value: do LASSO that satisfies the budget)
			lambda = -1;

			% ============ select only household variables?
			selectHH0 = 1;
	end;

	% report ratios of RMSE?
	reportRatios = 0;

	% order covariates by decreasing correlation with outcome?
	orderCovariates = 1;

	% perform budget comparisons?
	budgetComparison = 1;


% ------------------------------------------------------------------------
% define scenarios
% ------------------------------------------------------------------------
	

	% decay pattern for the regression coefficients ('exp', 'lin-sparse', 'lin-exp', ...; see getCoefficients.m for possible specifications)
	% model.coefficients = { 'exp' };

	% scale factor for regression coefficients
	model.scale = { 0, 0.3, 0.7, 1 };

	% correlation pattern of the covariates
	model.correlations = 'empirical';

	% treatment effect
	model.beta = 0.18656;


% ------------------------------------------------------------------------
% load data
% ------------------------------------------------------------------------

	switch application
		case 'childcare'
			outcomes = { 'test1', 'test4' };
			outcomeLabels = { 'cognitive test', 'health assessment' };

			dataFilename = 'data/childcare.txt';
			descriptionFilename = 'data/childcare_description.xlsx';
			descriptionCells = 'A2:A44';
			priceCells = 'B5:B44';

			covariateCells = 4:43;
			outcomeCells = 1:2;

			selectHH0 = 0;

		case 'schoolgrants'
			outcomes = { 'w1_testf', 'w1_testm', 'o1_total' };
			outcomeLabels = { 'French test', 'math test', 'oral test' };

			dataFilename = 'data/schoolgrants.txt';
			descriptionFilename = 'data/schoolgrants_description.xlsx';
			descriptionCells = 'A2:A263';
			priceCells = 'B8:B263';

			covariateCells = 7:262;
			outcomeCells = 3:5;
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

	% in schoolgrants application: if selectHH0=1, then only use the household variables starting with 'hh0'
	if (selectHH0 & strcmp(application,'schoolgrants'))
		selIndicators = strncmpi('hh0', Xorigvarnames, 3);
		selXorig = Xorig(:,selIndicators);
		Xvarnames = Xorigvarnames(selIndicators);
	else
		selXorig = Xorig;
		Xvarnames = Xorigvarnames;
	end;

	% load price data
	[~,~,priceOrig]=xlsread(descriptionFilename, 1, priceCells);
	priceOrig = cell2mat(priceOrig);

	% find missing values in X and Y
	missingX = any(isnan(selXorig),2);
	outcomeIndicator = strcmpi(outcome, Yvarnames);
	selYorig = Yorig(:,outcomeIndicator);
	missingY = isnan(selYorig);
	notmissing = (~missingY) & (~missingX);

	% select observations without any missing values in Y or X
	X = selXorig(notmissing,:);
	Y = selYorig(notmissing);

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

	% sample size
	N = length(Y);

	% order the covariates by their correlation with the outcome
	if orderCovariates
		cct = corrcoef([Y X]);
		[cc I] = sort(cct(2:end,1),'descend');
		X = X(:,I);
	end;


% ------------------------------------------------------------------------
% final initializations
% ------------------------------------------------------------------------

	% other initializations
	M = size(X,2);
	scriptN = floor(linspace(minN, maxN, nN));
	reportVars = { 'scale','M','method','$\hat{n}$','$|\hat{I}|$','cost/B', '$\sqrt{\widehat{MSE}_{\hat{n},N}(\hat{\mathbf{f}})}$', '$bias(\hat{\beta})$', '$sd(\hat{\beta})$', 'RMSE($\hat{\beta}$)', 'EQB' }; nVar = length(reportVars);
	nCoeff = length(model.coefficients); nScale = length(model.scale); nM = length(M);
	nScen = nCoeff*nScale*nM; scen=0; SummaryTable = zeros(nCoeff, nScale, nM, nMethods, nVar);

	nodisplay = 1;


% ------------------------------------------------------------------------
% perform Monte Carlo simulations
% ------------------------------------------------------------------------

	for coeff=1:nCoeff

		parfor scale=1:nScale

			for Mind=1:nM

				modelParameters = struct('coefficients', model.coefficients{coeff}, 'correlations', model.correlations, 'beta', model.beta, 'scale', model.scale{scale});
				[gammahat, nSelCov, fhat, nhatInd, R, Roracle, MSE, betahat, sebetahat, m, err, methodLabels, cost, eqb] = simulate(nMC, N, M(Mind), modelParameters, methods, lambda, scriptN, price, c, B, 'nodisplay', nodisplay, 'outcomeData', Y, 'covariateData', X);
				
				% keep only those MC samples for which there were no errors
				eqb(:,1) = cost(:,1);
				nhat=zeros(nMC,nMethods);
				if (sum(err)>0)
					fprintf('\n In %i out of %i MC samples some error occurred!\n\n', sum(err), nMC);
					noerr = find(1-err);
					fhat=fhat(noerr,:,:); gammahat=gammahat(noerr,:,:); R=R(noerr,:); Roracle=Roracle(noerr,:); nhatInd=nhatInd(noerr,:); cost=cost(noerr,:); nSelCov=nSelCov(noerr,i);
					nhat=nhat(noerr,:); betahat=betahat(noerr,:); eqb=eqb(noerr,:);
				end;

				% write results into table
				for k=1:nMethods
					if (sum(err)<(length(err)-2))
						nSelCov(:,k) = sum(1-(gammahat(:,:,k)==0),2);
						if (strcmp(methods{k}, 'experiment'))
							nhat(:,k) = ones(nMC,1)*N;
						else
							nhat(:,k) = scriptN(nhatInd(:,k));
						end;
						costToB = mean(cost(:,k)/B);
						SummaryTable(coeff,scale,Mind,k,:) = [ scale M(Mind) k mean(nhat(:,k)) mean(nSelCov(:,k)) costToB mean(sqrt(R(:,k)./nhat(:,k))) mean(betahat(:,k)-model.beta) std(betahat(:,k)) sqrt(mean((betahat(:,k)-model.beta).^2)) mean(eqb(:,k))];
					else
						SummaryTable(coeff,scale,Mind,k,:) = [ NaN ];
					end;
				end;

			end;

		end;

	end;



% ------------------------------------------------------------------------
% save results
% ------------------------------------------------------------------------


disp(' '); disp([ 'simulations based on empirical distribution of covariates in ' application ' example']);
disp(' '); disp([ 'there are ' num2str(size(X,2)) ' covariates']);

% print the set scriptN and the corresponding max no. of covariates
disp(' '); disp('set of potential sample sizes:'); disp(scriptN);
disp(' '); disp('price per covariate:'); disp(price');

% construct labels for methods
methodLabelsLambda = { };
for i=1:length(methods)
	if strcmp(methods{i}, 'LASSO') & (length(lambda)>1)
		for j=1:length(lambda)
			methodLabelsLambda = [ methodLabelsLambda, strcat('LASSO ($\lambda=', num2str(lambda(j)), '$)') ];
		end;
	else
		methodLabelsLambda = [ methodLabelsLambda, methods{i} ];
	end;
end;

for coeff=1:nCoeff
	disp(model.coefficients{coeff});
	outputTable = { }; row = 1;
	
	for scale=1:nScale
		firstTime = 1;
		for Mind=1:nM
			for i=1:nMethods

				% change R_n and MSE columns to ratios
				if (reportRatios)
					if i==1
						send = SummaryTable(coeff,scale,Mind,1,end); SummaryTable(coeff,scale,Mind,i,end) = 1;
						send2 = SummaryTable(coeff,scale,Mind,1,end-2); SummaryTable(coeff,scale,Mind,i,end-2) = 1;
					else
						SummaryTable(coeff,scale,Mind,i,end) = SummaryTable(coeff,scale,Mind,i,end) / send;
						SummaryTable(coeff,scale,Mind,i,end-2) = SummaryTable(coeff,scale,Mind,i,end-2) / send2;
					end;
				end;
				
				% assemble the table
				summCell = num2cell(squeeze(SummaryTable(coeff,scale,Mind,i,4:end)))';
				summCell = cellfun(@num2str, summCell, 'uniformoutput', 0);
				if firstTime
					outputTable(row,:) = [ num2str(model.scale{scale}) num2str(M(Mind), '%i') methodLabelsLambda{i} summCell ];
					firstTime = 0;
				else
					outputTable(row,:) = [ NaN num2str(M(Mind), '%i') methodLabelsLambda{i} summCell ];
				end;
				row = row + 1;
			end;
		end;
	end;

	disp(' '); disp('--------------- Results ---------------');
	disp(' '); disp(dataset({squeeze(outputTable) 'scale','M','method','n','selection','costratio','sqrtR_n','bias','sd','RMSE','EQB'})); disp(' ');
	disp('where (all numbers are averaged across MC samples):');
	disp(' ');
	disp('     n:           selected sample size');
	disp('     selection:   no. of selected covariates #(gamma)');
	disp('     costratio:   the ratio of cost divided by the budget B');
	disp('     sqrtR_n:     sqrt of prediction error divided by n');
	disp('     bias:        bias of betahat');
	disp('     sd:          standard deviation of betahat');
	disp('     RMSE:        RMSE of betahat');
	disp('     EQB:         equivalent budget');
	disp(' ');

	
	% create directory for results
	rfolder = datestr(now,'yyyymmdd');
	mkdir(['results-sim'], rfolder)

	% save table to latex file
	rfilename = [ 'results-sim/' rfolder '/Table' tableno '.tex'];
	matrix2latex(outputTable, rfilename, 'headerRow', reportVars, 'alignment', 'r', 'caption', ['coefficients: ' model.coefficients{coeff}], 'size', 'tiny' );
	disp(['results are saved in ' rfilename]);
end;



% plot coefficients
figure('Position',[1,1,800,800]);
marker = {'+-','o-', 'x-','s-','d-','^-','v-','>-','<-','p-','h-'};

coeff0hat = (X'*X)^(-1)*X'*Y;
plot(1:M, coeff0hat, '-', 'Linewidth', 2, 'Color', [0 0 1])
hold on
for i=1:nCoeff
	coeffSim = coeff0hat + sign(coeff0hat).*getCoefficients(M, model.coefficients{i}, model.scale{4});
	plot(1:M, coeffSim, marker{i}, 'Linewidth', 1, 'Color', [0.2 0.2 0.2], 'MarkerSize', 6);
end;
plot(1:M, coeff0hat, '-', 'Linewidth', 2, 'Color', [0 0 1])

legend([ 'data' model.coefficients ]);
xlabel(' ');
title('regression coefficients');

% save plot
tightfig;
saveas(gcf, ['results-sim/' rfolder '/FigureS1.pdf'], 'pdf');



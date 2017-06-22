%% An Example of a Clustered Cost Function
% Daniel Wilhelm & Nicolas Cerkez
% University College London, 2017

function [c]=cluster_cfun(S, ss)
    % Input: 
        % S = a vector containing zeros or ones indicating which covariates to include
        % ss = sample size
    % Ouput: 
        % c = the cost

    % cost features of the experiment (to be set by researcher) 
    TS_exp = 120;           % assume the total survey cost is 120 (in units of time)
    n_exp = 1466;           % assume sample size is 1466
    cAdmin_exp = 10000;     % assume administrative costs are 10'000
    cTrain_exp = 25000;     % assume training costs are 25'000 
    cInterv_exp = 630000;   % assume interview costs are 630'000
    c_exp = 350;            % assume 350 clusters
    nc_exp = 24;            % assume sample size per cluster is 24

    % assumed parameter values
    alpha = 0.4;        % see paper (to be set by researcher)
    eta = 200;          % see paper (to be set by researcher)
    tau0 = 1;           % define cost of collecting outcome variable (to be set by researcher)
    tau = ones(1, 10);  % define a 1xM vector with the jth element being the cost of the jth covariate (to be set by reseacher) (suppose it is one per covariate)
    lambda = 1/7;       % see paper (to be set by researcher)

    % calibrated cost function
    phi = cAdmin_exp/TS_exp^alpha;              % see paper
    muc = lambda*c_exp;                         % see paper 
    mun = @(ss) 150*(ss>0 & ss<=1400) + (cTrain_exp/TS_exp)*(ss>1400 & ss<=3000) + 250*(ss>3000 & ss<=4500) + 300*(ss>4500 & ss<=6000) + 350*(ss>6000); % see paper (the various cut off points are to be set by the researcher)
    mu = abs(muc*mun(ss));                      % see paper
    kappa = @(ss) 150*(mu>0 & mu<=1400) + (cTrain_exp/TS_exp)*(mu>1400 & mu<=3000) + 250*(mu>3000 & mu<=4500) + 300*(mu>4500 & mu<=6000) + 350*(mu>6000); % see paper (the various cut off points are to be set by the researcher)
    T = @(S) tau0+tau*S;                        % see paper (survey time cost)
    p = (cInterv_exp/n_exp - eta)/TS_exp;       % see paper 

    cadmin = @(S) phi*T(S).^alpha;              % administrative costs
    ctrain = @(S,ss) kappa(ss)*T(S);            % training costs
    cinterv = @(S,ss) c_exp*eta + ss*p*T(S);    % interview costs (I here assume that psi(c) = c; this assumption may be changed)
    
    
    c =  cadmin(S) + ctrain(S,ss) + cinterv(S,ss);
        
end
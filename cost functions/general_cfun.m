%% An Example of a General Cost Function
% Daniel Wilhelm & Nicolas Cerkez
% University College London, 2017

function [c] = general_cfun(S, ss)
    % Create a function with the following input and output arguments: 
    % Inputs:  S = indices indexing which covariates to include
    %          ss = number of observations (sample size)
    % Output:  c = the cost function to be used
    
    % This cost function is of the following linear form: 
    %               c(S, ss) = c_admin + c_train + c_interv
    % where c_admin are administration costs, c_train are training costs,
    % and c_interv are interview costs. 
    % We start by specifying survey time costs (T). let tau(j) = 1..M be
    % the costs of collecting variable j for one individual, measured in
    % units of survey time. let tau0 denote the costs of collecting the
    % outcome variable, measured in units of suvery time. Then, the total
    % time costs of surveying one individual to elicit the variables
    % indicate by S are: T = tau0 + sum(tau(j)S(j))
    M = length(S);
    for i=1:M
        T = tau0 + tau(i)*S(i);
    end
    %kappa = ss;
    c_admin = phi*T^alpha;
    c_train = kappa * T^alpha;
    c_interv = ss(eta + p*T);
    
    c = c_admin + c_train + c_interv;
    
end

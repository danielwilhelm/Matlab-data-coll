clear all; close all; clc;
application='schoolgrants';
outcomeInd=2;
tableno = '8d';
minN=200;
maxN=10000;
nN=300;
budgetComparison=1;

baselineOutcomes=0;
dropAllHighCostCovariates = 0;
forceOwnBaselineOutcome = 0;
highCostPriceScalingFactor = 1;
outcomeCorrScalingFactor = 1.3;
K=1;

empEstimation

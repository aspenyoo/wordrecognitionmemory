function [bestFitParam, nLL, startTheta, Output] = paramfit_patternbayes(modelname, binningfn, nnew_part, nold_part, fixparams,nStartVals, nConf)
% 
% paramfit_patternbayes(MODELNAME) uses bayesian patternsearch (one of Luigi
% Acerbi's optimization algorithms) to find the best fit parameters (and
% coinciding negative loglikelihoods (nLLs) for MODELNAME of one subjects
% data.
%
% ===== INPUT VARIABLES =====
% MODELNAME: 'FP','FPheurs','VP','VPheurs',or 'uneqVar'
% ISLOGBINNING: 1: logistic. 0: linear
% NNEW_PART: 1x20 vector of responses for new distribution (total 150)
% NOLD_PART: 1x20 vector of responses for old distribution (total 150)
% FIXPARAMS: a 2xn matrix in which first row corresponds to index of which
% parameter to fix and second row corresponds to the value of that
% parameter.
% NSTARTVALS: number of starting thetas for any given iteration over M
% (default 10)
% ITERMS: vector of M values that will be tested
%
% ===== OUTPUT VARIABLE =====
% BESTFITPARAM: 1xn vector, where n is the number of unfixed parameters
% NLL: negative log likelihood of best fitting parameter values
% STARTTHETA: the starting theta value of the best fit parameter (useful
% when investigating if fminsearch is exploring very far or staying local).
%
%
% Aspen Yoo -- January 28, 2016
%

if nargin < 5; fixparams = []; end
if nargin < 6; nStartVals = 1; end
if nargin < 7; nConf = 20; end
nX = 30;
nS = 20;

% random number generator
rng('shuffle');

options = bps('defaults');              % Default options
options.UncertaintyHandling = 1;        % Activate noise handling
options.NoiseSize = 1;                  % Estimated noise magnitude

% for iM = 1:length(iterMs);
%     M = iterMs(iM);
%     if (M); display(M); end
%
%     if ~isempty(fixparams) && ~sum(iterMs==0);
%         fixparams = [1; M];
%     end

startTheta = genStartTheta;
bfp = nan(size(startTheta));
nll = nan(nStartVals,1); exitflag = nll;
for istartval = 1:nStartVals;
    obj_func = @(x) nLL_approx_vectorized(modelname, x, binningfn, nnew_part, nold_part, fixparams, nX, nS, nConf );
    [bfp(istartval,:) ,nll(istartval), exitflag(istartval), outputt{istartval}] = bps(obj_func,startTheta(istartval,:),lb,ub,plb,pub,options);
end

nll(exitflag < 0) = Inf;

nLL = min(nll);
bestFitParam =  bfp((nll == min(nll)),:);
Output = outputt{nll == min(nll)};
if ~isempty(fixparams);
    nParams = size(fixparams,2) + length(bestFitParam);
    unfixparams = 1:nParams;
    unfixparams(fixparams(1,:)) = [];
    bestFitParam = nan(1,nParams);
    bestFitParam(unfixparams) = bfp((nll == min(nll)),:);
    bestFitParam(fixparams(1,:)) = fixparams(2,:);
    st = startTheta;
    startTheta = nan(nStartVals,nParams);
    startTheta(:,unfixparams) = st;
    startTheta(:,fixparams(1,:)) = fixparams(2,:);
end
    

    function [starttheta] = genStartTheta
        switch modelname
            case 'FP'
                starttheta = [ randi(30,nStartVals,1) 1e-3+rand(nStartVals,1).*3 1e-3+rand(nStartVals,1).*3 rand(nStartVals,1)-.5];
                lb = [ 1 exp(-6) 1e-3 -10];
                ub = [100 6 10  10];
                plb = [1 exp(-6) 1e-3 -5];
                pub = [50 exp(1) 5 5];
            case 'FPheurs'
                starttheta = [ randi(45,nStartVals,1) 1e-3+rand(nStartVals,1).*3 1e-3+rand(nStartVals,1).*3  -1e-3-rand(nStartVals,1).*7];
                lb = [ 1 exp(-6) 1e-3 -50];
                ub = [100 6 10  -1e-3];
                plb = [1 exp(-6) 1e-3 -20];
                pub = [50 exp(1) 5 -1e-3];
            case 'uneqVar'
                starttheta = [3*rand(nStartVals,1) 3*rand(nStartVals,1)+1e-3 1e-3+rand(nStartVals,1).*3  rand(nStartVals,1)-.5];
                lb = [0 1e-3 1e-3 -50];
                ub = [10 10 10 50];
                plb = [0 1e-3 1e-3 -5];
                pub = [2 5 5 5];
            case 'REM'
                starttheta = [5+randi(50,nStartVals,1) rand(nStartVals,3) randi(15) 1e-3+rand(nStartVals,1).*3 rand(nStartVals,1)-.5];
                lb = [ 1 1e-3 1e-3 1e-3 1 1e-3 -10];
                ub = [50 1 1 1 50 10 10];
                plb = [1 1e-3 1e-3 1e-3 1 1e-3 -5];
                pub = [50 1 1 1 15 5 5];
        end
        
        % setting last two parameters
        switch binningfn
            case 0 % linear binning
                starttheta = [starttheta(:,1:end-2) -rand(nStartVals,1)*10, rand(nStartVals,1)*10 ];
                lb(end-1:end) = [-Inf 0];
                ub(end-1:end) = [0 Inf];
                plb(end-1:end) = [-10 0];
                pub(end-1:end) = [0 10];
            case {2,3} % log binning
                starttheta = [starttheta(:,1:end-2) rand*20 -5+rand*10 -5+rand*10 rand*10];
                lb = [lb(1:end-2) 0 -100 -100 0];
                ub = [ub(1:end-2) 100 100 100 100];
                plb = [plb(1:end-2) 0 -10 -10 0];
                pub = [pub(1:end-2) 10 10 10 10];
        end
        
        %deleting fix parameter values
        if ~isempty(fixparams)
            starttheta(:,fixparams(1,:)) = [];
            lb(fixparams(1,:)) = [];
            ub(fixparams(1,:)) = [];
            plb(fixparams(1,:)) = [];
            pub(fixparams(1,:)) = [];
        end
    end
end
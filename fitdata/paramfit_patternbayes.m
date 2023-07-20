function [bestFitParam, nLL, startTheta, Output] = paramfit_patternbayes(modelname, binningfn, nnew_part, nold_part, fixparams, nStartVals, nConf)
%
% paramfit_patternbayes(MODELNAME) uses bayesian patternsearch (one of Luigi
% Acerbi's optimization algorithms) to find the best fit parameters (and
% coinciding negative loglikelihoods (nLLs) for MODELNAME of one subjects
% data.
%
% ===== INPUT VARIABLES =====
% MODELNAME: 'FP','VP','VPheurs',or 'UVSD', or 'REM'
% BINNINGFN: 0: linear, 1: logistic, 2: log, 3: power, 4: weibull
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
nX = 300;
nS = 50;

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
for istartval = 1:nStartVals
    if ~isempty(fixparams);
        if fixparams(1,1) == 1; % if M is a fixed parameter
            obj_func = @(x) nLL_approx_vectorized(modelname, x, binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf );
        else
            obj_func = @(x) nLL_approx_vectorized_Minterp(modelname, x, binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf );
        end
    else
        obj_func = @(x) nLL_approx_vectorized_Minterp(modelname, x, binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf );
    end
    [bfp(istartval,:) ,nll(istartval), exitflag(istartval), outputt{istartval}] = bps(obj_func,startTheta(istartval,:),lb,ub,plb,pub,options);
end

nll(exitflag < 0) = Inf;

nLL = min(nll);
bestFitParam =  bfp((nll == min(nll)),:);
Output = outputt{nll == min(nll)};
if ~isempty(fixparams)
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
        
        % getting model specific parameters
        switch modelname
            case 'FP'
                %  M, sigma
                lb = [1 1e-3];
                ub = [30 6 ];
                plb = [1 1e-3 ];
                pub = [30 3 ];
                logflag = [1 1];
            case 'UVSD'
                % mu, sigma
                lb = [0 1e-3];
                ub = [10 10];
                plb = [0 1e-3];
                pub = [2 5];
                logflag = [0 1];
            case 'REM'
                %  M g ustar c m
                lb = [ 1 1e-3 1e-3 1e-3 1];
                ub = [30 1 1 1 50];
                plb = [1 1e-3 1e-3 1e-3 1];
                pub = [30 1 1 1 15];
                logflag = [1 0 0 0 0];
        end

        % setting binnfn parameters
        switch binningfn
            case {0,1} % linear or logistic mapping
                % k
                % slope y-int
                lb = [lb 1e-3];
                ub = [ub 10];
                plb = [plb 1e-3];
                pub = [pub 5];
                logflag = [logflag 0];
            case 2      % log
                % a, b
                lb = [lb 0 -100];
                ub = [ub 100 100];
                plb = [plb 0 -10];
                pub = [pub 10 10];
                logflag = [logflag 0 0];
            case 3      % power law
                % a, b, gamma
                lb = [lb 0 -100 -30];
                ub = [ub 100 100 30];
                plb = [plb 0 -10 -5];
                pub = [pub 10 10 5];
                logflag = [logflag 0 0 0];
            case 4 % weibull binning
                % scale, shift, a, b
                lb = [lb 0 0 0 -10];
                ub = [ub 10 10 10 10];
                plb = [plb 0 0 0 -3];
                pub = [pub 10 10 3 3];
                logflag = [logflag 0 0 0 0];
        end
        
        % d0, sigma_mc
        lb = [lb -30 1e-6];
        ub = [ub 30 10];
        plb = [plb -3 1e-6];
        pub = [pub 3 3];
        logflag = [logflag 0 1];

        %deleting fix parameter values
        if ~isempty(fixparams)
            lb(fixparams(1,:)) = [];
            ub(fixparams(1,:)) = [];
            plb(fixparams(1,:)) = [];
            pub(fixparams(1,:)) = [];
            logflag(fixparams(1,:)) = [];
        end
        logflag = logical(logflag);
        
        lb(logflag) = log(lb(logflag));
        ub(logflag) = log(ub(logflag));
        plb(logflag) = log(plb(logflag));
        pub(logflag) = log(pub(logflag));
        starttheta = bsxfun(@plus,bsxfun(@times,rand(nStartVals,size(lb,2)),pub-plb),plb);
 
    end
end
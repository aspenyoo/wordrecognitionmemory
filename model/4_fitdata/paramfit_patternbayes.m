function [bestFitParam, nLL, startTheta, Output] = paramfit_patternbayes(modelname, binningfn, memstrengthvar, nnew_part, nold_part, fixparams,nStartVals, nConf)
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

if nargin < 6; fixparams = []; end
if nargin < 7; nStartVals = 1; end
if nargin < 8; nConf = 20; end
nX = 30;
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
for istartval = 1:nStartVals;
    obj_func = @(x) nLL_approx_vectorized(modelname, x, binningfn, memstrengthvar, nnew_part, nold_part, fixparams, nX, nS, nConf );
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
            case {'FP','FPheurs'}
                %                 starttheta = [ randi(30,nStartVals,1) 1e-3+rand(nStartVals,1).*3 1e-3+rand(nStartVals,1).*3 rand(nStartVals,1)-.5];
                lb = [1 1e-3];
                ub = [100 6 ];
                plb = [1 1e-3 ];
                pub = [50 3 ];
            case 'UVSD'
                %                 starttheta = [3*rand(nStartVals,1) 3*rand(nStartVals,1)+1e-3 1e-3+rand(nStartVals,1).*3  rand(nStartVals,1)-.5];
                lb = [0 1e-3];
                ub = [10 10];
                plb = [0 1e-3];
                pub = [2 5];
            case 'REM'
                %                 starttheta = [5+randi(50,nStartVals,1) rand(nStartVals,3) randi(15) 1e-3+rand(nStartVals,1).*3 rand(nStartVals,1)-.5];
                lb = [ 1 1e-3 1e-3 1e-3 1];
                ub = [50 1 1 1 50];
                plb = [1 1e-3 1e-3 1e-3 1];
                pub = [50 1 1 1 15];
        end
        
        % setting binnfn parameters
        switch binningfn
            case {0,1} % linear or logistic mapping
                % k, d0
                if strcmp(modelname,'FPheurs') % FP heurs model spans (-Inf, 0], so -d0 should be be negative
                    lb = [lb 1e-3 0];
                    ub = [ub 10 100];
                    plb = [plb 1e-3 0];
                    pub = [pub 2 50];
                else
                    lb = [lb 1e-3 -50];
                    ub = [ub 10 50];
                    plb = [plb 1e-3 -10];
                    pub = [pub 5 10];
                end
            case {2,3} % log or power law mapping
                % a, b, d0
                %                 starttheta = [starttheta(:,1:end-2) rand*20 -5+rand*10 -5+rand*10 rand*10];
                lb = [lb 0 -100 -30];
                ub = [ub 100 100 30];
                plb = [plb 0 -10 -1];
                pub = [pub 10 10 1];
                if memstrengthvar == 2 % if 1/(1-p(correct))
                    % a
                    plb(end-2) = -10;
                    pub(end-2) = 0;
                end
                if binningfn == 3 % power law binning
                    % lambda
                    lb = [lb -30];
                    ub = [ub 30];
                    plb = [plb -5];
                    pub = [pub 5];
                end
            case 4 % weibull binning
                % scale, shift, d0
                lb = [lb 0 0 -10];
                ub = [ub 10 10 10];
                plb = [plb 0 0 -3];
                pub = [pub 10 10 3];
        end
        
        % setting sigma_mc parameters
%         if ~strcmp(modelname,'UVSD')
            lb = [lb 0];
            ub = [ub 10];
            plb = [plb 0];
            pub = [pub 3];
%         end

        starttheta = [bsxfun(@plus,randi(pub(1)-plb(1),nStartVals,1),plb(1)) bsxfun(@plus,bsxfun(@times,rand(nStartVals,size(lb,2)-1),pub(2:end)-plb(2:end)),plb(2:end))];
        if strcmp(modelname,'UVSD')
            starttheta(:,1) = rand(size(starttheta,1),1)*pub(1);
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
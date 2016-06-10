function [bestFitParam, nLL, startTheta, outputt] = paramfit_patternsearch_correct(modelname,nnew_part, nold_part, fixparams,nStartVals, nConf)
% TAKEN FROM PREVIOUS MODEL 2015-09-13
%
if nargin < 4; fixparams = []; end
if nargin < 5; nStartVals = 1; end
if nargin < 6; nConf = 20; end
islogbinning = 1;
nX = 30;
nS = 20;

% paramfit_fminsearch(MODELNAME) uses genetic algorithm to find the best fit parameters (and
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
% last updated July 24, 2015 by Aspen Yoo
%
% adapted from paramfit_fmincon.m

% if nargin < 5; fixparams = []; end
% if nargin < 6; nStartVals = 10; end
% if nargin < 7; iterMs = 0; end
% if isempty(iterMs); iterMs = 0; end
% if (iterMs); fixparams = [1; nan]; end % iterates fixed M from 1 to 50.
% if strcmp(modelname,'uneqVar'); fixparams = []; nStartVals = nStartVals*length(iterMs); iterMs = 0; end % no M to iterate over, so never iterate over uneqVar model.

% random number generator
rng('shuffle');

options = psoptimset('Display','iter'); %,'MaxIter',5,'MaxFunEvals',2000*4);
% options = [];

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
    obj_func = @(x) nLL_approx_vectorized(modelname, x, islogbinning, nnew_part, nold_part, fixparams, nX, nS, nConf );
%     if fixparams(1,1) == 1;
%         [bfp(istartval,:) ,nll(istartval), exitflag(istartval), outputt{istartval}] = patternsearch(obj_func,startTheta(istartval,:),A,b,Aeq,beq,lb,ub,[],options);
%     else
        [bfp(istartval,:) ,nll(istartval), exitflag(istartval), outputt{istartval}] = patternsearch(obj_func,startTheta(istartval,:),A,b,Aeq,beq,lb,ub,nonlcon,options);
%     end
end



nll(exitflag < 0) = Inf;


nLL = min(nll);
bestFitParam =  bfp((nll == min(nll)),:);
if ~isempty(fixparams);
    nParams = size(fixparams,2) + length(bestFitParam);
    unfixparams = 1:nParams;
    unfixparams(fixparams(1,:)) = [];
    bestFitParam = nan(1,nParams);
    bestFitParam(unfixparams) = bfp((nll == min(nll)),:);
    bestFitParam(fixparams(1,:)) = fixparams(2,:);
end
    

    function [starttheta] = genStartTheta
        switch modelname
            case 'FP'
                starttheta = [ randi(30,nStartVals,1) 1e-3+rand(nStartVals,1).*3 1e-3+rand(nStartVals,1).*3 rand(nStartVals,1)-.5];
                % A = [0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 0]; b = [0; 0; 0; 0];
                A = []; b = [];
                Aeq = [];
                beq = [];
                lb = [1 exp(-6) 1e-3 -5];
                ub = [50 exp(1) 5 5];
                nonlcon = @mycon;
            case 'FPheurs'
                starttheta = [ randi(45,nStartVals,1) 1e-3+rand(nStartVals,1).*3 1e-3+rand(nStartVals,1).*3  -1e-3-rand(nStartVals,1).*7];
                % A = [0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 0]; b = [0; 0; 0; 0];
                A = []; b = [];
                Aeq = [];
                beq = [];
                lb = [1 exp(-6) 1e-3 -20];
                ub = [50 exp(1) 5 -1e-3];
                nonlcon = @mycon;
        end
        
        % setting last two parameters
        if ~(islogbinning)
            starttheta = [starttheta(:,1:end-2) -rand(nStartVals,1)*10, rand(nStartVals,1)*10];
            lb(end-1:end) = [-Inf 0];
            ub(end-1:end) = [0 Inf];
            
            %             if strcmp(modelname(1),'V')
            %                 A(4,:) = []; b(4) = [];
            %             else
            %                 A(3,:) = []; b(3) = [];
            %             end
        end
        
        %deleting fix parameter values
        if ~isempty(fixparams)
            starttheta(:,fixparams(1,:)) = [];
            %             A(:, fixparams(1,:)) = [];
            %             A(fixparams(1,:),:) = [];
            %             b(fixparams(1,:)) = [];
            lb(fixparams(1,:)) = [];
            ub(fixparams(1,:)) = [];
        end
    end

    function [c,ceq] = mycon(x)
        c = [];
        ceq = round(x(1))-x(1);
    end
end
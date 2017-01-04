function [ responses, trueParam ] = simulate_resp(modelname, nSubj, fixparams)
% simulate_resp simulates the number of responses for each confidence
% value.
%
% [RESPONSES, TRUEPARAM] = simulate_resp(MODELNAME,ISLOGBINNING) simulates
% the 20 confidence RESPONSES for model MODELNAME and parameter combination
% TRUEPARAM with logistic or linear binning (ISLOGBINNING)
%
% [RESPONSES, TRUEPARAM] = simulate_resp(MODELNAME, ISLOGBINNING, NSUBJ)
% gives a NSUBJ x 20 matrix of responses for different parameter values
%
%
% RESPONSES: struct containing values for new and old words (150 total for
% each). each struct value (old, new) contains a nSubj x 20 matrix.
% TRUEPARAM: nSubj x nParams matrix of true parameter values
%
% MODELNAME should be 'FP', 'VP', 'FPheurs', 'VPheurs', or 'uneqVar'
% ISLOGBINNING = 1 for log binning and 0 for lin binning.
% NSUBJ: number of simulated subject. default: 1.
% FIXPARAMS is a 2 x k matrix, where k is the number parameters you want to
% fix. first row should be the position of parameter you wish to fix (e.g.,
% \sigma is the second position for all models) and second row should be
% the fixed value you would like. default: []; for current version, the
% second row for the first nParams-2 parameters should be the INDEX of the
% parameter vector grid you want, not the actual VALUE! (e.g., if you want
% exp(-6) for sigma, you would write one as the "value")
% SAVEDATA: 1 - yes. 0 - no (default).
%
%
% Aspen Yoo -- February 2, 2016

% random number generator
% rng('shuffle');

if nargin < 2; nSubj = 1; end
if nargin < 3; fixparams = []; end

% switch modelname
%     case 'FP21', 
        nParams = 6;
% end


% setting variable sizes
responses.new = nan(nSubj,20);
responses.old = nan(nSubj,20);
trueParam = nan(nSubj,nParams);

% getting best fit params for reference
load(['paramfit_patternbayes_' modelname '.mat'])
MU = mean(bestFitParam);
SIGMA = cov(bestFitParam);

fprintf('\n simulating %s responses for participant...', modelname)
for isimsubj = 1:nSubj
    fprintf('...%d \n',isimsubj);
    
    nnew_part = [150 zeros(1,19)]; nold_part = nnew_part; % setting to arbitrary values so while loop starts
    while  sum(abs(nnew_part - nold_part)) > 1.8*150 || sum(abs(nnew_part- nold_part)) < .2*150 || (nnew_part(1)+ nnew_part(end)) >= 0.8*150 || (nold_part(1) + nold_part(end)) >= 0.8*150;
        % if the distribution is not being cut up by the bins, or are
        % overlapping
        
        % setting theta values
        theta = mvnrnd(MU,SIGMA);
        
%         switch modelname
%             case 'FP21'; 
                modelname = 'FP';
                binningfn = 2; 
                memstrengthvar = 1;
                
                theta(1) = min([abs(round(theta(1))) 75]);
                theta(5) = 0;
                theta(2) = abs(theta(2));
                theta(6) = abs(theta(6));
%         end
theta        

        if ~isempty(fixparams); theta(fixparams(1,:)) = fixparams(2,:); end % if some parameters are fixed
        
        % simulating number of response for each confidence value
        try
            [pnew, pold] = nLL_approx_vectorized( modelname, theta, binningfn, memstrengthvar, nnew_part, nold_part);
            nnew_part = round(pnew.*150);
            nold_part = round(pold'.*150);
        end
    end
    
    
    % saving trueParam and response values
    trueParam(isimsubj,:) = theta;
    if ~isequal(size(nnew_part),[1 20])
        nnew_part = nnew_part';
        nold_part = nold_part';
    end
    responses.new(isimsubj,:) = nnew_part;
    responses.old(isimsubj,:) = nold_part;
    
end


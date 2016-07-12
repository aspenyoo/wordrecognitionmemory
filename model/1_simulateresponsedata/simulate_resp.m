function [ responses, trueParam ] = simulate_resp(modelname,islogbinning, nSubj, fixparams, savedata)
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

if nargin < 3; nSubj = 1; end
if nargin < 4; fixparams = []; end
if nargin < 5; savedata = 0; end % won't save data automatically

switch modelname
    case 'FP'; nParams = 4;
    case 'FPheurs'; nParams = 4;
    case 'VP'; nParams = 5;
    case 'VPheurs'; nParams = 5;
    case 'uneqVar'; nParams = 4;
    case 'REM'; nParams = 7;
end

% setting variable sizes
responses.new = nan(nSubj,20);
responses.old = nan(nSubj,20);
trueParam = nan(nSubj,nParams);

fprintf('\n simulating %s responses for participant...', modelname)
for isimsubj = 1:nSubj;
    fprintf('...%d \n',isimsubj);
    
    nnew_part = 150; nold_part = nnew_part; % setting to arbitrary values so while loop starts
    while  sum(abs(nnew_part - nold_part)) > 1.8*150 || sum(abs(nnew_part- nold_part)) < .2*150 || (nnew_part(1)+ nnew_part(end)) >= 0.8*150 || (nold_part(1) + nold_part(end)) >= 0.8*150;
        % if the distribution is not being cut up by the bins, or are
        % overlapping
        
        % setting theta values
        switch modelname
            case 'FP'; theta = [5+randi(45) exp(-6)+rand*3];
            case 'VP'; theta = [5+randi(45) rand*3 rand*3];
            case 'FPheurs'; theta = [randi(50) rand*3];
            case 'VPheurs'; theta = [randi(50) rand*3 rand*3];
            case 'uneqVar'; theta = [rand*3-.5 rand*3+.5];
            case 'REM'; theta = [5+randi(45) rand(1,3) randi(10)];
        end
        
        if (islogbinning)
            theta = [theta rand*3 rand-.5];
            if strcmp(modelname(end),'s');
                theta(end) = -rand*15;
            end
        else% if linear binning
            theta = [theta rand*-4 rand*4];
        end
        
        if ~isempty(fixparams); theta(fixparams(1,:)) = fixparams(2,:); end % if some parameters are fixed
        
        % simulating number of response for each confidence value
        try
        [nnew_part, nold_part] = simulate_data(modelname, theta, islogbinning);
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

% saving the data
if (savedata)
    save(['simresp_' char(modelname) '.mat'],'filename','modelname','simresp','trueParam')
end



function nLL = nLL_approx_vectorized_Minterp(modelname, theta, binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf )
% calculates the log-likelihood of a paarameter combination with a
% non-integer M value. 
% used in paramfit_patternbayes
if nargin < 6; 
    switch modelname
        case 'FP'
            logflag = [1 1];
        case 'REM'
            logflag = [1 0 0 0 0];
        case 'UVSD'
            logflag = [0 1];
    end
    switch binningfn
        case 3 % power law
            logflag = [logflag 0 0 0];
        case 4 % cumulative weibull
            logflag = [logflag 0 0 0 0];
    end
    logflag = [logflag 0 1];
    logflag = logical(logflag);
end
if nargin < 7; fixparams = []; end
if nargin < 8; nX = 30; end
if nargin < 9; nS = 100; end
if nargin < 10; nConf = 20; end

if ~isempty(fixparams); assert(fixparams(1,1) ~= 1);end
if ~isempty(fixparams) && (nargin < 6); logflag(fixparams(1,:)) = [];end

M = exp(theta(1));
lowM = floor(M);
highM = ceil(M);
logflag(1) = 0;

lownLL = nLL_approx_vectorized( modelname, [lowM theta(2:end)], binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf );
highnLL =  nLL_approx_vectorized( modelname, [highM theta(2:end)], binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf );

slopee = (highnLL-lownLL)/(highM-lowM);

nLL = slopee*(M-lowM) + lownLL;

end
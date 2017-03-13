function nLL = nLL_approx_vectorized_Minterp(modelname, theta, binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf )
% calculates the log-likelihood of a paarameter combination with a
% non-integer M value. 
% used in paramfit_patternbayes

if ~isempty(fixparams); assert(fixparams(1,1) ~= 1);end

M = exp(theta(1));
lowM = floor(M);
highM = ceil(M);
logflag(1) = 0;

lownLL = nLL_approx_vectorized( modelname, [lowM theta(2:end)], binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf );
highnLL =  nLL_approx_vectorized( modelname, [highM theta(2:end)], binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf );

slopee = (highnLL-lownLL)/(highM-lowM);

nLL = slopee*(M-lowM) + lownLL;

end
function clusterwrapper_parforsubj(modelname,binningfn,nStartVals,fixparams,optimMethod)

if nargin < 3; nStartVals = 1; end
if nargin < 4; fixparams = []; end
if nargin < 5; optimMethod = 'patternbayes'; end

parfor isubj = 1:14;
	nStartVal = max([nStartVals-countnum(modelname,isubj,fixparams(2,1)) 0]);
	fitdata_cluster(isubj,modelname,binningfn,optimMethod, fixparams,[],[],nStartVal);
end
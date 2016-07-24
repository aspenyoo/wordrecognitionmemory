function clusterwrapper_parforsubj2(modelname,binningfn,nStartVals,fixparams,optimMethod)

if nargin < 3; nStartVals = 1; end
if nargin < 4; fixparams = []; end
if nargin < 5; optimMethod = 'patternbayes'; end

parfor isubj = 1:14;
	nStartVal = max([nStartVals-countnum2(modelname,binningfn, isubj, fixparams)+countnum2(modelname,binningfn, isubj, [fixparams [6;0]]) 0]);
	fitdata_cluster(isubj,modelname,binningfn,optimMethod, fixparams,[],[],nStartVal);
end
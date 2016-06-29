%% checking to see SD of nLLs for different model
% 6/29/2016

modelname = 'FP';
nSamples = 100;

% simulate data
[resp, tp] = simulate_resp(modelname,1,1,[])

paramtry = [12 1.5 .81 .3];
nLLVec = nan(1,nSamples);
for isample = 1:nSamples
    isample
    nLLVec(isample) = nLL_approx_vectorized(modelname,paramtry,1,resp.new,resp.old);
end

hist(nLLVec,15); defaultplot
std(nLLVec)
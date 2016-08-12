function [d_new, d_old] = calculate_d(M, sigma, nS, Nnew, Nold, SNew, SOld, X)
% calculate log odds
% function used as a basis for Luigi to code up C code!
% 
% ================ OUTPUT VARIABLES ==================
% D_NEW: log odds of new trials. Nnew*nS x 1. (double)
% D_OLD: log odds of old trials. Nnew*nS x 1. (double)
%
% ================ INPUT VARIABLES ====================
% M: number of features. scalar. (integer)
% SIGMA: memory noise. scalar. (double)
% NS: number of samples of SNew and SOld. scalar. (integer)
% NNEW: number of new words. scalar. (integer)
% NOLD: number of old words. scalar. (integer)
% SNEW: new words across S simulations. Nnew*nS x M (double)
% SOLD: old words across S simulations. Nold*nS x M (double)
% X: noisy memories. Nold x M (double)
% 
% 
% Hi, Luigi. This is everything you need to calculate d for old and new
% word trials. If forloops are faster than large matrix operations, let me
% know because this can be computed over nS samples rather than being a
% large matrix. 
% Also, if generating samples from a normal distribution is quicker in C,
% then we should generate SNew and SOld and maybe X in C as well. 

J = 1/sigma.^2+1;
d_new = M/2*log(1+1/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
    sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew*nS])).^2,2)))));
d_old = M/2*log(1+1/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
    sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold*nS])).^2,2)))));
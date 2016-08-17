%% EXAMPLE2 -- Log odds for REM model

M = 40;
g = 0.8;
c = 0.8;
nS = 30;
Nnew = 150;
Nold = 149; 

ustar = 0.8; 
m = 10;

p0 = (1-ustar)^m;               % probability of x_ij= 0
pM = (1-p0)*geocdf(m-1,c);      % probability of x_ij = s_ij (match)
pQ = 1 - p0 - pM;               % probability drawn randomly (mismatch)

X = binornd(1,1-p0,[Nold M]).*(geornd(g,[Nold M])+1);
% generating new and old test words
SNew = geornd(g,[1 M Nnew*nS])+1; % new words
idx = logical(binornd(1,pQ,[Nold*nS M]) + repmat((X == 0),[nS 1])); % indices of randomly drawn features
SOld = (1-idx).*repmat(X,[nS 1]) + idx.*(geornd(g,[Nold*nS M]) + 1); % old words from noisy memories

            
matchoddsVec = (c+(1-c).*(g.*(1-g).^(1:30)))./(g.*(1-g).^(1:30));

tic; [d_new, d_old] = calculate_d_REM (M, g, c, nS, Nnew, Nold, SNew, SOld, X, matchoddsVec); t1 = toc; 
tic; [d_new2,d_old2] = calculate_d_REM_mex (M, g, c, nS, Nnew, Nold, SNew, SOld, X, matchoddsVec); t2 = toc;
            
rmse = sqrt(mean(sum((d_new(:) - d_new2(:)).^2)));

fprintf('=========================\n\t TESTING MAT VS MEX CODE:\n');
fprintf('\tRMSE: %g.\n\tMATLAB time: %.1f ms.\n\tMEX time: %.1f ms.\n\tSpeed gain: %.1f.\n', ...
    rmse, t1*1e3, t2*1e3, t1/t2);
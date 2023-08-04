

%% SIMULATE DATA / CALCULATE LL OF DATA

clear all

modelname = 'FP';
binningfn = 3;

switch modelname
    case 'FP'
        theta = [5 .1]; % M sigma
    case 'VP'
         
    case 'REM'
        
end

switch binningfn
    case 3 % power law
        theta = [theta 1 0 .1 0 .1]; % a b gamma d0 sigma_mc
    case 4 % weibull
        theta = [theta 1 0 .1 .1 0 .1]; % a b shape scale d0 sigma_mc
end

% number of old and new ratings for participant (only relevant aspect when
% simulating data is that sum(nnew_part) = sum(nold_part) = 150)
[nnew_part, nold_part] = deal([150 zeros(1,19)]); 


[ pnew, pold ] = calc_nLL_approx_vectorized( modelname, theta, binningfn, nnew_part, nold_part,[],10,10);%, logflag, nConf )

% plot pnew pold
figure;
plot(1:20,pnew,1:20,pold); 
defaultplot


%% CALCULATE LL















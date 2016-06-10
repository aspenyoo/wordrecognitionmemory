function makesimdata(modelname,islogbinning,nSimSubj)
% simulates subjects and saves it into a struct in data file 'subjdata.mat'
%
% this is so that data can be saved and retrieved an fit in a similar way
% as real data on the cluster    

if nargin < 3; nSimSubj = 20; end
    load('subjdata.mat')
    
    [ responses, trueParam ] = simulate_resp(modelname,islogbinning, nSimSubj);
    nParams = size(trueParam,2);
    nConf = 20;
    
    simdata.(modelname).nnew = [nan(14,nConf); responses.new];
    simdata.(modelname).nold = [nan(14,nConf); responses.old];
    simdata.(modelname).trueparam = [nan(14,nParams); trueParam];
    
    
    save('subjdata.mat','nNew_part','nOld_part','simdata')
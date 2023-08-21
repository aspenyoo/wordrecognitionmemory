function [nnew_part, nold_part] = loadsubjdata(isubj,modelname,nConf)
if nargin < 2; modelname = []; end
if nargin < 3; nConf = 20; end

% subjects 2-6 and 9 don't have exactly 150 trials on both old and new.

load('subjdata.mat')
if isubj <= 14
    nnew_part = nNew_part(isubj,:);
    nold_part = nOld_part(isubj,:);
else
    nnew_part = simdata.(modelname).nnew(isubj,:);
    nold_part = simdata.(modelname).nold(isubj,:);
end

if nConf ~= 20
    num = round(20/nConf);
    for i = 1:nConf-1;
        nnew_parttemp(i) = sum(nnew_part(num*(i-1)+1:num*i));
        nold_parttemp(i) = sum(nold_part(num*(i-1)+1:num*i));
    end
    nnew_parttemp(nConf) = sum(nnew_part((num*(nConf-1)+1):end));
    nold_parttemp(nConf) = sum(nold_part((num*(nConf-1)+1):end));
    nnew_part = nnew_parttemp;
    nold_part = nold_parttemp;
end
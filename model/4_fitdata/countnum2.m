function count = countnum2(testmodelname,binningfn,isubj,fixparams, truemodelname, optimMethod)
if nargin < 5; truemodelname = [testmodelname num2str(binningfn)]; end
if nargin < 6; optimMethod = 'patternbayes'; end

filename = ['paramfit_' optimMethod '_' testmodelname num2str(binningfn) '_subj' num2str(isubj) '.txt'];
if isubj > 14
    filename = ['modelrecovery_patternbayes_' testmodelname num2str(binningfn) '_' truemodelname 'subj' num2str(isubj) '.txt'];
end

% ensuring that there are no space in the txt to fuck up the reading of the file
% removetxtspaces(testmodelname,isubj,optimMethod)

filepath = 'model/4_fitdata/BPSfits/';
% try
    nFixedParams = size(fixparams,2);
    data = dlmread([filepath filename]);
    count = ones(size(data,1),1); % logicals
    for iparam = 1:nFixedParams;
        count = count & data(:,fixparams(1,iparam)) == fixparams(2,iparam);
    end
    count = sum(count);
% catch
%     count = 0;
% end


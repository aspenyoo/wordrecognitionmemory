% function [bestFitParam, nLL_est] = concatparamfit_correct

%% combine .mat files into one .mat file

clear

c = clock;

filename = 'paramfit_patternbayes';
modelname = 'FP';

nSubj = 25;


filestart = ['4_fitdata/separatefiles/' modelname];

for isubj = 1:nSubj;
    filenamee = ['paramfit_patternbayes*_subj' num2str(isubj) '_'];
    varnames = {'bestFitParam', 'nLL_est','output'};
    %     filepath = [filestart 'subj' num2str(isubj)];
    [concatvars] = concatcode(filestart,filenamee,varnames);
    
    %     output = nan(length(concatvars.outputt),1);
    %     for i = 1:length(concatvars.outputt);
    %         output(i) = concatvars.outputt(i).fsd;
    %     end
    output = concatvars.output;
    
    blah = [concatvars.bestFitParam, concatvars.nLL_est output];
    blah = sortrows(blah);
    bestFitParam = blah(:,1:end-2);
    nLL_est = blah(:,end-1);
    nLL_SD = blah(:,end);
    
    save([filename '_' modelname '_subj' num2str(isubj) '_' num2str(c(2)) num2str(c(3)) num2str(c(1)) '.mat'],'bestFitParam','nLL_est','output');
    
    minnLL = min(nLL_est);
    bestFitParamm(isubj,:) = bestFitParam(nLL_est == minnLL,:);
    nLL_estt(isubj) = minnLL;
    nLL_SDD(isubj) = output(nLL_est == minnLL);
    
end

bestFitParam = bestFitParamm;
nLL_est = nLL_estt;
nLL_SD = nLL_SDD;

save([filename '_' modelname '_' num2str(c(2)) num2str(c(3)) num2str(c(1)) '.mat'],'nLL_est','bestFitParam','nLL_SD')

%% combine .txt and .mat files into one .txt file
clear

filename = 'paramfit_patternbayes';
modelname = 'FP';

nSubj = 25;

filestart = ['4_fitdata/separatefiles/' modelname];
formatSpec = ['%2.4f \t %2.4f \t %2.4f \t %4.4f \t %5.4f \t' ...
    '%2.4f \t %2.4f \t %2.4f \t %4.4f \t %4.4f \r\n'];

for isubj = 1:nSubj;
    filenamee = ['paramfit_patternbayes_' modelname '_subj' num2str(isubj) '_'];
    varnames = {'bestFitParam', 'nLL_est','output'};
    %     filepath = [filestart 'subj' num2str(isubj)];
    [concatvars] = concatcode(filestart,filenamee,varnames);

    txtfilename = [filestart '/paramfit_patternbayes_' modelname '_subj' num2str(isubj) '.txt'];
    blah = [concatvars.bestFitParam, concatvars.nLL_est nan(size(concatvars.output,1),4) concatvars.output]';
    fileID = fopen(txtfilename,'a+');
    fprintf(fileID, formatSpec, blah); % save stuff in txt file
    fclose(fileID);
                
end


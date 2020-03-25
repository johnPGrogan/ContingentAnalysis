% ContingentDrugs

clear; close all;
load('./ContingentAnalysis.mat','doHC', 'on_off', 'cond_names', 'cg', 'vr', 'dmvr','nPD');


%% drug stuff - LDE & drug type

load('ContingentDemographics.mat','ppID','pdTable');

LDE = pdTable.LDE;


plotargs = {'pearson',0,'plot_ci',1,'text',0, 'plotline', 2, 'showzero',0}; % use spearman correlation?

figure();
for j = 1:4
    for i = 1:2
        subplot(4,2,(j-1)*2+i)
        scatterRegress(LDE, sq(nanmean(vr(:,1:nPD,j,i),1)), plotargs{:});
        xlabel('LDE')
        ylabel(cond_names{j});
        if j==1; title(on_off{i});end
    end
end
%
figure();
for j = 1:2
    for i = 1:2
        subplot(2,2,(j-1)*2+i)
        scatterRegress(LDE, dmvr(1:nPD,j,i), plotargs{:});
        xlabel('LDE')
        ylabel(cg{j});
        if j==1; title(on_off{i});end
    end
end
%
figure();
for j = 1:2
    for i = 1:1
        subplot(1,2,j)
        scatterRegress(LDE, -diff(dmvr(1:nPD,j,1:2),[],3), plotargs{:});
        xlabel('LDE')
        ylabel(['ON-OFF ' cond_names{j}]);
        if j==1; title(on_off{i});end
    end
end


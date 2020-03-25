% ContingentUPDRS

clear; close all;
load('./ContingentAnalysis.mat','doHC', 'on_off', 'cond_names', 'cg', 'vr', 'dmvr');
%% updrs correls

load('./ContingentDemographics.mat', 'ppID','pdTable');

updrsTotal = [pdTable.UPDRS_III_OFF, pdTable.UPDRS_III_ON];

nPD = size(updrsTotal,1);

%% plot updrs again residual velocity

plotargs = {'pearson',0,'plot_ci',1,'text',1, 'plotline', 2, 'showzero',0}; % use spearman correlation?

figure();
for i = 1:4
%     subplot(2,2,i)
%     plot( updrsTotal(:,1:2),sq(nanmean(vr(:,:,i,1:2),1)),'x');
    for j = 1:2
        subplot(4,2,(i-1)*2+j)
        scatterRegress( updrsTotal(:,j),sq(nanmean(vr(:,1:nPD,i,j),1)) , plotargs{:});
        hold on;
        xlabel('updrs')
        ylabel(cond_names{i});
        if i==1; title(on_off{j});end
    end
%     lsline;

end
legend(on_off(1:2),'Location','Best');

%% rew effects on resid velocity vs updrs

figure();
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j)
        scatterRegress( updrsTotal(:,j),sq(nanmean(-diff(vr(:,1:nPD,i*2-1:i*2,j),[],3),1)), plotargs{:});
        hold on;
        xlabel('updrs')
        ylabel(cg{i});
        if i==1; title(on_off{j});end
    end

end
legend(on_off(1:2),'Location','Best');


%% change in updrs (ON-OFF) vs rew effects on resid vel

updrsChange = updrsTotal(:,2) - updrsTotal(:,1);

figure();
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        scatterRegress(updrsChange, dmvr(1:nPD,i,j), plotargs{:});
        xlabel('updrs change');
        ylabel(cg{i});
        if i==1; title(on_off{j});end
    end
end

figure();
for i = 1:2
   
    subplot(1,2,i);
    scatterRegress(updrsChange, diff(dmvr(1:nPD,i,1:2),[],3), plotargs{:});
    xlabel('updrs change');
    ylabel([cg{i} ' change']);
    if i==1; title(on_off{j});end
   
end

% ContingentCorrels

load('./ContingentDemographics.mat')
%%
vars = pdTable.Properties.VariableNames(3:end);
vars = vars([2 7 10]);
nVars = length(vars);
velMean = sq(nanmean(vr(:,1:nPD,:,:)));

plotargs = {'pearson',0,'plot_ci',2,'text',1, 'plotline', 1, 'showzero',0}; % use spearman correlation?

%% mot effects

figure();
for i = 1:nVars
    for j = 1:2
        subplot(2, nVars,(j-1)*nVars+i)
        scatterRegress(repmat(pdTable.(vars{i}),1,1), velMean(:,[3],j) - velMean(:,[4],j),plotargs{:});
        hold on;
        
        if i==1, ylabel(on_off{j}); end
        if j==2, xlabel(vars{i}); else, title(vars{i}); end
        
    end
end

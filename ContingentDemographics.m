function ContingentDemographics(file)
% ContingentDemographics(file)
% Load up the demographic data, split into PD/HC, give summary statistics
% Requires ./ContingentAnalysis.mat and ./ContingentPDData.mat
%   file = path&filename of demographics file 
%           (e.g. './Contingent_PD_Demographics.xlsx')
% 
% John P Grogan, 2020.
% 
%% load
if ~exist('file','var')
    file = '../../../../OneDrive - Nexus365/Comp_Neurology_Lab/PD_ON_OFF/Contingent_PD_Demographics.xlsx';
end
% change this path
[num, ~, raw] = xlsread(file,'Demographics');

%%

varNames = raw(1,:);

ppID = raw(2:end,1); % ID numbers
nPP = sum(~cellfun(@(x) any(isnan(x)), ppID)); 


demoTable = array2table(raw(2:nPP+2, [1 3]), 'VariableNames', varNames([1 3])); % ppID and gender are not in num so take as cell

demoTable2= array2table(num(:,[1 3:(length(varNames)-1)]),'VariableNames', varNames([2 4:end]));
demoTable = horzcat(demoTable, demoTable2);

toRemove = nancat(cellfun(@(x) any(isnan(x)), demoTable.ID, 'UniformOutput', 0));

demoTable(toRemove, :) = []; % remove rows full of NaNs
ppID = demoTable.ID;


%% only keep pps with data 

load('./ContingentAnalysis.mat','bad_sub','bad_sub2');
load('./ContingentPDData.mat','ppIDs')

pdToKeep = ppIDs(~any(cellfun(@isempty, ppIDs(:,1:2)),2),1);
% pdToKeep(bad_sub) = [];

hcToKeep = ppIDs(~any(cellfun(@isempty, ppIDs(:,3)),2),3);
% hcToKeep(bad_sub2(1:length(hcToKeep))) = [];

[~, pps] = ismember([pdToKeep; hcToKeep], ppID);
demoTable = demoTable(pps, :);

% set bad_sub to NaN
bad_subs = ismember(demoTable.ID, [pdToKeep(bad_sub(1:length(pdToKeep))); hcToKeep(bad_sub2(1:length(hcToKeep)))]);
demoTable(bad_subs,2) = deal({''});
demoTable(bad_subs,3:end) = deal({NaN});

ppID = demoTable.ID;

%% do some stats

cellregexpi = @(cellArray, pattern) ~cellfun(@isempty, regexpi(cellArray, pattern));

isPD = cellregexpi(ppID, 'PD');
nPD = sum(isPD);

pdTable = demoTable(isPD,:);
hcTable = demoTable(~isPD,:);

stats(1) = summary(pdTable);
stats(2) = summary(hcTable);



%% put means & SD into table
v = pdTable.Properties.VariableNames([3 4 8:end]);

statsTable = array2table(NaN(5,length(v)));
statsTable.Properties.VariableNames = v;
statsTable.Properties.RowNames = {'PD mean', 'PD SD', 'HC mean', 'HC SD', 'p'};
for i = 1:length(v)
    statsTable.(v{i})(1) = round(nanmean(pdTable.(v{i})),2); % PD mean
    statsTable.(v{i})(2) = round(nanstd(pdTable.(v{i})),2); % PD SD
    
    statsTable.(v{i})(3) = round(nanmean(hcTable.(v{i})),2); % HC mean
    statsTable.(v{i})(4) = round(nanstd(hcTable.(v{i})),2); % HC SD
   
    [~, statsTable.(v{i})(5)] = ttest2(pdTable.(v{i}), hcTable.(v{i})); % p val from t-test
    statsTable.(v{i})(5) = round(statsTable.(v{i})(5), 4);
end


%%

save('ContingentDemographics.mat')

end
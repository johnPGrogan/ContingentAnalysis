function DrawSaccadeFigures(tHC, dv,  varNames, condNames, motNames, saveFigs)
% function DrawSaccadeFigures(tHC, dv, varNames, condNames, motNames, saveFigs)
% Draw figures for main saccade measures.
% Inputs:
%   tHC: structure with fields named by dv
%   dv: cell array of tHC field names
%   varNames: cell array of variable names for ylabel
%   condNames: cell array of condition names
%   motNames: cell array of motivational effect names
%   saveFigs: 1 = save, 0 = not
% dependencies:
%   matlib [errorBarPlot.m, sq.m]

if ~exist('saveFigs','var'), saveFigs = 0; end

col = @(x) x(:); % func to columnise variables

n = 3;

nPP = size(tHC.(dv{1}), 2);

scatterargs = {'MarkerEdgeAlpha', .4, 'Marker', '.',}; % args for scatter plots

f = figure();
cols = get(gca,'ColorOrder');  % colour order
close(f);

width = .1/(4-n); % distance between groups
sep = width/2; % width of jitter
x = [1 2 NaN 3 4]' + linspace(-width, width, n); % x points

for i = 1:length(dv)
%     subplot(3,2,i+(i>1))
    f = figure();
    
    data = sq(nanmean(tHC.(dv{i})(:,:,:,1:n))); % get measure
    hold on;
    % jitter x
    jitters = permute(repmat(x([1 2 4 5],:),1,1,30),[3,1,2]) + rand(size(data))*sep - sep/2;
    
    if i==1, n=2; else, n=3; end
    for j = 1:n
        scatter(col(jitters(:,:,j)), col(data(:,:,j)), 120,  cols(j,:), scatterargs{:});
    end

    % error bar plot
    set(gca,'ColorOrder',cols(n:-1:1,:));
    h = errorBarPlot([data(:,1:2,n:-1:1) NaN(30,1,n) data(:,3:4,n:-1:1)],'type','line','plotargs',{'LineWidth',2},'xaxisvalues',x(:,n:-1:1));

    n = 3;
    set(gca,'XTick',1:4, 'XTickLabel',condNames);
%     xlabel(xTitle);
    ylabel(varNames{i});
    xlim([.5 4.5]);
    if i==1, ylim([-80 80]); yline(0, ':k'); end
    if saveFigs; saveas(f, ['ContingentMain_' dv{i} '.svg']); end
end

% plot difference effect

% subplot(3,2,2);
% f = figure();
% data = sq(-diff(reshape(nanmean(tHC.(dv{1})),nPP,2,2,3),[],2));
% 
% hold on;
% jitters = permute(repmat(x([1 2],:),1,1,30),[3,1,2]) + rand(size(data(:,:,1:n)))*sep - sep/2;
% 
% 
% for j = 1:n
%     scatter(col(jitters(:,:,j)), col(data(:,:,j)), 120,  cols(j,:), scatterargs{:});
% end
% 
% set(gca,'ColorOrder',cols(n:-1:1,:));
% h = errorBarPlot(data(:,:,n:-1:1),'type','line','plotargs',{'LineWidth',2},'xaxisvalues',x(1:2,n:-1:1));
% 
% yline(0,':k');
% set(gca,'XTick',1:2, 'XTickLabel', motNames);
% xlabel('Motivational effect');
% ylabel(['\Delta ' varNames{1}])
% xlim([.6 2.4]);
% % ylim([-250 250])
% 
% if saveFigs; saveas(f, ['ContingentMain_Eff_' dv{1} '.svg']); end


% vel resids just HC
f = figure();
i = 1;
data = sq(nanmean(tHC.(dv{i})(:,:,:,1:n))); % get measure
hold on;
% jitter x
jitters = permute(repmat(x([1 2 4 5],:),1,1,30),[3,1,2]) + rand(size(data))*sep - sep/2;

for j = 3
    scatter(col(jitters(:,:,j)), col(data(:,:,j)), 120,  cols(j,:), scatterargs{:});
end

% error bar plot
set(gca,'ColorOrder',cols(n:-1:1,:));
h = errorBarPlot([data(:,1:2,3) NaN(30,1,1) data(:,3:4,3)],'type','line','plotargs',{'LineWidth',2},'xaxisvalues',x(:,3));

set(gca,'XTick',1:4, 'XTickLabel',condNames);
%     xlabel(xTitle);
ylabel(varNames{i});
xlim([.5 4.5]);
ylim([-200 200]);
yline(0,':k');
if saveFigs; saveas(f, ['ContingentMain_HC_' dv{i} '.svg']); end

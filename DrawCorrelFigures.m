function DrawCorrelFigures(dmvr, on_off)
% function DrawCorrelFigures(dmvr, on_off)
% Correlate velocity residual effects between conds/meds
% Inputs:
%   dmvr = matrix of mean effects of motivation on velocity residuals [pp, motivation type, group]
%   on_off = cell array of group names
% dependencies:
%   matlib [sq.m, scatterRegress.m]


    figure();
    cols = get(gca,'ColorOrder');
    
    axes = [-100 100 -100 100; -100 100 -100 100; -150 150 -150 150;];

    plotargs = {'pearson',0,'plot_ci',0,'text',2, 'plotline', 2, 'showzero',0}; % use spearman correlation?

    set(0,'DefaultAxesColorOrder',[0 0 0; 0 0 0]);
    for i = 1:3
        subplot(2,3,i)
        x = dmvr(:,1:2,i); x(all(isnan(x),2),:) = [];
        withinPlot();
    end


    % now plot differences
    subplot(2,3,4)
    x = sq(dmvr(:,1,1:2)); x(all(isnan(x),2),:) = [];
    betweenPlot({'Contingent Effect PD ON'; 'Contingent Effect PD OFF'})

    subplot(2,3,5)
    x = sq(dmvr(:,2,1:2)); x(all(isnan(x),2),:) = [];
    betweenPlot({'Guaranteed Effect PD ON'; 'Guaranteed Effect PD OFF';})

    subplot(2,3,6)
    x = [dmvr(:,1,1)-dmvr(:,1,2), dmvr(:,2,1)-dmvr(:,2,2)]; x(all(isnan(x),2),:) = [];
    betweenPlot({'Contingent Effect: PD ON - OFF'; 'Guaranteed Effect: PD ON - OFF'})
    
    % makeSubplotScalesEqual(2,3);



    set(0,'DefaultAxesColorOrder', 'factory'); % restore

    function withinPlot()
        [~, ~, ~, ~, ~, h] = scatterRegress( x(:,1), x(:,2)  , plotargs{:});
        hold on; yline(0, 'k:'); xline(0,':k');
        xlabel(['Contingent Effect: ' on_off{i}]);  ylabel(['Guaranteed Effect: ' on_off{i}]);
        set(gca,'XTick',-100:100:100,'YTick',-100:100:100);
        axis(axes(i,:))
        h(1).MarkerEdgeColor = cols(i,:);
        h(1).Marker = 'x';
        h(1).SizeData = 24;
        axis('square')
        plot(nanmean(x(:,1)), nanmean(x(:,2)),  'Color', 'k', 'Marker','.', 'MarkerSize',20);
    end
    
    function betweenPlot(labels)
        [~, ~, ~, ~, ~, h] = scatterRegress( x(:,1), x(:,2)  , plotargs{:});hold on; yline(0, 'k:'); xline(0,':k');
        xlabel(labels{1}); ylabel(labels{2});
        set(gca,'XTick',-100:100:100,'YTick',-100:100:100);
        axis([-100 100 -100 100])
        h(1).MarkerEdgeColor = 'g';
        h(1).Marker = 'x';
        h(1).SizeData = 24;
        axis('square')
        plot(nanmean(x(:,1)), nanmean(x(:,2)),  'Color', 'k', 'Marker','.', 'MarkerSize',20);
        
    end

end
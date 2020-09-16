%% PlotResidDiag
close all;
clear;
set(0,'DefaultAxesFontSize',14);

load('ContingentAnalysis.mat', 'info')
%%


inds = [8 23 1]; % pps to use for each cond
clf 

% scatter
for i = 1:3
    h(i) = plot(info.amp(:,inds(i),i), info.vel(:,inds(i),i),'x','LineWidth',1.5); 
    hold on;
end

% main seq lines
l = lsline; 
[l(:).LineWidth] = deal(1.5);
[l(:).LineStyle] = deal('--');

% show residuals
for i = 1:3
    PlotResidLines(info.amp(:,inds(i),i), info.vel(:,inds(i),i), '-', 'Color', h(i).Color);
end

% plot invis thing for legend
h(4) = plot(NaN, 'k--');
h(5) = plot(NaN, 'k-');

legend(h([4 5 1 2 3]), {'main sequence', 'residual', 'PD ON', 'PD OFF', 'HC'}, 'Location','NorthWest');

xlabel('saccade amplitude (deg)')
ylabel('saccade velocity (deg/s)')
box off

%% dummy data
rng(1)
n = 20;
x = randn(n, 3) * 1 + [10 7 13];

y = x .* [20 20 20] + [150 100 200] + randn(size(x))*30;

figure();
h = plot(x, y, 'x', 'LineWidth', 1.5);
hold on;

% main seq lines
l = lsline; 
[l(:).LineWidth] = deal(1.5);
[l(:).LineStyle] = deal('--');

% show residuals
for i = 1:3
    PlotResidLines(x(:,i), y(:,i), '-', 'Color', h(i).Color);
end

% plot invis thing for legend
h(4) = plot(NaN, 'k--');
h(5) = plot(NaN, 'k-');

legend(h([1 2 3]), {'PD ON', 'PD OFF', 'HC'}, 'Location','NorthWest');

xlabel('saccade amplitude (deg)')
ylabel('saccade velocity (deg/s)')
box off


%%
n = 20;
goodIdx = info.amp > 5 & info.amp < 14 & info.vel > 200 & info.vel < 600;
figure();
cols = get(gca,'ColorOrder');
for i = 1:30
    subplot(5,6,i)
    hold on;
    for j=1:3
        points = find(goodIdx(:,i,j), n);
        if length(points) == n
            plot(sq(info.amp(points,i,j)), sq(info.vel(points, i,j)),'x','LineWidth',1.5, 'Color', cols(j,:)); 
            hold on;
            PlotResidLines(sq(info.amp(points,i,j)), sq(info.vel(points, i,j)), '-', 'Color', cols(j,:));
        end
    end
    lsline;
end

%%


clf;
cols = get(gca,'ColorOrder');
goodIdx = info.amp > 5 & info.amp < 14 & info.vel > 200 & info.vel < 600; % good points to plot
inds = [16 16 23; 16 16 24; 27 27 3];
n = 20; % take first 20 good points

for i = 1:3
    points = find(goodIdx(:,inds(1,i),i),n);
    h(i) = plot(info.amp(points,inds(1,i),i), info.vel(points,inds(1,i),i),'.','LineWidth',1.5, 'Color', cols(i,:),'MarkerSize',10); 
    hold on;
    PlotResidLines(info.amp(points,inds(1,i),i), info.vel(points,inds(1,i),i), '-', 'Color', [cols(i,:) .5], 'LineWidth',1);

end

% main seq lines
l = lsline; 
[l(:).LineWidth] = deal(1.5);
[l(:).LineStyle] = deal('-');
for i = 1:3, l(i).Color(4) = 1; end

% plot invis thing for legend
% h(4) = plot(NaN, 'k--');
% h(5) = plot(NaN, 'k-');

% legend(h([4 5 1 2 3]), {'main sequence', 'residual', 'PD ON', 'PD OFF', 'HC'}, 'Location','SouthEast');
legend(h([1 2 3]), {'PD ON', 'PD OFF', 'HC'}, 'Location','SouthEast');
xlabel('saccade amplitude (deg)')
ylabel('saccade velocity (deg/s)')
box off
yticks(200:100:600)
function DrawVelFigures(dmvAf, dmAccel, motNames, on_off, doHC)
% function DrawVelFigures(dmvAf, dmAccel, motNames, on_off, doHC)
% Draw velocity/acceleration traces with permutation testing
% Inputs:
%   dmvAf: matrix of mean velocities over time 
%   dmAccel: matrix of mean accelerations over time
%   motNames: motivational effect names
%   on_off: group names
%   doHC: 1 = include HC, 0 = not.
% dependencies:
%   ContingentAnalysis git repo [pbar.m]
%   matlib [permutationOLS.m, errorBarPlot.m, sq.m]

if exist('doHC','var') && doHC, n=3; else n=2; end

nPP = size(dmvAf,1);

figure();
doPerm = n; % 2=PD ON and OFF sep vs 0, 3 = HC also. 0 = none.
useClust = 1;
nPerms = 5000;
c = get(gca,'ColorOrder');
c = c(n:-1:1,:); % reverse colour order, so PD ON on top
args = {'area',1,'plotIndividuals',0, 'doStats', 0, 'alpha', .4};

% plot velocity
for i = 1:2
    subplot(2,2,i*2-1);
    set(gca,'ColorOrder',c);
    
    h = errorBarPlot(sq(dmvAf(:,:,i,n:-1:1)), args{:}); % plot
    
    hold on; yline(0, 'k:');
    box off;
    ylabel([motNames{i} ' effect']);
    xlim([0 50])
    ylim([-25 25])

    % permutation testing 
    if doPerm
        y1 = diff( sq(dmvAf(:,:,i,1:2)),[],3); % PD ON vs OFF
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
        if doPerm > 1 % test each group against zero
            for j=1:doPerm
                y1 = sq(dmvAf(:,:,i,j));
                x = ones(nPP,1);
                [~,pVel(j,:,i)]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
                hold on;
                pbar(pVel(j,:,i), 'yVal', min(ylim) + (diff(ylim)/30)*j, 'plotargs', {'Color', c(1+n-j,:),'LineWidth',5});
            end
        end
    end
   if i==1
       title('Velocity (deg/s)');
       set(gca,'XTick',[]);
   else
       xlabel('% of time through saccade')
   end
   set(gca,'XTick',0:25:50,'XTickLabel',0:50:100, 'YTick', -20:20:20);

end

% now do accel
for i = 1:2
    subplot(2,2,i*2);
    set(gca,'ColorOrder',c);

    h = errorBarPlot(sq(dmAccel(:,:,i,n:-1:1)), args{:});
    hold on; yline(0, 'k:');
    box off;
    xlim([0 50])
    ylim([-5 5])


    if doPerm
        y1 = diff( sq(dmAccel(:,:,i,1:2)),[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
        if doPerm > 1 % sep bars for each cond
            for j=1:doPerm
                y1 = sq(dmAccel(:,:,i,j));
                x = ones(nPP,1);
                [~,pAccel(j,:,i)]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
                hold on;
                pbar(pAccel(j,:,i), 'yVal', min(ylim) + (diff(ylim)/30)*j, 'plotargs', {'Color', c(1+n-j,:),'LineWidth',5});
            end
        end
    end
   if i==1
       title('Acceleration (deg/s^2)');
       set(gca,'XTick',[]);
   else
       xlabel('% of time through saccade')
   end
   set(gca,'XTick',0:25:50,'XTickLabel',0:50:100, 'YTick', -20:20:20);

end

legend(fliplr([h{:,1}]), on_off, 'Location','Best')

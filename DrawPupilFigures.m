function DrawPupilFigures(pupRew, xt, on_off, motNames, doHC)
% function DrawPupilFigures(pupRew, xt, on_off, motNames, doHC)
% Draw motivational effects on pupil dilation
% Inputs:
%   pupRew: matrix with baselined pupil dilatation [nPP, timepoints, motivations, groups]
%   xt: vector of timepoints to plot
%   on_off: cell array of group labels
%   motNames: names of motivational effects
%   doHC: 1 = include HC, 0 = not
% dependencies:
%   ContingentAnalysis git repo [pbar.m, supertitle.m]
%   matlib [permutationOLS.m, errorBarPlot.m, sq.m, makeSubplotScalesEqual.m]

if exist('doHC','var') && doHC, n=3; else, n=2; end

nPP = size(pupRew,1);

figure();
c = get(gca,'ColorOrder');
c = c(n:-1:1,:); % reverse order so PD ON on top

doPerm=3; useClust=1; nPerms=5000; 

for i = 1:2 % for each effect
    
    subplot(2,1,i);
    set(gca,'ColorOrder',c);
    
    % plot
    h = errorBarPlot(sq(pupRew(:,xt,i,n:-1:1)), 'area',1,'xaxisvalues',xt*20, 'doStats', 0, 'alpha', .4);
    
    yline(0,':k');
    
    if doPerm % permutation tests
        y1 = -diff(pupRew(:,xt,i,1:2),[],4); % ON vs OFF
        [~,p1]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p1, 'yVal', min(ylim));
        
        if doPerm > 1 % test each group against zero
            for j=1:doPerm
                y1 = sq(pupRew(:,xt,i,j));
                x = ones(nPP,1);
                [~,pPup(j,:,i)]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
                hold on;
                pbar(pPup(j,:,i), 'yVal', min(ylim) + (diff(ylim)/30)*j, 'plotargs', {'Color', c(1+n-j,:),'LineWidth',5});
            end
        end

    end
    xlabel('time (ms)')
    ylabel([motNames{i} ' effect']);
    xlim([0 1400]);
    box off
end
makeSubplotScalesEqual(2,1);
SuperTitle('Pupil dilatation (a.u.)');

legend(fliplr([h{:,1}]), on_off,'Location','Best')
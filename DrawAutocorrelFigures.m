function DrawAutocorrelFigures(corDiff, corDiffMeans, groups, limits)
% function DrawAutocorrelFigures(corDiff, corDiffMeans, groups, limits)
% Draw autocorrelation figures, with contours outlining significant
% clusters
% Inputs:
%   corDiff: matrix of autocorrelation coeffs [pp, time, time, cond, group]
%   corDiffMeans: matrix of mean autocor coefficients across pps
%   groups: group names
%   limits: imagesc limits to use
% dependencies:
%   matlib [permutationOLS.m, plotn.m]

usePermOLS = 1;
n = 2; % 2 = PD only, 3 = HC too
contourColour = [0 0 0];
contourWidth = 1.5;

figure();
titlefun = @(i,j) char(groups{j}*(2-i));

% imagesc each matrix
plotn( corDiffMeans, 'plotfun',@(x) imagesc(x, limits),'fixshape',1,'titlefun',titlefun);

lab = '% time through saccade';
subplot(2,4,1); xlabel(lab); ylabel(lab); % put labels on corner
subplot(2,4,5); xlabel(lab); ylabel(lab);

% do permutation tests against zero 
for i=1:2
    for j=1:3
        y = reshape(corDiff(:,:,:,i,j),[],50*50); y(isnan(y))=0;
        if usePermOLS
            [~,p]=permutationOLS( y, [],[],[],'cluster',1,'clusterdims',[50,50],'clustermethod','mean','two_tailed',true);
        else
          [~,p]=ttest( y );
        end
        subplot(2,4,(i-1)*4 + j + 1); 
        hold on;
        contour(reshape(p,50,50)<.05,1, 'Color',contourColour,'LineWidth',contourWidth); % outline sig clusters
        set(gca, 'XTick', [1 25 50], 'XTickLabel', 0:50:100); % xticks
    end
end
%

% permutation tests ON-OFF vs zero
for i=1:2
    % ON - OFF
    y = reshape(-diff(corDiff(:,:,:,i,1:2),[],5),[],50*50); y(isnan(y)) = 0;
    if usePermOLS
        [~,p]=permutationOLS( y, [],[],[],'cluster',1,'clusterdims',[50,50],'clustermethod','mean','two_tailed',true);
    else
      [~,p]=ttest( y ); % only does paired
    end
    
    subplot(2,4,i*4-3); 
    hold on;
    contour(reshape(p,50,50)<.05,1, 'Color',contourColour,'LineWidth',contourWidth);
    p1(i,:) = p;
    set(gca, 'YTick', [1 25 50], 'YTickLabel', 0:50:100);
    set(gca, 'XTick', [1 25 50], 'XTickLabel', 0:50:100); 

end


c = colorbar('Location','North');
c.Position = [0.3754 0.5024 0.3000 0.0340]; % put colorbar in centre
cbar_handle = findobj(gcf,'tag','Colorbar');
set(cbar_handle, 'YAxisLocation','bottom')

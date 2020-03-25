% ContingentCovarCorrel

clear; close all;
load('./ContingentAnalysis.mat','doHC', 'on_off', 'tt', 'nPD', 'xlabs', 'vr', 'dmvr', 'rw', 'nConds', 'ok_g', 'cond_names', 'bad_sub','bad_sub2');

set(0, 'DefaultFigureColormap', othercolor('PRGn11'));%crameri('broc'))
%% snip raw saccade trajectories
clear A
if doHC, saccFile = './saccades_parsed_stretched_HC.mat';
else,    saccFile= './saccades_parsed_stretched.mat';
end
if exist(saccFile,'file')
    load(saccFile,'Af');
else
    load('./ContingentPDData.mat','result','resultHC');
    result = nancat(2, result, resultHC);
    for i=1:size(result,1)
      for j=1:nConds
        if ~isempty(result{i,j})
        s = [result{i,j}.presnip];
        fprintf('snipping %g:%g\n',i,j);
        % grab the first saccade amplitude > 1 degree, after the target
        [Aij, rawInfoij] =snipSaccades(s,'starttarget_t','saccadeaccepted_t','saccadeonly',1, 'minsize',33,'stretch',50 ); 
        A{i,j} = Aij;
        rawInfo{i,j} = structfun(@(x) x(:,1), rawInfoij, 'UniformOutput', 0); % only first saccade is given, so just use that
        end
      end
    end
    A=sq(nancat(A)); % trial, time, sub, drg
    A=groupMeans(A, 1, repmat(permute(tt,[1,4,2,3]),1,50) ,'dim');
    % A ( condition, time, subject, drug, trial )
    A=A-repmat(A(:,1,:,:,:),1,50); % subtract starting point
    % now do some filtering of bad trials:
    angA = wrap(2*angle(A(:,end,:,:,:))); % sacc end angle where 0 is horizontal, +/- pi is vertical.
    badA = any(isnan(A),2)      | any(abs(A)>400,2)      | any( abs(imag(A))>75 ,2)   | abs(angA)>1 ;
    % A filtered:
    Af = A + bool2nan(repmat(badA,[1,50])); % remove bad trials
    Af(:,:,bad_sub,1:2,:) = nan; % remove bad subjects
    Af(:,:,bad_sub2,3,:) = nan; % remove bad subjects

    rawInfo = transpIndex(sq(nancat(rawInfo))); % combine
    
    save(saccFile,'Af','rawInfo');
end
%%  run this only once! it changes Af.
sgnA = (Af(:,end,:,:,:)>0)*2-1; % negative if the saccade goes left
Af_sgn = Af; % keep signed version if needed
Af = Af .* repmat(sgnA, [1,50]); % flip LR
NS = size(Af,3);


%% correlations
corAf = apply([2,5], @(x)corr(x','rows','pairwise'), real(Af) ); % time-time corr ( cond, t1, sub, drug, t2 )
corAf = permute(corAf, [3,2,5,1,4]); % ( sub, t1, t2, cond, drug )

% fisher transform the correlations
corAf2 = corAf;
corAf = atanh(corAf2);

% figure();
% titlefun = @(i,j) ['pos correl: ' cond_names{i} ', ' on_off{j}];
% plotn( sq(nanmean(corAf,1)) ,'plotfun',@imagesc,'titlefun',titlefun)


%% use permutationOLS or uncorrected t-tests
usePermOLS = 1;
n = 2; % 2 = PD only, 3 = HC too
contourColour = [0 0 0];
contourWidth = 1.5;
%% cor  motivated minus not-motivated

% difference in corr/cov: dc ( sub, t1, t2, cont/guar, drug ) 
dc = diff(corAf,[],4); dc= -dc(:,:,:,[1 3],:); % differences reward minus no-reward
% rows are cont/guar, cols are on/off
dc2 = sq(nanmean(dc(:,:,:,:,1:2)));
dc2(:,:,:,3) = -diff(dc2,[],4);

conds = {'ON','OFF','ON - OFF'};
titlefun = @(i,j) char(conds{j}*(2-i));
figure();
plotn( dc2, 'plotfun',@(x) imagesc(x, [-.41 .41]),'fixshape',1,'titlefun',titlefun);
subplot(2,3,1); ylabel('Contingent effect');
subplot(2,3,4); ylabel('Guaranteed effect');

cg = {'Contingent','Guaranteed'}; hilo = {'hi mot','lo mot'};
a = zeros(size(dc));
for i=1:50, for j=1:50, if j>i, a(:,i,j,:,:)=1;end;end;end;
for i=1:2
    for j=1:2
    y = reshape(dc(:,:,:,i,j),[],50*50); y(isnan(y))=0;
    if usePermOLS
        [~,p]=permutationOLS( y, [],[],[],'cluster',1,'clusterdims',[50,50],'clustermethod','mean','two_tailed',true);
    else
      [~,p]=ttest( y );
    end
    subplot(2,3,(i-1)*3 + j); 
    hold on;
    contour(reshape(p,50,50)<.05,1, 'Color',contourColour,'LineWidth',contourWidth);
    end
end


% corr motiv/conting ON - OFF
dc = diff(corAf,[],4); dc= -dc(:,:,:,[1 3],:); % differences reward minus no-reward
dc1 = [];
dc1(:,:,:,:,:,1) = dc(:,:,:,:,1:2);
dc1(:,:,:,:,:,2) = dc(:,:,:,:,[1 3]);
dc1(:,:,:,:,:,3) = dc(:,:,:,:,2:3);
dcDiff = sq(-diff(nanmean(dc1,1),[],5));
%
% figure();
% plotn( dcDiff(:,:,:,1), 'plotfun',@(x) imagesc(x, [-.31 .31]),'fixshape',1,'titlefun',titlefun);

for i=1:2
    for j=1
        if j==1 % if paired
            y = reshape(-diff(dc1(:,:,:,i,:,j),[],5),[],50*50); y(isnan(y)) = 0;
            if usePermOLS
                [~,p]=permutationOLS( y, [],[],[],'cluster',1,'clusterdims',[50,50],'clustermethod','mean','two_tailed',true);
            else
              [~,p]=ttest( y ); % only does paired
            end
        else % if unpaired
            y = reshape(sq(dc1(:,:,:,i,:,j)),[],50*50,2); 
            [~,p] = ttest2(y(:,:,1), y(:,:,2));
        end
        subplot(2,3,i*3); 
        hold on;
        contour(reshape(p,50,50)<.05,1, 'Color',contourColour,'LineWidth',contourWidth);
        p1(i,:) = p;
    end
end
xlabel('timepoint in saccade'); ylabel('timepoint in saccade');
c = colorbar('Location','North');
c.Position = [0.3754 0.5024 0.3000 0.0340]; % put colorbar in centre
cbar_handle = findobj(gcf,'tag','Colorbar');
set(cbar_handle, 'YAxisLocation','bottom')

%% just plot on-off contingent + guaranteed
figure(); 
q=1;
for i = 1:2
    subplot(2,1,i)
    x = dcDiff(:,:,i,1);x(isnan(x)) = 0;
    imagesc(x(q:end,q:end), [-.41 .41]);
    hold on;
    x = reshape(p1(i,:),50,50)<.05;
    contour(x(q:end,q:end),1, 'Color', [0 0 0],'LineWidth',contourWidth);
    ylabel([cg{i} ' effect']);
    axis('square')
end
xlabel('timepoint in saccade')
c = colorbar('Location','East');
c.Position =  [0.6753 0.1095 0.0381 0.8167]; % put colorbar in centre
SuperTitle('Autocorrelation: ON - OFF')
xlim([0 50]); ylim([0 50]);
set(gca,'XTick',0:10:50, 'YTick', 0:10:50);

% %% cor  ON minus OFF
% % dc ( sub, t1, t2, hi/lo, cont/guar )
% dc = -diff(corAf(:,:,:,:,1:2),[],5); dc = reshape(dc, size(dc,1), 50,50,2,2);
% % rows are hi vs lo motivation, columns are contingent/guaranteed
% titlefun = @(i,j) ['ON-OFF: ' cg{j} ', ' hilo{i}];
% figure();      
% plotn( sq(nanmean(dc)), 'plotfun',@(x) imagesc(x, [-.31 .31]),'fixshape',1,'titlefun',titlefun);
% 
% a = zeros(size(dc));
% for i=1:50, for j=1:50, if j>i, a(:,i,j,:,:)=1;end;end;end;
% for i=1:2
%     for j=1:2
%     y = reshape(dc(:,:,:,i,j),[],50*50); y(isnan(y))=0;
%     if usePermOLS
%       [~,p]=permutationOLS( y, [],[],[],'cluster',1,'clusterdims',[50,50],'clustermethod','mean','two_tailed',true);
%     else
%       [~,p]=ttest( y );
%     end
%     subplot(2,2,j+(i-1)*2); 
%     hold on
%     [h,c] = contour(reshape(p,50,50)<.05,1, 'Color',contourColour,'LineWidth',contourWidth);
%     end
% end
% 
% c = colorbar('Location','East');
% c.Position = [.52 .33 .0336 .401]; % put colorbar in centre

%% covariance
covAf = apply([2,5], @(x)nancov(x'), real(Af) ); % time-time corr ( cond, t1, sub, drug, t2 )
covAf = permute(covAf, [3,2,5,1,4]);

% use signed log covar for subtractions
covAf2 = covAf;
covAf = real(log(covAf2)) .* sign(covAf2);

% figure();
% titlefun = @(i,j) ['pos covar: ' cond_names{i} ', ' on_off{j}];
% plotn( sq(nanmean(covAf,1)) ,'plotfun',@imagesc, 'titlefun',titlefun)


%% cov  motivated minus not-motivated
% difference in corr/cov: dc ( sub, t1, t2, cont/guar, drug ) 
% change covAf to corAf to look at correlation
dc = diff(covAf,[],4); dc= -dc(:,:,:,[1 3],:); % differences reward m0inus no-reward

dc2 = sq(nanmean(dc(:,:,:,:,1:2)));
dc2(:,:,:,3) = -diff(dc2,[],4);

conds = {'ON','OFF','ON - OFF'};
titlefun = @(i,j) char(conds{j}*(2-i));
figure;
plotn( dc2, 'plotfun',@(x) imagesc(x, [-1.51 1.51]),'fixshape',1,'titlefun',titlefun);
subplot(2,3,1); ylabel('Contingent effect');
subplot(2,3,4); ylabel('Guaranteed effect');

for i=1:2
    for j=1:2
    y = reshape(dc(:,:,:,i,j),[],50*50); y(isnan(y))=0;
    if usePermOLS
        [~,p]=permutationOLS( y, [],[],[],'cluster',1,'clusterdims',[50,50],'clustermethod','mean','two_tailed',true);
    else
      [~,p]=ttest( y );
    end
    subplot(2,3, (i-1)*3+j); 
    hold on;
    contour(reshape(p,50,50)<.05,1, 'Color',contourColour,'LineWidth',contourWidth);
    end
end

dc = diff(covAf,[],4); dc= -dc(:,:,:,[1 3],:); % differences reward minus no-reward
dc1 = [];
dc1(:,:,:,:,:,1) = dc(:,:,:,:,1:2);
dc1(:,:,:,:,:,2) = dc(:,:,:,:,[1 3]);
dc1(:,:,:,:,:,3) = dc(:,:,:,:,2:3);
dcDiff = sq(-diff(nanmean(dc1,1),[],5));

for i=1:2
    for j=1
        if j==1 % if paired
            y = reshape(-diff(dc1(:,:,:,i,:,j),[],5),[],50*50); y(isnan(y)) = 0;
            if usePermOLS
              [~,p]=permutationOLS( y, [],[],[],'cluster',1,'clusterdims',[50,50],'clustermethod','mean','two_tailed',true);
            else
              [~,p]=ttest( y ); % only does paired
            end
        else % if unpaired
            y = reshape(sq(dc1(:,:,:,i,:,j)),[],50*50,2); 
            [~,p] = ttest2(y(:,:,1), y(:,:,2));
        end
        subplot(2,3,i*3); 
        hold on;
        contour(reshape(p,50,50)<.05,1, 'Color',contourColour,'LineWidth',contourWidth);
    end
end

c = colorbar('Location','North');
c.Position = [0.3754 0.5024 0.3000 0.0340]; % put colorbar in centre
cbar_handle = findobj(gcf,'tag','Colorbar');
set(cbar_handle, 'YAxisLocation','bottom')

% %% cov  ON minus OFF
% % dc ( sub, t1, t2, hi/lo, cont/guar )
% dc = -diff(covAf(:,:,:,:,1:2),[],5); dc = reshape(dc, size(dc,1), 50,50,2,2);
% % rows are hi vs lo motivation, columns are contingent/guaranteed
% titlefun = @(i,j) ['ON-OFF: ' cg{j} ', ' hilo{i}];
% figure();
% plotn( sq(nanmean(dc)), 'plotfun',@imagesc,'fixshape',1,'titlefun',titlefun);
% 
% for i=1:2
%     for j=1:2
%     y = reshape(dc(:,:,:,i,j),[],50*50); y(isnan(y))=0;
%     if usePermOLS
%       [~,p]=permutationOLS( y, [],[],[],'cluster',1,'clusterdims',[50,50],'clustermethod','mean','two_tailed',true);
%     else
%       [~,p]=ttest( y );
%     end
%     subplot(2,2,j+(i-1)*2); 
%     hold on
%     [h,c] = contour(reshape(p,50,50)<.05,1, 'Color',contourColour,'LineWidth',contourWidth);
%     end
% end


%%

save('ContingentCovarCorrel.mat', 'Af');
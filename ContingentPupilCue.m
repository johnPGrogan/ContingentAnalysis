% ContingentPupilCue

clear; close all;

load('./ContingentAnalysis.mat','doHC', 'on_off', 'cg', 'tt', 'nPD', 'xlabs', 'vr', 'dmvr', 'rw', 'nConds','nPP');

%% parse pupil size for each trial
clear P Q
if doHC, pupilFile = './pupil_parsed_downsampled_HC.mat';
else,    pupilFile = './pupil_parsed_downsampled.mat';
end
if exist(pupilFile,'file')
    load(pupilFile,'P','Q');
else
    load('./ContingentPDData.mat','result','resultHC')
    result = nancat(2, result, resultHC);

    DECIMATE = false; 
    for i=1:size(result,1)
      for j=1:size(result,2)
        if ~isempty(result{i,j})
        s = [result{i,j}.presnip];
        fprintf('snipping %g:%g\n',i,j);
        % pupil after cue
        Pij=snipSaccades(s,'startcue_t','endoftrial_t','pupil',1,'speedfilter',1); 
        % now do pupil size after saccade
        Q_e_ij=snipSaccades(s,'saccadeaccepted_t','endoftrial_t', 'pupil',1,'speedfilter',0);
        Q_p_ij=snipSaccades(s,'startITI_t',       'starttarget_t','pupil',1,'speedfilter',0); 
        Q_ij  = [Q_e_ij [Q_p_ij(2:end,:); nan(1,size(Q_p_ij,2))]];
        for tr=1:size(Pij,1) % trial
          if DECIMATE
            P{i,j}(tr,:) = decimate( Pij( tr,1:2400), 20 ); % 20 ms downsampling (SLOW)
            Q{i,j}(tr,:) = decimate( Q_ij(tr,1:2400), 20 );
          else % manual downsampling
            for b=1:2400/20 % bin index
              P{i,j}(tr,b) = nanmean( Pij( tr,1+(b-1)*20:b*20) ); % 20 ms downsampling (SLOW)
              Q{i,j}(tr,b) = nanmean( Q_ij(tr,1+(b-1)*20:b*20) );
            end
          end
        end
        end
      end
    end
    % P { subject, onoff } ( trial, timepoint )
    % P { N_subj,  2     } ( 384  , 120       )
    % P after cue, Q after reward
    save (pupilFile,'P', 'Q');
end
%% plot pupil cue responses - differences by condition
p=sq(nancat(P)); % concatenate subjects and sessions
tt_p = repmat(permute(tt,[1,4,2,3]),[1,size(p,2),1,1]); % trial types for each sample
% now p ( trial, time, sub, drug) 
p=groupMeans(p,1, tt_p ,'dim'); % split by condition
% now p ( cond, time, sub, drug, trial )
p=permute(p, [5,3,2,1,4]); % ( trial, sub, time, cond, drug );
basep = repmat( p(:,:,1,:,:), [1,1,size(p,3)]); % baselie pupil (at time zero)
p=p-basep;  % subtract baseline
% different methods of baseline correction - they don't seem to make much
% difference actually. 
switch 'divide-sub' % type of baseline correction 
  case 'subtract' % divide by baseline on each trial?
  case 'divide-trial' 
    p=p./basep;
  case 'divide-sub' % divide by baseline for each subject 
    p=bsxfun(@rdivide, p, nanmean(basep(:,:,:),3 )) * 100;
end
p = p + repmat( bool2nan( any(isnan(p),3) ), [1 1 size(p,3)] ); % remove trials with any nans

% example single subject
% plotn(permute(sq(p(:,4,:,:,:)),[2,1,3,4]),'fixshape',1)
pTr = p;
p=sq(nanmean(p)); % mean across trials
p=smoothn(2,p,20); % 400 ms smooth

t = 0:400:2400;


%%
n=2;
figure();
c = get(gca,'ColorOrder');
c = c(n:-1:1,:);

subplot(2,1,1);
set(gca,'ColorOrder',c);
h = errorBarPlot( sq(-diff(p(:,:,[1 2],n:-1:1),[],3)) , 'area',1,'withinsubjecterror',0, 'doStats', 0); 
xlim([0 70]);
set(gca, 'XTick',[]);%0:20:120,'XTickLabel',t);
hold on; 
plot(xlim,[0 0],'--k','linewidth',1); hold off
% xlabel('ms after cue')
ylabel('Contigent effect')
box off

subplot(2,1,2);
set(gca,'ColorOrder',c);
errorBarPlot( sq(-diff(p(:,:,[3 4],n:-1:1),[],3)) , 'area',1,'withinsubjecterror',0, 'doStats', 0);
xlim([0 70]);
set(gca, 'XTick',0:20:120,'XTickLabel',t);
% title '10p minus 0p';
hold on; 
plot(xlim,[0 0],'--k','linewidth',1); hold off
xlabel('ms after cue')
ylabel('Guaranteed effect')
box off

legend([h{n:-1:1,1}],on_off,'Location','Best'); 


makeSubplotScalesEqual(2,1);
SuperTitle('Pupil dilatation (a.u.)');
drawnow

%% 
pupRew = p(:,:,[1 3],:) - p(:,:,[2 4],:); % rew effects
xt = 1:70;
figure();
doPerm=1; useClust=0; nPerms=5000;
for i = 1:3
    subplot(2,2,i);
    h = errorBarPlot(pupRew(:,xt,:,i), 'area',1,'xaxisvalues',xt*20, 'doStats', 0);
    yline(0,':k');
    if doPerm
        y1 = diff(pupRew(:,xt,:,i),[],3);
%         y1 = reshape(permute(pupRew(:,xt,:,i),[1,3,2]),[],length(xt));
%         x = kron([1:2]',ones(nPP,1));
        [~,p1]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
%         [~,p1] = ttest(y1);
        hold on;
        pbar(p1, 'yVal', min(ylim));
    end
    xlabel('time (ms)')
    ylabel('rew effect on pupil')
    title(on_off{i});
    xlim([0 1400]);
end
makeSubplotScalesEqual(2,2,1:3);

legend([h{:,1}], cg,'Location','Best')

%%
% region-of-interest based anova
% 30:50 is 600 to 1000 ms.
% 50:70 is 1000 to 1400 ms
roi = 50:70; 
P_roi = reshape(nanmean(p(:,roi,:,1:2),2),[size(p,1), 2,2,2]); 
[~,Tab]=anovanTable( P_roi, ...
  'varnames', {'sub','mot','cont','drg'},'display',0 , 'model', [
  1 0 0 0
  0 1 0 0
  0 0 1 0
  0 0 0 1
  0 1 1 0
  0 1 0 1
  0 0 1 1 
  0 1 1 1
  ])
pupAnova = rmanova( P_roi, {'sub','mot','cont','drg'},'categorical',4 )

% then via permutation t-test
py0  = reshape( p(:,:,:,1:2), size(p,1),size(p,2), 2,2,2 );
% [~,pval]=permutationOLS(py0);
py   = reshape( permute(py0,[1 3 4 5 2]), [],size(p,2) ); 
px0  = [ flat( bsxfun(@times, permute( [1:size(py0,1)]', [1,2]       ), ones(size(py0(:,1,:,:,:))) ) )  ...
         flat( bsxfun(@times, permute( [1:size(py0,3)]', [2,3,1]     ), ones(size(py0(:,1,:,:,:))) ) )  ...
         flat( bsxfun(@times, permute( [1:size(py0,4)]', [2,3,4,1]   ), ones(size(py0(:,1,:,:,:))) ) )  ...
         flat( bsxfun(@times, permute( [1:size(py0,5)]', [2,3,4,5,1] ), ones(size(py0(:,1,:,:,:))) ) )
       ];
px = x2fx( nanzscore(px0(:,[2:4])), fullfact([2,2,2])-1);
[ ~,pval ] = permutationOLS( py, px , eye(size(px,2)), px0(:,1) );

% subplot(1,2,1);
% hold on;
% pbar(pval(5,:),'yVal',-.2);
% subplot(1,2,2);
% hold on;
% pbar(pval(6,:),'yVal',-.2);
%% plot p values
figure();
imagep(pval, {'c','mot','cont','mot*cont','drg','mot*drg','cont*drg','m*c*d'})

%% 2-ways

PRoiCont = reshape(nanmean(p(:,roi,1:2,1:2),2),[size(p,1), 2,2]); 
pAnTab{1} = rmanova( PRoiCont, {'sub','cont','drg'} ,'categorical',3);
py0  = reshape( p(:,:,1:2,1:2), size(p,1),size(p,2), 2,2 );
% [~,pval]=permutationOLS(py0);
py   = reshape( permute(py0,[1 3 4 2]), [],size(p,2) ); 
px0  = [ flat( bsxfun(@times, permute( [1:size(py0,1)]', [1,2]       ), ones(size(py0(:,1,:,:,:))) ) )  ...
         flat( bsxfun(@times, permute( [1:size(py0,3)]', [2,3,1]     ), ones(size(py0(:,1,:,:,:))) ) )  ...
         flat( bsxfun(@times, permute( [1:size(py0,4)]', [2,3,4,1]   ), ones(size(py0(:,1,:,:,:))) ) )  ...
       ];
px = x2fx( nanzscore(px0(:,[2:3])), fullfact([2,2])-1);
[ ~,pvalCont ] = permutationOLS( py, px , eye(size(px,2)), px0(:,1) );
figure();
subplot(2,1,1)
imagep(pvalCont, {'c','cont','drg','cont*drg'})


PRoiMot = reshape(nanmean(p(:,roi,3:4,1:2),2),[size(p,1), 2,2]); 
pAnTab{2} = rmanova( PRoiMot, {'sub','mot','drg'} ,'categorical',3);

py0  = reshape( p(:,:,3:4,1:2), size(p,1),size(p,2), 2,2 );
% [~,pval]=permutationOLS(py0);
py   = reshape( permute(py0,[1 3 4 2]), [],size(p,2) ); 
px0  = [ flat( bsxfun(@times, permute( [1:size(py0,1)]', [1,2]       ), ones(size(py0(:,1,:,:,:))) ) )  ...
         flat( bsxfun(@times, permute( [1:size(py0,3)]', [2,3,1]     ), ones(size(py0(:,1,:,:,:))) ) )  ...
         flat( bsxfun(@times, permute( [1:size(py0,4)]', [2,3,4,1]   ), ones(size(py0(:,1,:,:,:))) ) )  ...
       ];
px = x2fx( nanzscore(px0(:,[2:3])), fullfact([2,2])-1);
[ ~,pvalMot ] = permutationOLS( py, px , eye(size(px,2)), px0(:,1) );
subplot(2,1,2)
imagep(pvalMot, {'c','mot','drg','mot*drg'})


%% does pupil dilation differ between trials
% % get max change in pupil size per trial - 1200:1400
pupDiff = sq(nanmean(pTr(:,:,60:70,:,:), [3]));
pupDiffMean = sq(nanmean(pupDiff));
% 
% figure();
% errorBarPlot([pupDiffMean(:,1:2,:), NaN(nPP,1,size(pupDiffMean,3)), pupDiffMean(:,3:4,:)], 'plotargs',{'LineWidth',2.5});
% ylabel('Pupil Change (a.u.)')
% set(gca,'xtick',1:5,'xticklabel',xlabs);
% xlim([.5 5.5])
% box off

% %%
% n=3;
% figure();
% errorBarPlot(pupDiffMean(:,[1 3],1:n) - pupDiffMean(:,[2 4],1:n), 'plotargs',{'LineWidth',2.5});
% ylabel('Rew Pupil Change (a.u.)')
% set(gca,'xtick',1:2,'xticklabel',cg)
% xlim([.5 2.5])
% box off
% 
% 
% pdmstats = rmanova(reshape(pupDiffMean(:,:,1:n), nPP,2,2,n), {'pp','rew','cont','drug'},'categorical',4);



%% get the diffs of cont/mot
pupDiffConds = [pupDiffMean(:,1,:) - pupDiffMean(:,2,:), pupDiffMean(:,3,:) - pupDiffMean(:,4,:)];
plotargs = {'pearson',0,'plot_ci',2,'text',2, 'plotline', 2, 'showzero',0}; % use spearman correlation?

cols = get(gca,'ColorOrder');
% plot those
figure();
for i = 1:2
    for j = 1:3
        subplot(2,3,(i-1)*3+j)
    %     for j = 1:2
        [~,~,~,~,~,h] = scatterRegress(sq(pupDiffConds(:,i,j)), sq(dmvr(:,i,j)), plotargs{:});
        h(1).CData = cols(j,:);
        hold on; yline(0,':k'); xline(0, ':k');
    %     lsline
        if i==2; xlabel('Pupil change');end
        if j==1; ylabel([cg{i} ' velocity effect']);end
        if i==1; title(on_off{j});end
    end
end
makeSubplotScalesEqual(2,3);



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%

% after reward

% %% pupil after reward
% p=sq(nancat(Q));
% p=groupMeans(p,1, repmat(permute(tt + rw,[1,4,2,3]),[1,size(p,2),1,1]),'dim'); % split by condition
% % 1:4 is no reward, 5 is contingent rew, 
% % 6 is rand rew, 7 is guaranteed reward.
% p=permute(p, [5,3,2,1,4]); % ( trial, sub, time, cond, drug );
% base = repmat( p(:,:,1,:,:), [1,1,size(p,3)]);
% p=p-base; % subtract baseline
% % remove trials with more than 10 nan samples
% p=p + repmat( bool2nan( sum(isnan(p),3)>10 ), [1 1 size(p,3)] ); 
% p(:,:,:,3,:) = []; % get rid of unrewarded "guaranteed reward" trials 
% % (there are only about 100 such trials across everyone)
% % now we have 1:3=rewarded cont/rand/guar, 4:6=unrewarded.
% p=sq(nanmean(p)); % mean across trials
% 
% %%
% figure();
% colourMap(jet(6))
% p=smoothn(2,p,20); % 400 ms smoothing
% for i=1:nConds
%   subplot(1,nConds,i)
%   h = errorBarPlot(p(:,:,1:3,i)-p(:,:,4:6,i),'area',1,'dostats',0,'withinsubjecterror',0);
%   %plot(sq(nanmean(p(:,:,1:3,i) - p(:,:,4:6,i) )),'linewidth',2)
%   title(on_off{i}); plotZeroLines;
%   xlim([0 120]);set(gca, 'XTick',0:40:120,'XTickLabel',t(1:2:end));
%   xlabel('time after reward (ms)')
% end
% legend(h(:,1),{'cont','rand','guar'}); ylabel 'reward effect after outcome'
% % tmp( sub, cond, outcome, drug )
% tmp = reshape(nanmean(p(:,90:end,:,1:2),2), [size(p,1),3,2,2]);
% rmanova( tmp, {'s','cond','rew','drg'},'categorical',4 )
% makeSubplotScalesEqual(1,nConds);
% 
% %%
% figure();
% for i = 1:3
%     subplot(3,3,i*3-2)
%     h = errorBarPlot(-diff(p(:,:,i,1:2)-p(:,:,i+3,1:2),[],4),'area',1,'dostats',0,'withinsubjecterror',0);
%     hold on; plot(xlim,[0 0],':k','linewidth',2); hold off
%     if i==1; title('ON - OFF'); end
%     if i==3; xlabel('time step post reward');end
%     ylabel(conds{i})
%     
%     subplot(3,3,i*3-1)
%     h = errorBarPlot(-diff(p(:,:,i,[1 3])-p(:,:,i+3,[1 3]),[],4),'area',1,'dostats',0,'withinsubjecterror',0);
%     hold on; plot(xlim,[0 0],':k','linewidth',2); hold off
%     if i==1; title('ON - HC'); end
%     if i==3; xlabel('time step post reward');end
% %     ylabel(conds{i})
%     
%     subplot(3,3,i*3)
%     h = errorBarPlot(-diff(p(:,:,i,2:3)-p(:,:,i+3,2:3),[],4),'area',1,'dostats',0,'withinsubjecterror',0);
%     hold on; plot(xlim,[0 0],':k','linewidth',2); hold off
%     if i==1; title('OFF - HC'); end
%     if i==3; xlabel('time step post reward');end
% %     ylabel(conds{i})
% end
% makeSubplotScalesEqual(3,3)
% %% raw traces of above
% figure();
% for i = 1:nConds
%     subplot(1,nConds,i)
%     errorBarPlot(p(:,:,:,i),'area',1,'dostats',0); title(on_off{i})
% end
% makeSubplotScalesEqual(1,nConds)

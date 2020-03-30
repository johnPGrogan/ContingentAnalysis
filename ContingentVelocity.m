% ContingentVelocity

clear; close all;
load('./ContingentAnalysis.mat','doHC', 'on_off','cg', 'tt','nPD', 'xlabs', 'nConds', 'vr','nPP','rtCutoffs','dpp');
if doHC, saccFile = './saccades_parsed_stretched_HC.mat';
else,    saccFile= './saccades_parsed_stretched.mat';
end
load(saccFile', 'Af','rawInfo');

%% remove bad trials

Af2 = Af;

% already has minsize of 33pix ~ 1 deg
MAX_AMP           = 400;    % the target is 300 px
EXCLUDE_BAD_VEL   = true;   % whether to discard velocity < 100 deg/s 
EXCLUDE_SUB_MIN_TRIALS = 10;   % put zero to include everyone (ok for mixed model)

ok = rawInfo.sAmpl <= MAX_AMP ...
    & rawInfo.sRT >= rtCutoffs(1) ...     
    & rawInfo.sRT <= rtCutoffs(2) ...
    & rawInfo.sBlink==0 ...
    ;
ok2 = rawInfo.sAmpl <= MAX_AMP & rawInfo.sRT >= rtCutoffs(1) & rawInfo.sRT <= rtCutoffs(2) & rawInfo.sBlink==0;
ok(:,3,:) = ok2(:,3,:);
if EXCLUDE_BAD_VEL
  ok = ok & rawInfo.sSpd >= (2.5) & rawInfo.sSpd <= (80);
end

% split by tt
ok_1 = permute(groupMeans(ok, 3, permute(tt,[2,3,1]), 'dim'), [1,3,2,4]);
bad_sub = sq(any(nansum(ok_1,4) < EXCLUDE_SUB_MIN_TRIALS, 2));
bad_sub(:,1:2) = repmat(any(bad_sub(:,1:2),2),[1 2]);

Af(repmat(permute(ok_1==0,[2,5,1,3,4]), 1, 50)) = NaN;
Af(repmat(permute(bad_sub==1, [3,4,1,2,5]), [4,50,1,1,120])) = NaN;

Af = Af * dpp; % convert from pix to deg

%% look at average velocity 

vAf = permute(real(diff(Af,[],2)), [3,2,1,4,5]); % take horizontal velocity
vAf = vAf .* 1000; % from deg/ms -> deg/sec
vAf = vAf .* repmat(sign(vAf(:,10,:,:,:,:)),[1 49 1 1 1 1]); % flip leftwards

vAf(repmat(any(abs(vAf)>(40*dpp*1000),2),[1 49 1 1 1])) = NaN; % remove outlying vel due to blinks

% smooth
vAf = movmean(vAf, 3, 2); % smooth by 3 time points

vAf = reshape(vAf, nPP,49,2,2,3,120); % split into conds

mvAf = nanmean(vAf, 6); % mean across trials

dmvAf = sq(-diff(mvAf,[],3)); % rew effects



% %% plot each cond
% n=2;
% figure();
% for i = 1:2
%     for j = 1:2
%         subplot(2,2,(i-1)*2+j);
%         errorBarPlot(sq(mvAf(:,:,i,j,:)),'area',1);
%         xlabel('timepoint in saccade');
%         ylabel('velocity')
%         title(xlabs{ (i-1)*2+j });
%     end
% end

% %% normalise amplitude
% 
% saccAmpl = (Af(:,end,:,:,:));
% saccAmpl = repmat(reshape(permute(saccAmpl,[3,2,1,4,5]), [nPP,1,2,2,3,120]),[1,49]);
% 
% % normAf = Af ./ saccAmpl;
% 
% % vAfNorm = permute(real(diff(normAf,[],2)), [3,2,1,4,5]);
% vAfNorm = vAf ./ abs(real(saccAmpl)); % normalise
% 
% % smooth
% vAfNorm = movmean(vAfNorm, 3, 2); % smooth by 3 time points
% 
% vAfNorm = reshape(vAfNorm, nPP,49,2,2,3,120);
% mvAfNorm = nanmean(vAfNorm,6);
% dmvAfNorm = sq(-diff(mvAfNorm,[],3));
% n = 2;
% figure();
% doPerm = 1;
% useClust = 1;
% nPerms = 5000;
% for i = 1:2
%     subplot(2,1,i);
%     
%     h = errorBarPlot(sq(dmvAfNorm(:,:,i,1:n)), 'area',1);
%     hold on; yline(0, 'k:');
%     box off;
%     xlabel('timepoint within saccade')
%     ylabel(cg{i});
% 
%     if doPerm
%         y1 = diff( sq(dmvAfNorm(:,:,i,1:2)),[],3);
%         x = ones(nPP,1);
%         [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
%         hold on;
%         pbar(p, 'yVal', min(ylim));
%     end
%     
% end
% 
% makeSubplotScalesEqual(2,1);
% SuperTitle('reward effects on normalised velocity over time');
% legend(h(:,1), on_off, 'Location', 'Best');

% %% look at indiv vel traces
% 
% figure();
% for i = 1:3
%     for j = 1:2
%         subplot(2,3,(j-1)*3+i)
%         plot(sq(vAf(1,:,1,j,i,:)),'-');
%         
%         box off;
%         title(on_off{i});
%         
%         if i==1; ylabel(cg{j});end
%         if j==2; xlabel('timepoint in saccade'); end
%     end 
% end
% 


%% acceleration
vAf2 = permute(real(diff(Af,[],2)), [3,2,1,4,5]); % horizontal acceleration
vAf2 = vAf2 .* 1000; % from deg/ms -> deg/sec
vAf2 = vAf2 .* repmat(sign(vAf2(:,10,:,:,:,:)),[1 49 1 1 1 1]); % flip those that went left

vAf2(repmat(any(abs(vAf2)>(40*dpp*1000),2),[1 49 1 1 1])) = NaN; % remove outlying vel due to blinks

vAf2 = movmean(vAf2,3,2); % smooth
vAccel = (diff(vAf2,[],2));


% smooth
accelSmooth = movmean(vAccel, 5, 2); % smooth by 3 time points

accelSmooth = reshape(accelSmooth, nPP,48,2,2,3,120); % split into conds

mAccel = nanmean(accelSmooth, 6); % mean across trials

dmAccel = sq(-diff(mAccel,[],3)); % rew effects
% %% plot each cond
% n=2;
% figure();
% for i = 1:2
%     for j = 1:2
%         subplot(2,2,(i-1)*2+j);
%         errorBarPlot(sq(mAccel(:,:,i,j,:)),'area',1);
%         xlabel('timepoint in saccade');
%         ylabel('velocity')
%         title(xlabs{ (i-1)*2+j });
%     end
% end

% %% normalise amplitude
% 
% saccAmpl = (Af(:,end,:,:,:));
% normAf = Af ./ saccAmpl;
% 
% vAfNorm = permute(abs(diff(normAf,[],2)), [3,2,1,4,5]);
% 
% vAfNorm(repmat(any(abs(vAfNorm)>0.4,2),[1 49 1 1 1])) = NaN; % remove outlying vel due to blinks
% 
% vAfNorm = movmean(vAfNorm,3,2); % smooth
% vAccelNorm = (diff(vAfNorm,[],2));
% 
% % smooth
% vAccelNorm = movmean(vAccelNorm, 5, 2); % smooth by 3 time points
% 
% vAccelNorm = reshape(vAccelNorm, nPP,48,2,2,3,120);
% mvAccelNorm = nanmean(vAccelNorm,6);
% dmvAccelNorm = sq(-diff(mvAccelNorm,[],3));
% n = 2;
% figure();
% doPerm = 1;
% useClust = 1;
% nPerms = 5000;
% for i = 1:2
%     subplot(2,1,i);
%     
%     h = errorBarPlot(sq(dmvAccelNorm(:,:,i,1:n)), 'area',1);
%     hold on; yline(0, 'k:');
%     box off;
%     xlabel('timepoint within saccade')
%     ylabel(cg{i});
% 
%     if doPerm
%         y1 = diff( sq(dmvAccelNorm(:,:,i,1:2)),[],3);
%         x = ones(nPP,1);
%         [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
%         hold on;
%         pbar(p, 'yVal', min(ylim));
%     end
%     
% end

%% vel & accel on one fig

n = 2;
figure();
doPerm = 2; % 2=PD ON and OFF sep vs 0.
useClust = 1;
nPerms = 5000;
c = get(gca,'ColorOrder');
c = c(n:-1:1,:);
for i = 1:2
    subplot(2,2,i*2-1);
    set(gca,'ColorOrder',c);
    
    h = errorBarPlot(sq(dmvAf(:,:,i,n:-1:1)), 'area',1,'plotIndividuals',1, 'doStats', 0);
    hold on; yline(0, 'k:');
    box off;
    ylabel([cg{i} ' effect']);
    xlim([0 50])
    ylim([-25 25])

    if doPerm
        y1 = diff( sq(dmvAf(:,:,i,1:2)),[],3); % PD ON vs OFF
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
        if doPerm > 1 % sep bars for each cond vs zero
            for j=1:doPerm
                y1 = sq(dmvAf(:,:,i,j));
                x = ones(nPP,1);
                [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
                hold on;
                pbar(p, 'yVal', min(ylim) + (diff(ylim)/30)*j, 'plotargs', {'Color', c(1+n-j,:),'LineWidth',5});
            end
        end
    end
   if i==1
       title('Velocity (deg/s)');
       set(gca,'XTick',[]);
   else
       xlabel('timepoint within saccade')
   end
   set(gca,'XTick',0:25:50, 'YTick', -20:20:20);

end

% now do accel
for i = 1:2
    subplot(2,2,i*2);
    
    h = errorBarPlot(sq(dmAccel(:,:,i,1:n)), 'area',1,'plotIndividuals',1, 'doStats', 0);
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
                [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
                hold on;
                pbar(p, 'yVal', min(ylim) + (diff(ylim)/30)*j, 'plotargs', {'Color', c(1+n-j,:),'LineWidth',5});
            end
        end
    end
   if i==1
       title('Acceleration (deg/s^2)');
       set(gca,'XTick',[]);
   else
       xlabel('timepoint within saccade')

   end
   set(gca,'XTick',0:25:50,'YTick',-5:5:5);

end

legend([h{:,1}], on_off, 'Location','Best')

%%

save('ContingentVelocity.mat')
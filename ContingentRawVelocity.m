% ContingentRawVelocity

clear; %close all;
load('ContingentAnalysis.mat', 'doHC', 'cg','on_off','xlabs','bad_sub','bad_sub2','nConds','tt','nPP','rtCutoffs','dpp');
%% velocities
% snip raw saccade trajectories
clear S
if doHC, saccFile = './saccades_raw_HC.mat';
else,    saccFile = './saccades_raw.mat';
end
if exist(saccFile,'file')
    load(saccFile,'Sf','rawInfo');
else
    load('./ContingentPDData.mat','result','resultHC');
    result = nancat(2, result, resultHC);
    for i=1:size(result,1)
      for j=1:nConds
        if ~isempty(result{i,j})
        s = [result{i,j}.presnip];
        fprintf('snipping %g:%g\n',i,j);
        % grab the first saccade amplitude > 1 degree, after the target
        [Sij, rawInfoij] =snipSaccades(s,'starttarget_t','saccadeaccepted_t','saccadeonly',1, 'minsize',33); 
        S{i,j} = Sij;
        rawInfo{i,j} = structfun(@(x) x(:,1), rawInfoij, 'UniformOutput', 0); % only first saccade is given, so just use that
        % get end angle
        s = alignRight(Sij);
        ang{i,j} = angle(s(:,end) - Sij(:,1));
        end
      end
    end
    save(saccFile,'S');
    S1=sq(nancat(S)); % trial, time, sub, drg
    S2=groupMeans(S1, 1, repmat(permute(tt,[1,4,2,3]),1,size(S1,2)) ,'dim');
    % A ( condition, time, subject, drug, trial )
    S3=S2-S2(:,1,:,:,:); % subtract starting point
    
    % need to nan after subtraction
    S3(isnan(real(S3))) = complex(NaN,NaN);

    % get the final angle
    ang1 = sq(nancat(ang));
    ang2 = permute(groupMeans(ang1, 1, tt,'dim'),[1,5,2,3,4]);
    
    
    % now do some filtering of bad trials:
    angS = wrap(2*ang2); % sacc end angle where 0 is horizontal, +/- pi is vertical.
    badS = all(isnan(S3),2)      | any(abs(S3)>400,2)      | any( abs(imag(S3))>75 ,2)   | abs(angS)>1 ;
    disp(nanmean(badS,'all'))
    % A filtered:
    Sf = S3 + bool2nan(repmat(badS,[1,size(S3,2)])); % remove bad trials
    Sf(:,:,bad_sub,1:2,:) = nan; % remove bad subjects
    Sf(:,:,bad_sub2,3,:) = nan; % remove bad subjects

    rawInfo = transpIndex(sq(nancat(rawInfo))); % combine

    save(saccFile,'Sf','rawInfo','-append');
end

%% remove bad trials

Sf2 = Sf;
% Af(permute(ok_g==0,[1,5,2,3,4])) = NaN;

% minsize = 33 pix ~ 1 deg
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
  ok = ok & rawInfo.sSpd >= (2.5)   & rawInfo.sSpd <= (80);
end

% split by trialtype
ok_1 = permute(groupMeans(ok, 3, permute(tt,[2,3,1]), 'dim'), [1,3,2,4]);
bad_sub = sq(any(nansum(ok_1,4) < EXCLUDE_SUB_MIN_TRIALS, 2));
bad_sub(:,1:2) = repmat(any(bad_sub(:,1:2),2),[1 2]);


% bad_sub = any( any ( (sum(~isnan(vr(:,:,:,1:2)))) < EXCLUDE_SUB_MIN_TRIALS , 3), 4 ); 
% bad_sub2 = any( any ( (sum(~isnan(vr(:,:,:,3)))) < EXCLUDE_SUB_MIN_TRIALS , 3), 4 ); 

Sf(repmat(permute(ok_1==0,[2,5,1,3,4]), 1, 50)) = complex(NaN,NaN);
Sf(repmat(permute(bad_sub==1, [3,4,1,2,5]), [4,size(Sf,2),1,1,120])) = complex(NaN,NaN);

Sf = Sf * dpp; % convert from pixels -> degrees

%% get velocities from raw traces

stretchPoints = 50;

raw_vels = real(diff(Sf,[],2)) .* 1000; % horizontal only - convert from deg/ms to deg/sec
raw_vels = raw_vels .* repmat(sign(raw_vels(:,10,:,:,:,:)),[1 size(raw_vels,2) 1 1 1 1]); % flip leftwards

smooth_vels = movmean(raw_vels, 3, 2);

% get acceleration
raw_accel = diff(raw_vels, [], 2);
smooth_accel = movmean(raw_accel, 5,2);

int_vels = NaN(4,50,nPP,3,120);
int_accel = NaN(4,50,nPP,3,120);
% interpolate to fit on same graph
for i=1:4
    for j = 1:nPP
        for k=1:3
            for l = 1:120
                npnts = max([0 find(~isnan(sq(smooth_vels(i,:,j,k,l))), 1,'last')]);
                if npnts > 10
                    int_vels(i,:,j,k,l) = interp1(sq(smooth_vels(i,1:npnts,j,k,l)), linspace(1,npnts,stretchPoints)');
                    int_accel(i,:,j,k,l) = interp1(sq(smooth_accel(i,1:npnts-1,j,k,l)), linspace(1,npnts,stretchPoints)');
                end
            end
        end
    end
end
%%
% 
mean_vels = (permute(nanmean(int_vels,5),[3,2,1,4]));
% 
% figure();
% for i = 1:4
%     subplot(2,2,i);
%     errorBarPlot(sq(mean_vels(:,:,i,:)),'area',1);
%     xlabel('time in movement');
%     ylabel('velocity');
%     title(xlabs{i})
% end


%%
dmvels = mean_vels(:,:,[1 3],:) - mean_vels(:,:,[2 4],:);

cg = {'Contingent', 'Guaranteed'};
n = 2;
figure();
doPerm = 1;
useClust = 1;
nPerms = 5000;
for i = 1:2
    subplot(2,1,i);
    
    h = errorBarPlot(sq(dmvels(:,:,i,1:n)), 'area',1, 'doStats', 0);
    hold on; yline(0,'k:');
    box off;
    xlabel('time')
    ylabel(cg{i});

    if doPerm
        y1 = diff( sq(dmvels(:,:,i,1:2)),[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    
end

makeSubplotScalesEqual(2,1)
SuperTitle('reward effects on raw velocity over time');
legend([h{:,1}], on_off,'Location','Best');

%% accel
% 
mean_accel = (permute(nanmean(int_accel,5),[3,2,1,4]));
% 
% figure();
% for i = 1:4
%     subplot(2,2,i);
%     errorBarPlot(sq(mean_accel(:,:,i,:)),'area',1);
%     xlabel('time in movement');
%     ylabel('acceleration');
%     title(xlabs{i})
% end

%%
dmAccel= mean_accel(:,:,[1 3],:) - mean_accel(:,:,[2 4],:);

% last few are huge, set to NaN
dmAccel(:,end-2:end,:,:) = NaN;

cg = {'Contingent', 'Guaranteed'};
n = 2;
figure();
doPerm = 1;
useClust = 1;
nPerms = 5000;
for i = 1:2
    subplot(2,1,i);
    
    h = errorBarPlot(sq(dmAccel(:,:,i,1:n)), 'area',1, 'doStats', 0);
    hold on; yline(0,'k:');
    box off;
    xlabel('time')
    ylabel(cg{i});

    if doPerm
        y1 = diff( sq(dmAccel(:,:,i,1:2)),[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    
end

makeSubplotScalesEqual(2,1)
SuperTitle('reward effects on raw acceleration over time');
legend([h{:,1}], on_off,'Location','Best');
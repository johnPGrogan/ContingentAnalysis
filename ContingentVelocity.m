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


%% vel & accel on one fig

DrawVelFigures(dmvAf, dmAccel, cg, on_off, doHC)

% source data
readme = "Data to create Figure 3. Please download the ContingentAnalysis GitHub repo and Matlib repo (links available in paper), and then run: DrawVelFigures(dmvAf, dmAccel, cg, on_off, 1);";
save('Figure3SourceData.mat', 'dmvAf', 'dmAccel','cg','on_off','doHC','readme');

%% plot indiv
figure();
cols = get(gca,'ColorOrder'); %cols = cols([3 1 2], :);
alpha = .5; 
n=3;
% vel
for i = 1:2
    subplot(2,2,i*2-1)
    
    for j=n:-1:1
        hold on
        plot(dmvAf(:,:,i,j)', '-', 'Color', [cols(j,:) alpha]);
    end

    if i==1
       title('Velocity (deg/s)');
       set(gca,'XTick',[]);
    else
       xlabel('% of time through saccade')
    end
    set(gca,'XTick',0:25:50,'XTickLabel',0:50:100);
    xlim([0 50]);
    ylabel([cg{i} ' effect']);

end

% accel
h = [];
for i = 1:2
    subplot(2,2,i*2)
    
    for j=n:-1:1
        hold on
        h(:,j) = plot(dmAccel(:,:,i,j)', '-', 'Color', [cols(j,:) alpha]);
    end

    if i==1
       title('Acceleration (deg/s^2)');
       set(gca,'XTick',[]);
    else
       xlabel('% of time through saccade')
    end
    set(gca,'XTick',0:25:50,'XTickLabel',0:50:100);

    xlim([0 50])
end

legend(h(1,:), fliplr(on_off), 'Location', [0.3220 0.3672 0.1482 0.1548]);
%%

save('ContingentVelocity.mat')
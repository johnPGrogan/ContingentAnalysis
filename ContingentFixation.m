% ContingentFixation
% requires edft package for fourier transforms

clear all; close all;
load('./ContingentAnalysis.mat','doHC', 'on_off','cg', 'tt','nPD', 'xlabs', 'nConds', 'vr','nPP','cond_names');
%% snip raw saccade trajectories
clear F
if doHC, fixFile = './fixation_parsed_HC.mat';
else,    fixFile = './fixation_parsed.mat';
end
load('./ContingentAnalysis.mat','nConds','nPP','tt')
if exist(fixFile,'file')
    load(fixFile); % load if already run
else % otherwise, snip and process (takes a while)
    load('./ContingentPDData.mat','result','resultHC')
    result = nancat(2, result, resultHC);
    for i=1:size(result,1)
      for j=1:nConds
        if ~isempty(result{i,j})
        s = [result{i,j}.presnip];
        fprintf('snipping %g:%g\n',i,j);
        % grab the fixation period between cue and target
        [F{i,j}, FInfo{i,j}] =snipSaccades(s,'startcue_t','starttarget_t'); 
        F{i,j} = F{i,j}(:,1:1400);
        
        fInf(i,j).amp   = nancat(1, FInfo{i,j}.sAmpl);
        fInf(i,j).blink = nancat(1, FInfo{i,j}.sBlink);
        fInf(i,j).vel = nancat(1, FInfo{i,j}.sSpd);
        fInf(i,j).rt = nancat(1, FInfo{i,j}.sRT);
        end
      end
    end
    
    F1=sq(nancat(F)); % trial, time, sub, drg
    F1=F1-F1(:,1,:,:,:); % subtract starting point
    F1(:,1,:,:,:) = complex(0,0);
    
    % F ( condition, time, subject, drug, trial )
    
    fInfo = catStruct(3,fInf,''); % info.amp ( trial, saccade, sub, drug )

    % only keep 1400ms
    tooLong = fInfo.rt > 1400;
    fnames = fieldnames(fInfo);
    for i = 1:length(fnames)
        fInfo.(fnames{i})(tooLong) = NaN;
    end
    
    % need to remove trials with debug or driftcorr
    
    % now do some filtering of bad trials:
%     angA = wrap(2*angle(A(:,end,:,:,:))); % sacc end angle where 0 is horizontal, +/- pi is vertical.
%     badA = any(isnan(A),2)      | any(abs(A)>400,2)      | any( abs(imag(A))>75 ,2)   | abs(angA)>1 ;
%     % A filtered:
%     Af = A + bool2nan(repmat(badA,[1,50])); % remove bad trials
%     Af(:,:,bad_sub,:,:) = nan; % remove bad subjects

    save(fixFile,'F','F1','fInfo','FInfo');
end


%% FFT on eye pos, plot average spectrum for early/late by condition

if exist('FixationData.mat', 'file')
    load('FixationData.mat');
else % if not already run, process them
    O = cell(nPP,3);
    parfor i = 1:nPP
    %     disp(i)
        for j = 1:3
            disp([i,j])
            if ~isempty(F{i,j})
                cond = tt(1:size(F{i,j},1),i,j);
                O{i,j} = analyseFixationPeriod(F{i,j}, FInfo{i,j}, cond, [1 1400],0); % entire period, no FFT
                fft1 = analyseFixationPeriod(F{i,j}, FInfo{i,j}, cond, [200 700],1); % early FFT
                fft2 = analyseFixationPeriod(F{i,j}, FInfo{i,j}, cond, [700 1200],1); % early FFT

                O{i,j}.fft1 = fft1.fft;
                O{i,j}.fft2 = fft2.fft;
    %             O{i,j} = O1;
            end
        end
    end

    save('FixationData.mat','O');
end
%% extract

for i = 1:nPP
    for j=1:nConds
        if ~isempty(O{i,j})
            O2(i,j) = O{i,j};
        end
    end
end

[microDens, meanSpeedCond, fft1, fft2] = deal([]);
for j = 1:3
    microDens = nancat(4, microDens, nancat(3, O2(:,j).microDensity));
    meanSpeedCond = nancat(5, meanSpeedCond, nancat(4, O2(:,j).meanSpeedCond));
    fft1 = nancat(5, fft1, nancat(4, O2(:,j).fft1));
    fft2 = nancat(5, fft2, nancat(4, O2(:,j).fft2));
end

microDens = permute(microDens, [3,1,2,4]);
meanSpeedCond = permute(meanSpeedCond, [4,1,3,5,2]);
fft1 = permute(fft1, [4,1,3,5,2]);
fft2 = permute(fft2, [4,1,3,5,2]);

N = 500;
srate = 1000; %Hz
nyquist = srate/2;
freqs = linspace(0,nyquist,floor(N/2)+1);


%% stats - use permOLS

for i = 1:4
    [~,fixP(i,:)] = permutationOLS(-diff(microDens(:,:,i,1:2),[],4), [],[],[],'cluster',true,'clustermethod','mean','two_tailed',true);
end
    

%% plot
n=2;%2=PD ON vs OFF, 3=HC too
figure();
doPerm = 1;
useClust = 1;
nPerms = 5000;
plotfreq = 5:125;
for i = 1:4
    subplot(4,4,i)
    errorBarPlot(sq(microDens(:,:,i,1:n)),'area',1, 'doStats', 0);%,'xaxisvalues',xi);
    title(cond_names{i});
    if i==1; ylabel ('microsaccades'); xlabel('time in fixation (ms)'); end
    set(gca,'XTick',0:50:100,'XTickLabel',0:700:1400)
    if doPerm
        [~,p]=permutationOLS( -diff(microDens(:,:,i,1:2),[],4), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p);
    end
    
    subplot(4,4,4+i);
    errorBarPlot(sq(nanmean(meanSpeedCond(:,:,i,1:n,:),5)),'area',1, 'doStats', 0);
    if i==1; ylabel('drift speed');xlabel('time in fixation (ms)'); end
    set(gca,'XTick',0:700:1400,'YTick',0:.1:1);
    if doPerm
        [~,p]=permutationOLS( -diff(nanmean(meanSpeedCond(:,:,i,1:2,:),5),[],4), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p,'yVal',0.3);
    end
    
    subplot(4,4,8+i); 
%     plotfreq=50:250; %  don't report ultra-high or low frequencies
    errorBarPlot(sq(nanmean(log(fft1(:,plotfreq,i,1:n,:)),5)),'area',1, 'xaxisvalues',freqs(plotfreq), 'doStats', 0);
    if i==1; ylabel 'early log power'; xlabel 'frequency (Hz)'; end
%     set(gca,'XTick',50:100:250)
    if doPerm
        [~,p]=permutationOLS( -diff(nanmean(log(fft1(:,plotfreq,i,1:n,:)),5),[],4), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p);
    end
    
    subplot(4,4,12+i);
%     plotfreq=50:250; %  don't report ultra-high or low frequencies
    errorBarPlot(sq(nanmean(log(fft2(:,plotfreq,i,1:n,:)),5)),'area',1, 'xaxisvalues',freqs(plotfreq), 'doStats', 0);
    if i==1; ylabel 'late log power'; xlabel 'frequency (Hz)'; end
%     set(gca,'XTick',50:100:250)
    if doPerm
        [~,p]=permutationOLS( -diff(nanmean(log(fft2(:,plotfreq,i,1:n,:)),5),[],4), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p);
    end
end
drawnow

for i = 1:4, makeSubplotScalesEqual(4,4,[i*4-3:i*4]); end

subplot(4,4,16); 
h = findobj(gca,'Type','Line');
legend(flipud(h),[on_off(1:n) ],'Location','Best');

%% plot rew effects for each condition

s = 10; % movmean window size
s2 = 10; % window size for fft
n=2;
figure();
doPerm = 1;
useClust = 1;
motNames = {'Contingent', 'Guaranteed'};
nPerms = 5000;
plotfreq = 1:length(freqs);
for i = 1:2
    subplot(4,2,i)
    y = movmean(sq(-diff(microDens(:,:,i*2-1:i*2,1:n),[],3)),s,2);
    errorBarPlot(y,'area',1, 'doStats', 0);%,'xaxisvalues',xi);
    title([motNames{i} ' effect']);
    if i==1; ylabel ('microsaccades'); xlabel('time in fixation (ms)'); end
    set(gca,'XTick',0:50:100,'XTickLabel',0:700:1400)
    yline(0,':k');
    if doPerm
        y1 = diff(y,[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, x,[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    
    
    subplot(4,2,2+i); 
    y = movmean(sq(-diff(nanmean(meanSpeedCond(:,:,i*2-1:i*2,1:n,:),5),[],3)),s,2);
    errorBarPlot(y,'area',1, 'doStats', 0);%,'xaxisvalues',xi);
    if i==1; ylabel('drift speed');xlabel('time in fixation (ms)'); end
    yline(0,':k');
    xlim([0 1400]);
    set(gca,'YTick',-1:.1:1,'XTick',0:700:1400);
    if doPerm
        y1 = diff(y,[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, x,[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    
    subplot(4,2,4+i); 
    y = movmean(sq(nanmean(-diff(log(fft1(:,plotfreq,i*2-1:i*2,1:n,:)),[],3),5)),s2,2);
    errorBarPlot(y,'area',1, 'xaxisvalues',freqs(plotfreq), 'doStats', 0);
    if i==1; ylabel 'early log power'; xlabel 'frequency (Hz)'; end
    yline(0,':k');
    set(gca,'XTick',0:250:500)
    if doPerm
        y1 = diff(y,[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, x,[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    
    subplot(4,2,6+i);
    y = movmean(sq(nanmean(-diff(log(fft2(:,plotfreq,i*2-1:i*2,1:n,:)),[],3),5)),s2,2);
    errorBarPlot(y,'area',1, 'xaxisvalues',freqs(plotfreq), 'doStats', 0);
    if i==1; ylabel 'late log power'; xlabel 'frequency (Hz)'; end
    yline(0,':k');
    if doPerm
        y1 = diff(y,[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1, x,[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    set(gca,'XTick',0:250:500)
end
drawnow

for i = 1:4, makeSubplotScalesEqual(4,2,[i*2-1:i*2]); end

subplot(4,2,8); 
h = findobj(gca,'Type','Line');
legend(flipud(h),on_off(1:n),'Location', [0.7772 0.2584 0.1268 0.0726]);

%% mean microsaccade density across 1400ms
dens = permute(reshape(microDens(:,:,:,1:3),[nPP,100,2,2,3]),[1,3,4,5,2]);

meanDens = nanmean(dens, 5); % mean across time

vnames = {'pp','cont','mot','drug'};
meanDensStats = rmanova(meanDens(:,:,:,1:2), vnames, 'categorical',[2 3 4],'DummyVarCoding','effects');


figure();
cols = get(gca,'ColorOrder'); cols = repmat(cols(1:3,:), 10,1);
condCols = [1 2; 3 4];
xvals = [1:4]' + [0 .02 .05];
for j = 1:2
    set(gca,'ColorOrder',cols);
    errorBarPlot(sq(nanmean(microDens(:,:,condCols(j,:),1:2),2)),'xaxisvalues',xvals(j*2-1:j*2,1:2),'plotargs',{'LineWidth',2}, 'doStats', 0);
    hold on;
    errorBarPlot(sq(nanmean(microDens(:,:,condCols(j,:),3),2)),'xaxisvalues',xvals(j*2-1:j*2,3),'plotargs',{'LineWidth',2}, 'doStats', 0);
end
set(gca,'XTick',1:4, 'XTickLabel',cond_names);
xlim([.5 4.5])
ylabel('# microsaccades')
legend(on_off, 'Location','Best')
box off

%% mean drift speed across 1400ms

meanDrift = reshape(sq(nanmean(meanSpeedCond,[5,2])), [nPP,2,2,3]);

meanDriftStats = rmanova(meanDrift(:,:,:,1:2), vnames, 'categorical', [2 3 4],'DummyVarCoding','effects');


figure();
cols = get(gca,'ColorOrder'); cols = repmat(cols(1:3,:), 10,1);
condCols = [1 2; 3 4];
xvals = [1:4]' + [0 .02 .05];
for j = 1:2
    set(gca,'ColorOrder',cols);
    errorBarPlot(sq(meanDrift(:,:,j,1:2)),'xaxisvalues',xvals(j*2-1:j*2,1:2),'plotargs',{'LineWidth',2}, 'doStats', 0);
    hold on;
    errorBarPlot(sq(meanDrift(:,:,j,3)),'xaxisvalues',xvals(j*2-1:j*2,3),'plotargs',{'LineWidth',2}, 'doStats', 0);
end
set(gca,'XTick',1:4, 'XTickLabel',cond_names);
xlim([.5 4.5])
ylabel('mean ocular drift speed (deg/s)')
legend(on_off,'Location','Best')
box off
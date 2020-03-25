% ContingentAnalysis.m
% Run the main analysis, draws figures, does statistics.
% loads up ContingentPDData.mat
% Calls all the extra analyses at the bottom

close all;clc; clear all;

%% load up data
if exist('./ContingentPDData.mat','file')
    load('./ContingentPDData.mat','result','cond_names','nTr', 'resultHC', 'nTrHC');
else
    ContingentDataLoad();
    load('./ContingentPDData.mat','result','cond_names','nTr', 'resultHC', 'nTrHC');
end


%% combine PD + HC?
nPD = size(result,1); % num PD patients

doHC = 1; % include HC in analysis?
if doHC %combine HC with PD
    nHC = size(resultHC, 1);
    result = nancat(2, result, resultHC);
    nConds = 3;
    on_off = {'ON','OFF','HC'}; % condition names
    nTrAll = nancat(2, nTr, nTrHC); % number of trials
    nPP = max(nPD, nHC); % num participants
else
    nConds = 2;
    nTrAll = nTr;
    nPP = nPD;
    on_off = {'ON','OFF'};
end

cg = {'Contingent','Guaranteed'}; % condition labels
%
dpp = 0.0311; % degrees per pixel

%% concatenate data into a nd-array
% this facilitates removing particular saccades within a trial 
clear inf tt rw
for i=1:nPP % for each dataset
  for j=1:nConds
    if ~isempty(result{i,j})
    inf(i,j).amp = nancat( 1, result{i,j}.sAmpl  ).*dpp; % ok, bad idea to use inf as a variable, but...
    inf(i,j).vel = nancat( 1, result{i,j}.sSpd   ).*dpp;
    inf(i,j).rt  = nancat( 1, result{i,j}.sRT    );
    inf(i,j).ep  = nancat( 1, result{i,j}.sEndpt ).*dpp;
    inf(i,j).bln = nancat( 1, result{i,j}.sBlink );
    
    % trial types
    tt{i,j}  = nancat( 1, result{i,j}.trialType );
    rw{i,j}  = nancat( 1, result{i,j}.reward    ); 
    
    % reflect x-axis around 512, keep y unchanged
    targDir = repmat(nancat( 1, result{i,j}.targetDir ), 1, size(inf(i,j).amp,2));
    inf(i,j).epRefl = real(inf(i,j).ep);
    inf(i,j).epRefl(targDir==1) = - inf(i,j).epRefl(targDir==1) + (1024 .* dpp);
    inf(i,j).epRefl = complex(inf(i,j).epRefl, imag(inf(i,j).ep));
    end
  end
end
tt=sq(nancat(tt));   % trial_type ( trial, sub, drug ) = 1/2/3/4
rw=sq(nancat(rw));   % reward = how much they actually won = 0/10

%%
info = catStruct(3,inf,''); % info.amp ( trial, saccade, sub, drug )
% exclude bad trials
% note: the criteria can actually make a difference to the statistics.

% criteria options:
% remove subjects with <10 trials / condition
EXCLUDE_SUB_MIN_TRIALS = 10;   % put zero to include everyone (ok for mixed model)

% note: result_ij.result.params.targetSize      = 70
% and   result_ij.result.params.targetLocations = +/- 300
% so an amplitude of 400 is significant overshoot
MAX_AMP           = 400.*dpp;    % the target is 300 px
MIN_AMP           = 1; % 1 vis deg
EXCLUDE_BAD_VEL   = true;   % whether to discard velocity < 100 deg/s 

MODE = 'not-firstonly'; % which exclusion criteria mode to use (See switch below)
switch MODE
  case 'firstonly' 
    % only use first saccade of the trial; discard trial if it's not valid
    FIRST_ONLY=true; 
    ok = info.amp > MIN_AMP...
       & info.amp < MAX_AMP ...
       & info.rt > 100 ...     
       & info.rt < 600 ...
       ;
  case 'not-firstonly'
    FIRST_ONLY        = false;   % use the first valid saccade of the trial
    ok = info.amp > MIN_AMP...
       & info.amp < MAX_AMP ...
       & info.rt > 100 ...     
       & info.rt < 600 ...
       ;
end

% use 
ok2 = info.amp > MIN_AMP & info.amp < MAX_AMP & info.rt > 100 & info.rt < 600;
ok(:,:,:,3) = ok2(:,:,:,3);

if EXCLUDE_BAD_VEL
  ok = ok & info.vel > (2.5.*dpp)   & info.vel< (80.*dpp);
end

% make velocity of bad saccades nan
info.vel( ~ok ) = nan;
% which of the saccades in a trial is the 'first good' one?
good_idx = sq(findn(permute( ~isnan(info.vel), [2,1,3,4,5] ), 1 ));       

% actually remove the bad saccades:
info = stripnan( info, 2 ); 

fn=fieldnames(info);
for i=1:length(fn) % for each field, get first good saccade
  % select just first good saccade in each trial
  info.(fn{i}) = sq(info.(fn{i})(:,1,:,:,:));
  % exclude any trials where the accepted saccade wasn't the first one?
  if FIRST_ONLY
    info.(fn{i}) = info.(fn{i}) + bool2nan( good_idx ~= 1 );
  end
end

% Now group trials by condition.
% info groupings: ( cond, subj, drug, trial )
% infg = groupMeans(info, 1, repmat( permute(tt, [1 4 2 3]), [1,3] ) ,'dim');
infg = groupMeans(info,1, tt, 'dim'); 
ok_g = groupMeans(sq(ok(:,1,:,:)),  1, tt, 'dim'); % (cond,sub,drg,trial)
% save concatenated_saccades info tt

amp = permute(infg.amp,[4,2,1,3]); % amp ( trial, subj, cond4, onoff )
rt  = permute(infg.rt ,[4,2,1,3]);
vel = permute(infg.vel,[4,2,1,3]);

% endpoint deviance = absolute distance from centre of closest target
infg.epDiv = abs(infg.epRefl-repmat(nanmean(infg.epRefl,4),[1,1,1,size(infg.epRefl,4)]));
epvar2d = permute(infg.epDiv, [4,2,1,3]);

% calculate residuals relative to the main sequence
clear vel_res b_ij
for i=1:size(info.amp,2) % subj
  for j=1:nConds % group/condition
    [b_ij(i,j,:), ~, vel_res{i,j}] = regress( info.vel(:,i,j),  [ ones(size(info.amp,1),1)  info.amp(:,i,j) ] );
  end
end
vel_res = sq(nancat(vel_res)); % trial, sub, drug
% now group by condition
vr = permute( groupMeans(vel_res, 1, tt, 'dim'), [4,2,1,3]);
% vr ( trial, sub, cond, drug )]

infg.vr = permute(vr, [3,2,4,1]);
% exclude subjects with < 10 trials in any condition
bad_sub = any( any ( (sum(~isnan(vr(:,:,:,1:2)))) < EXCLUDE_SUB_MIN_TRIALS , 3), 4 ); 
bad_sub2 = any( any ( (sum(~isnan(vr(:,:,:,3)))) < EXCLUDE_SUB_MIN_TRIALS , 3), 4 ); 
fn=fieldnames(infg);
if any(bad_sub) % EXCLUDE subjects - note not required if using mixed model
  for i=1:length(fn)
    infg.(fn{i})(:,bad_sub,1:2,:) = nan;
    infg.(fn{i})(:,bad_sub2,3,:) = nan;
  end
%   vr(:,bad_sub,1:2,:) = nan;
%   vr(:,bad_sub2,3,:) = nan;

  fprintf('%g subjects excluded\n',sum(bad_sub + bad_sub2));
end


vr = permute(infg.vr, [4,2,1,3]);

prop_OK = sq(sum(~isnan(info.vel))) ./ nTrAll; % proportion OK, per participant
% show what proportion of trials were ok (just nonexcluded participants)
mean_OK = mean(mean(prop_OK([~bad_sub ~bad_sub ~bad_sub2])))  
%% stats
% put these into a Table and do stats: 
% make a table for fitlme:
t.ones = ones(size(vr(:,:,:,1:2)));
t.sub  = repmat( [1:size(vr,2)], [ size(vr,1), 1, size(vr,3),  size(vr,4)-1 ] );
t.mot  = repmat( permute([1;-1;1;-1], [2,3,1,4]) , [ size(vr,1), size(vr,2), 1,size(vr,4)-1 ] );
t.cont = repmat( permute([1;1;-1;-1], [2,3,1,4]) , [ size(vr,1), size(vr,2), 1, size(vr,4)-1 ] );
t.drg  = repmat( permute([1;-1], [2,3,4,1]) ,      [ size(vr,1), size(vr,2), size(vr,3), 1 ] );
t.vr   = vr(:,:,:,1:2);
t.vel  = permute( infg.vel(:,:,1:2,:), [4,2,1,3] );
t.amp  = permute( infg.amp(:,:,1:2,:), [4,2,1,3] );
t.rt   = permute( infg.rt(:,:,1:2,:) , [4,2,1,3] );
t.epv  = permute( infg.epDiv(:,:,1:2,:), [4,2,1,3] );
fn=fieldnames(t);
T=table();
for i=1:length(fn)
  T.(fn{i}) = t.(fn{i})(:);
end

dv = {'vr','amp','rt', 'epv', 'vel'};
for i = 1:length(dv)

    %%%%%
    % Statistical method 1: 
    % full linear mixed effects model using all trials
%     M{i} = fitlme( T, sprintf('%s ~ 1 + mot * cont * drg + (1|sub)',dv{i}) );

    %%%%%
    % Statistical method 2:
    % repeated measures ANOVA on per-condition means:
    anTab{i} = rmanova( reshape(nanmean(t.(dv{i})),nPP,2,2,2), {'s','mot','cont','drg'} , 'categorical', 4 );
    pvalues(1:8,i) = anTab{i}.pValue';
%     pvalues(9,i) = M{i}.Coefficients.pValue(8);
    
end


%% two two-way anova

for i = 1:length(dv)


    % contingency
        % full linear mixed effects model using all trials
%     M2{i,1} = fitlme( T, sprintf('%s ~ 1 + cont * drg + (1|sub)',dv{i}) );

    
    % repeated measures ANOVA on per-condition means:
    anTab2{i,1} = rmanova(reshape(nanmean(t.(dv{i})(:,:,1:2,:)),nPP,2,2), {'s','cont','drg'}, 'categorical', 3 );
    pvalues2(1:4,i,1) = anTab2{i,1}.pValue';
%     pvalues2(5,i,1) = M2{i,1}.Coefficients.pValue(4);
    
    % certain rewards 
    % full linear mixed effects model using all trials
    M2{i,2} = fitlme( T, sprintf('%s ~ 1 + mot * drg + (1|sub)',dv{i}) );
    
    % repeated measures ANOVA on per-condition means:
    anTab2{i,2} = rmanova(reshape(nanmean(t.(dv{i})(:,:,3:4,:)),nPP,2,2), {'s','mot','drg'}, 'categorical', 3 );
    pvalues2(1:4,i,2) = anTab2{i,2}.pValue';
%     pvalues2(5,i,2) = M2{i,2}.Coefficients.pValue(4);
    

end

%% two two-way anova on each drug

for i = 1:length(dv)


    % contingency
        % full linear mixed effects model using all trials
%     M3{i,1} = fitlme( T, sprintf('%s ~ 1 + mot * cont + (1|sub)',dv{i}) );

    
    % repeated measures ANOVA on per-condition means:
    anTab3{i,1} = rmanova(reshape(nanmean(t.(dv{i})(:,:,:,1)),nPP,2,2), {'s','mot', 'cont'}, 'categorical', 3 );
    pvalues3(1:4,i,1) = anTab3{i,1}.pValue';
%     pvalues3(5,i,1) = M3{i,1}.Coefficients.pValue(4);
    
    % certain rewards 
    % full linear mixed effects model using all trials
%     M3{i,2} = fitlme( T, sprintf('%s ~ 1 + mot * drg + (1|sub)',dv{i}) );
    
    % repeated measures ANOVA on per-condition means:
    anTab3{i,2} = rmanova(reshape(nanmean(t.(dv{i})(:,:,:,2)),nPP,2,2), {'s','mot', 'cont'}, 'categorical', 3 );
    pvalues3(1:4,i,2) = anTab3{i,2}.pValue';
%     pvalues3(5,i,2) = M3{i,2}.Coefficients.pValue(4);
    

end

%% HC

if doHC
    tHC.ones = ones(size(vr));
    tHC.sub  = repmat( [1:size(vr,2)], [ size(vr,1), 1, size(vr,3),  size(vr,4) ] );
    tHC.mot  = repmat( permute([1;-1;1;-1], [2,3,1,4]) , [ size(vr,1), size(vr,2), 1,size(vr,4) ] );
    tHC.cont = repmat( permute([1;1;-1;-1], [2,3,1,4]) , [ size(vr,1), size(vr,2), 1, size(vr,4) ] );
    tHC.cond  = repmat( permute([1;-1;2], [2,3,4,1]) ,      [ size(vr,1), size(vr,2), size(vr,3), 1 ] );
    tHC.vr   = vr;
    tHC.vel  = permute( infg.vel, [4,2,1,3] );
    tHC.amp  = permute( infg.amp, [4,2,1,3] );
    tHC.rt   = permute( infg.rt , [4,2,1,3] );
    tHC.epv  = permute( infg.epDiv, [4,2,1,3] );
    
    pvaluesHC = NaN(8,5,4);

    for i = 1:length(dv)

        %%%%%
        % Statistical method 1: 
        % full linear mixed effects model using all trials
    %     MHC{i} = fitlme( T, sprintf('%s ~ 1 + mot * cont * cond + (1|sub)',dv{i}) );

        %%%%%
        % Statistical method 2:
        % repeated measures ANOVA on per-condition means:
        % ON vs OFF vs HC
        anTabHC{i,1} = rmanova( reshape(nanmean(tHC.(dv{i})),nPP,2,2,3), {'s','mot','cont','group'}, 'categorical', 4 );
        pvaluesHC(1:8,i,1) = anTabHC{i,1}.pValue';
    %     pvaluesHC(9,i) = M2{i}.Coefficients.pValue(8);


        % ON vs HC
        anTabHC{i,2} = rmanova( reshape(nanmean(tHC.(dv{i})(:,:,:,[1 3])),nPP,2,2,2), {'s','mot','cont','group'},'categorical',4 );
        pvaluesHC(1:8,i,2) = anTabHC{i,2}.pValue';

        % OFF vs HC
        anTabHC{i,3} = rmanova( reshape(nanmean(tHC.(dv{i})(:,:,:,[2 3])),nPP,2,2,2), {'s','mot','cont','group'},'categorical',4 );
        pvaluesHC(1:8,i,3) = anTabHC{i,3}.pValue';
        
        
        % just HC
        anTabHC{i,4} = rmanova( reshape(nanmean(tHC.(dv{i})(:,:,:,[3])),nPP,2,2), {'s','mot','cont'} );
        pvaluesHC(1:4,i,4) = anTabHC{i,4}.pValue';
    end
end
%% %%%%%%% plot everything  

nDVs = 4; % how many DVs to plot

ylabs = {'Peak Velocity (deg/s)', 'Amplitude (deg) ','Saccadic RT (ms)',...
    'Endpoint Varibility (deg)','Peak Velocity (deg/s)'};
xlabs = {'Perform','Random','+10p','0p'};
xTitle = {'Contingent              Reward', '    Motivation           Expectation'};
plotargs = {'LineWidth',2.5};
condCols = [1 2; 3 4];
cols = get(gca,'ColorOrder'); %cols = cols([3 1 2], :);
xvals = [1:4]' + [-.1 0 .1];

figure();
for i = 1:nDVs
    subplot(ceil(nDVs/2),2,i)
    if doHC % get variable
        data = sq(nanmean(tHC.(dv{i})));
    else
        data = sq(nanmean(t.(dv{i})));
    end
    
    % plot each pair of conditions
    for j = 1:2
        set(gca,'ColorOrder',cols); % keep colour order
        if doHC
            h(3) = errorBarPlot( data(:,condCols(j,:),3), 'xaxisvalues', xvals(j*2-1:j*2,3), 'plotargs',plotargs, 'doStats', 0);
            hold on;
        end
        h(1:2) = errorBarPlot( data(:,condCols(j,:),1:2), 'xaxisvalues', xvals(j*2-1:j*2,1:2), 'plotargs',plotargs, 'doStats', 0);
        for k = 1:2+doHC
            h(k).Color = cols(k,:);
        end
    end
    ylabel(ylabs{i})
    set(gca,'xtick',1:4,'xticklabel',xlabs);
    xlim([.75 4.25])
    box off
    xlabel(xTitle);
    
end
legend(h, on_off, 'location','Best')

%% plot with individual points

plotargs = {'LineWidth',2.5,'Marker','x','MarkerSize',4};
xTitle = {'Contingent              Reward', '    Motivation           Expectation'};
cols = get(gca,'ColorOrder'); %cols = cols([3 1 2], :);
xvals = [1:4]' + [-.15 0 .15];
figure();
errbars = 0; %0= none, 1=SEM, 2=SD
clear h;
for i = 1:nDVs
    subplot(ceil(nDVs/2),2,i)
    if doHC
        data = sq(nanmean(tHC.(dv{i})));
    else
        data = sq(nanmean(t.(dv{i})));
    end
        
    if ~errbars
    %%%    plot means - no error bars
        hold on;
        for j = 1:2
            h(:,j) = plot(xvals(j*2-1:j*2 ,:), sq(nanmean(data(:,j*2-1:j*2,:))),'-','LineWidth',2.5);
            for k=1:(2+doHC)
                h(k,j).Color = cols(k,:);
            end
        end
        plot(xvals, sq(nanmean(data)),'k.','MarkerSize',12);
    else
    %%%have error bars still
        for j = 1:2
            set(gca,'ColorOrder',cols);
            if doHC
                h(3,1) = errorBarPlot( data(:,condCols(j,:),3),'standardError',errbars, 'xaxisvalues', xvals(j*2-1:j*2,3), 'plotargs',plotargs, 'doStats', 0);
                hold on;
            end
            h(1:2,1) = errorBarPlot( data(:,condCols(j,:),1:2),'standardError',errbars, 'xaxisvalues', xvals(j*2-1:j*2,1:2), 'plotargs',plotargs, 'doStats', 0);
            for k = 1:2+doHC
                h(k).Color = cols(k,:);
                h(k).MarkerEdgeColor = [0 0 0];
            end
        end
    end
    
    for k=1:(2+doHC)
        for j=1:4
            h1 = scatter(repmat(xvals(j,k),1,nPP) + rand(1,nPP)*.05-.025,data(:,j,k)', 36, cols(k,:));
            h1.Marker = 'x';
            alpha(h1, 0.5);
        end
%         h1(k,:) = plot(repmat(xvals(:,k),1,nPP) + rand(4,nPP)*.05-.025,data(:,:,k)','x','Color',[0.1 0.1 0.1]);%cols(k,:));
%         hold on;
    end
    ylabel(ylabs{i})
    set(gca,'xtick',1:4,'xticklabel',xlabs);
    xlim([.75 4.25])
    box off
    xlabel(xTitle);
end
legend(h(:,1), on_off, 'location','Best')

%% Between subject correlations
figure()
% run this after the previous section to get vr.
% calculate motivation effect = difference in mean peak velocity 
% mvr = mean vel resid ( sub, mot_level, cont/guar, drug )
mvr = reshape(nanmean(vr(:,:,:,:)),[],2,2,nConds);
% dmvr = difference in mvr with motivation ( subj, cont/guar, on/off )
dmvr = sq( diff( mvr, [],2 ) );

plotargs = {'pearson',0,'plot_ci',0,'text',2, 'plotline', 2, 'showzero',0}; % use spearman correlation?

set(0,'DefaultAxesColorOrder',[0 0 0; 0 0 0]);
for i = 1:3
    subplot(2,3,i)
    x = dmvr(:,1:2,i); x(all(isnan(x),2),:) = [];
    [~, ~, ~, ~, ~, h] = scatterRegress( x(:,1), x(:,2)  , plotargs{:});
    hold on; yline(0, 'k:'); xline(0,':k');
    xlabel(['Contingent Effect: ' on_off{i}]);  ylabel(['Guaranteed Effect: ' on_off{i}]);
    set(gca,'XTick',-.1:.1:.1,'YTick',-.1:.1:.1);
    axis([-.15 .15 -.15 .15])
    h(1).MarkerEdgeColor = cols(i,:);
    h(1).Marker = 'x';
    axis('square')
end


% now plot differences 
subplot(2,3,4)
x = sq(dmvr(:,1,1:2)); x(all(isnan(x),2),:) = [];
[~, ~, ~, ~, ~, h] = scatterRegress( x(:,1), x(:,2)  , plotargs{:});hold on; yline(0, 'k:'); xline(0,':k');
xlabel 'Contingent Effect PD ON'; ylabel 'Contingent Effect PD OFF';
set(gca,'XTick',-.1:.1:.1,'YTick',-.1:.1:.1);
axis([-.15 .15 -.15 .15])
h(1).MarkerEdgeColor = 'k';
axis('square')

subplot(2,3,5)
x = sq(dmvr(:,2,1:2)); x(all(isnan(x),2),:) = [];
[~, ~, ~, ~, ~, h] = scatterRegress( x(:,1), x(:,2)  , plotargs{:});hold on; yline(0, 'k:'); xline(0,':k');
hold on; yline(0, 'k:'); xline(0,':k');
xlabel 'Guaranteed Effect PD ON'; ylabel 'Guaranteed Effect PD OFF';
set(gca,'XTick',-.1:.1:.1,'YTick',-.1:.1:.1);
axis([-.15 .15 -.15 .15])
h(1).MarkerEdgeColor = 'k';
axis('square')

subplot(2,3,6)
x = [dmvr(:,1,1)-dmvr(:,1,2), dmvr(:,2,1)-dmvr(:,2,2)]; x(all(isnan(x),2),:) = [];
[~, ~, ~, ~, ~, h] = scatterRegress( x(:,1), x(:,2)  , plotargs{:});hold on; yline(0, 'k:'); xline(0,':k');
hold on; yline(0, 'k:'); xline(0,':k');
xlabel 'Contingent Effect: PD ON - OFF'; ylabel 'Guaranteed Effect: PD ON - OFF';
set(gca,'XTick',-.1:.1:.1,'YTick',-.1:.1:.1);
axis([-.15 .15 -.15 .15])
h(1).MarkerEdgeColor = 'k';
axis('square')
makeSubplotScalesEqual(2,3);


set(0,'DefaultAxesColorOrder', 'factory'); % restore


  
%% save
save('ContingentAnalysis.mat','infg','info','t','tt','rw','dmvr','vr','doHC',...
    'nConds','on_off','nTrAll','cond_names','dpp','nPD','bad_sub', 'cg',...
    'xlabs', 'ok_g','bad_sub2','nPP','anTab','anTabHC','anTab2','anTab3','nConds',...
    'tHC','dv');


%% extra analyses
% 
ContingentRawVelocity
ContingentPupilCue
ContingentCovarCorrel
ContingentVelocity
ContingentFixation

ContingentDemographics
ContingentUPDRS
ContingentDrugs

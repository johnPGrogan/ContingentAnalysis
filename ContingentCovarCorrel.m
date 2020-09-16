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
    Af(:,:,bad_sub,1:2,:) = complex(NaN, NaN); % remove bad subjects
    Af(:,:,bad_sub2,3,:) = complex(NaN,NaN); % remove bad subjects

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


%% use permutationOLS or uncorrected t-tests
usePermOLS = 1;
n = 2; % 2 = PD only, 3 = HC too
contourColour = [0 0 0];
contourWidth = 1.5;
%% cor  motivated minus not-motivated

% difference in corr/cov: dc ( sub, t1, t2, cont/guar, drug ) 
corDiff = diff(corAf,[],4); 
corDiff(isnan(corDiff)) = 0; % set nan to zero
corDiff = -corDiff(:,:,:,[1 3],:); % keep only mot effects
% rows are cont/guar, cols are on/off
corDiffMeans = sq(nanmean(corDiff)); % take mean across people
corDiffMeans = cat(4, -diff(corDiffMeans(:,:,:,1:2),[],4), corDiffMeans); % put ON-OFF first

conds = {'ON-OFF','ON','OFF','HC'};
cg = {'Contingent','Guaranteed'}; hilo = {'hi mot','lo mot'};

% draw figures
DrawAutocorrelFigures(corDiff, corDiffMeans, conds, [-.41 .41])

% save source data
readme = "Data to create Figure 4. Please download the ContingentAnalysis GitHub repo and Matlib repo (links available in paper), and then run: DrawAutocorrelFigures(corDiff, corDiffMeans, conds, [-.41 .41]);";
save('Figure4SourceData.mat', 'corDiff','corDiffMeans','conds','readme');

%% indiv plots

corDiff(:,:,:,:,4) = -diff(corDiff(:,:,:,:,1:2),[],5); % on-off

% change order + shape for plotn
corDiffInd = permute(corDiff,[2,3,1,4,5]);
corDiffInd = corDiffInd(:,:,:,:,[4,1,2,3]);

corDiffInd(:,:,:,:,1:3) = corDiffInd(:,:,[sort(find(~bad_sub)), sort(find(bad_sub))],:,1:3); % move missing pps to end
corDiffInd(:,:,:,:,4) = corDiffInd(:,:,[sort(find(~bad_sub2)), sort(find(bad_sub2))],:,4); % move missing pps to end

corDiffInd = reshape(corDiffInd,50,50,5,6,2,4);


xlabs = {'Contingent Effect', 'Guaranteed Effect'};
conds = {'ON - OFF', 'ON','OFF','HC'};
titlefun = @(i,j) ((j-1)*5+i);
eachfun = @(x) set(gca,'YTick', [1 50], 'YTickLabel', [0 100],'XTick', [1 50], 'XTickLabel', [0 100]);
% f = figure();
for i = 1:2
    for j = 1:4
%         clf
        f = figure(); 
%         f.WindowState = 'maximized';
        plotn( corDiffInd(:,:,:,:,i,j), 'plotfun',@(x) imagesc(x, [-1 1]),'fixshape',1,'titlefun',titlefun,'eachfun',eachfun);
                
        % colorbar
        subplot(5,6,1)
        c = colorbar('Location','West');
        c.Position = [0.0446 0.1476 0.0411 0.6921];
        c.Ticks = [-1 0 1];
        
        SuperTitle([xlabs{i} ': ' conds{j}]);
        saveas(f, sprintf('./Fig4_S2_%d.pdf', (i-1)*4+j));
    end
end



%% covariance
covAf = apply([2,5], @(x)nancov(x'), real(Af) ); % time-time corr ( cond, t1, sub, drug, t2 )
covAf = permute(covAf, [3,2,5,1,4]);

% use signed log covar for subtractions
covAf2 = covAf;
covAf = real(log(covAf2)) .* sign(covAf2);


%% cov  motivated minus not-motivated
% difference in corr/cov: dc ( sub, t1, t2, cont/guar, drug ) 
% change covAf to corAf to look at correlation
covDiff = diff(covAf,[],4); 
covDiff(isnan(covDiff)) = 0; % set nan to zero
covDiff = -covDiff(:,:,:,[1 3],:); % differences reward minus no-reward
% rows are cont/guar, cols are on/off
covDiffMeans = sq(nanmean(covDiff(:,:,:,:,1:3)));
covDiffMeans = cat(4, -diff(covDiffMeans(:,:,:,1:2),[],4), covDiffMeans);

DrawAutocorrelFigures(covDiff, covDiffMeans, conds, [-2, 2])


%%

save('ContingentCovarCorrel.mat', 'Af');
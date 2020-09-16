% ContingentModalityAnalysis

clc; clear all; close all;

%% get modality

load('./ContingentPDData.mat','result', 'resultHC')
result = nancat(2, result, resultHC);

load('ContingentAnalysis.mat')

%% extract modality 

for i=1:nPP % for each dataset
  for j=1:nConds
    if ~isempty(result{i,j})
        modality{i,j}  = nancat( 1, result{i,j}.subjectFeedbackCounterbalance);
    end
  end
end

modality = sq(nancat(modality));   % 1 = audio speed, 2 = visual speed

isVisSpeed = modality == 2;

%% anova comparing

% same dims as t.vr
isVisTT = permute(groupMeans(isVisSpeed,1,tt,'dim'), [4,2,1,3]);

%% anova with mot/cont/drug/modal

etaSqP = @(x) (x.FStat.*x.DF1) ./ (x.FStat.*x.DF1 + x.DF2); % function for partial eta squared

dv = {'vr','amp','rt', 'epv', 'vel'};
for i = 1:length(dv)

    % [pp, mot, cont, drug, modal]
    x = reshape(permute(groupMeans(t.(dv{i}), 1, isVisTT(:,:,:,1:2)), [2,3,4,1]),[nPP,2,2,2,2]);

    anTabModality{i} = rmanova( x, {'s','mot','cont','drug', 'modal'}, 'categorical', [2 3 4 5] ,'DummyVarCoding','effects');
    anTabModality{i}.etaSqP = etaSqP(anTabModality{i});

    % reshape into [pp*2, mot, cont, modal]
    x = reshape(permute(x,[1,4,2,3,5]), [], 2,2,2);
    anTabModality2{i} = rmanova( x, {'s','mot','cont', 'modal'}, 'categorical', [2 3 4] ,'DummyVarCoding','effects');
    
end

%% mot/cont/modal

for i = 1:length(dv)

    % [pp, mot, cont, drug, modal]
    x = reshape(permute(groupMeans(t.(dv{i}), 1, isVisTT(:,:,:,1:2)), [2,3,4,1]),[nPP,2,2,2,2]);

    % reshape into [pp*2, mot, cont, modal]
    x = reshape(permute(x,[1,4,2,3,5]), [], 2,2,2);
    anTabModality2{i} = rmanova( x, {'s','mot','cont', 'modal'}, 'categorical', [2 3 4] ,'DummyVarCoding','effects');
    
end


%% modal only

% get info.vr and info.epv

for i=1:size(info.amp,2) % subj
  for j=1:nConds % group/condition
    [~, ~, vel_res{i,j}] = regress( info.vel(:,i,j),  [ ones(size(info.amp,1),1)  info.amp(:,i,j) ] );
  end
end
vel_res = sq(nancat(vel_res));
info.vr = vel_res; 

% deviance from mean
info.epv = abs(info.epRefl - nanmean(info.epRefl,1));

for i = 1:length(dv)
    x = permute(groupMeans(info.(dv{i}), 1, isVisSpeed),[2,1,3]);
    
    anTabModality1{i} = rmanova( x, {'s','modal','grp'}, 'categorical', [1 2 3] ,'DummyVarCoding','effects');
end
    
    
    
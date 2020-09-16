% ContingentOrderEffects

clc; clear all; close all;

%%

load('./ContingentPDData.mat','timestamps')

timestamps = nancat(timestamps(:,1:2)); % reformat

t = nancat(cellfun(@datenum, timestamps,'UniformOutput',0)); % make into number

onIsFirst = t(:,1) < t(:,2); % order

% remove missing
onIsFirst(any(isnan(t),2)) = [];
timestamps(any(isnan(t),2),:) = [];

%% swap order for anova

load('ContingentAnalysis.mat')

%% anova with mot/cont/order

etaSqP = @(x) (x.FStat.*x.DF1) ./ (x.FStat.*x.DF1 + x.DF2); % function for partial eta squared

dv = {'vr','amp','rt', 'epv', 'vel'};
for i = 1:length(dv)
    x = reshape(nanmean(t.(dv{i})),nPP,2,2,2); % shape [pp, mot, cont, drug]
    x(~onIsFirst,:,:,:) = x(~onIsFirst,:,:,[2 1]); % swap on/off

    anTabOrder{i} = rmanova( x, {'s','mot','cont','order'}, 'categorical', [2 3 4] ,'DummyVarCoding','effects');
    anTabOrder{i}.etaSqP = etaSqP(anTabOrder{i});

end


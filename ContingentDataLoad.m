function ContingentDataLoad(folder)
% This was run by the authors already, on data that had been loaded via
% snipsaccades.m (snipSaccades(read,'starttarget_t','startfeedback_t');)
% So this is provided for reference only.
% 
% 
% load up the result data from the contingent task
% remove patients with only one session
% 

if ~exist('folder','var')
    folder = '../../../../OneDrive - Nexus365/Comp_Neurology_Lab/Data/PD_ON_OFF/Converted Contingent/';
end
%% load data
file_template = {
  [folder,  'PD_%03d_CTN_ON_data'];
  [folder,  '/PD_%03d_CTN_OFF_data'];
}  ;

ppIDs = cell(31,3);
timestamps = ppIDs;
for i=1:31 % for each subject
  for j=1:2 % on and off
    fn = sprintf(file_template{j},i); % create filename
    try
      result_ij=load(fn);  % load the file
      result{i,j} = result_ij.result.data; % extract the per-trial data
      timestamps{i,j} = result_ij.result.startTimes; % store timestamp for order
      ppIDs{i,j} = sprintf('PD_%03d',i);
    catch mexp
      fprintf('unable to read %s\n',fn);
    end
  end
end
% missing: 
%  13,  20, 21
% low trial numbers: 7 on; 23 off
bad_sub = any(cellfun(@isempty, result),2);
% exclude altogether the subjects with only one session.
% (thought there might be a way to include them using a mixed model)
result(bad_sub,:) = [];  
% useful labels
on_off     = {'ON','OFF'};
cond_names = {'cont','rand','guar+','guar0'};
nTr = cellfun(@length,result); % num trials ( sub, onoff )

%% load up HC too


files = what(folder);
files.mat = files.mat(~cellfun(@isempty,regexp(files.mat, 'EHC')));


for i=1:length(files.mat) % for each subject
  for j=1:1 % on and off
    fn = fullfile(files.path, files.mat{i}); % create filename
    try
      result_ij=load(fn);  % load the file
      resultHC{i,j} = result_ij.result.data; % extract the per-trial data
      timestamps{i,3} = result_ij.result.startTimes; % store timestamp for order
      ppIDs{i,3} = sprintf('EHC_%03d',str2num(files.mat{i}(regexp(files.mat{i}, '\d'))));
    catch mexp
      fprintf('unable to read %s\n',fn);
    end
  end
end
% missing: 
%  13,  20, 21
% low trial numbers: 7 on; 23 off
bad_sub = any(cellfun(@isempty, resultHC),2);
% exclude altogether the subjects with only one session.
% (thought there might be a way to include them using a mixed model)
resultHC(bad_sub,:) = [];  
% useful labels
nTrHC = cellfun(@length,resultHC); % num trials ( sub, onoff )

%% save (result is large variable)

save('./ContingentPDData.mat','-v7.3','result','on_off','cond_names','nTr', 'resultHC', 'nTrHC', 'ppIDs','timestamps');

end
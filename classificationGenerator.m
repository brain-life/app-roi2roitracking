function [] = classificationGenerator()

if ~isdeployed
    disp('loading path')
    addpath(genpath('/N/u/hayashis/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/encode'))
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/spm'))
    addpath(genpath('/N/u/brlife/git/wma'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
end

% Load configuration file
config = loadjson('config.json')

% Set tck file path/s
rois=dir('*.tck');

roiPair = config.roiPair;
    
for ii = 1:length(rois); 
    fgPath{ii} = rois(ii).name;
end

% Create classification structure
[mergedFG, classification]=bsc_mergeFGandClass(fgPath)

% Amend name of tract in classification structure
if isnumeric(roiPair(1))
    for ii = round((1:length(roiPair))/2)
        classification.names{ii} = strcat('ROI_',num2str(roiPair((2*ii) - 1)),'_ROI_',num2str(roiPair((2*ii))));
    end
else
    roiPair = split(roiPair);
    for ii = round((1:length(roiPair))/2)
        classification.names{ii} = strcat('ROI_',roiPair{(2*ii) - 1},'_ROI_',roiPair{(2*ii)});
    end
end

fgWrite(mergedFG, 'track/track.tck', 'tck')

% Create fg_classified structure
fg_classified = bsc_makeFGsFromClassification_v4(classification,mergedFG);

if ~exist('wmc', 'dir')
    mkdir('wmc')
end
if ~exist('wmc/tracts', 'dir')
    mkdir('wmc/tracts')
end

% Save output
save('wmc/classification.mat','classification','fg_classified','-v7.3');

% Create structure to generate colors for each tract
tracts = fg2Array(fg_classified);

% Make colors for the tracts
%cm = parula(length(tracts));
cm = distinguishable_colors(length(tracts));
for it = 1:length(tracts)
   tract.name   = strrep(tracts{it}.name, '_', ' ');
   all_tracts(it).name = strrep(tracts{it}.name, '_', ' ');
   all_tracts(it).color = cm(it,:);
   tract.color  = cm(it,:);

   %tract.coords = tracts(it).fibers;
   %pick randomly up to 1000 fibers (pick all if there are less than 1000)
   fiber_count = min(1000, numel(tracts{it}.fibers));
   tract.coords = tracts{it}.fibers(randperm(fiber_count)); 
   
   savejson('', tract, fullfile('wmc','tracts', sprintf('%i.json',it)));
   all_tracts(it).filename = sprintf('%i.json',it);
   clear tract
end

% Save json outputs
savejson('', all_tracts, fullfile('wmc/tracts/tracts.json'));

% Create and write output_fibercounts.txt file
for ii = 1 : length(fg_classified)
    name = fg_classified{ii}.name;
    num_fibers = length(fg_classified{ii}.fibers);
    
    fibercounts(ii) = num_fibers;
    tract_info{ii,1} = name;
    tract_info{ii,2} = num_fibers;
end

T = cell2table(tract_info);
T.Properties.VariableNames = {'Tracts', 'FiberCount'};

writetable(T, 'output_fibercounts.txt');


exit;
end




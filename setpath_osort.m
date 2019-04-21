%
% Run this file to initialize the paths to be able to run osort v4.
%

%=== Windows
basepathCode=[ '/Users/huimin/Documents/Matlab/osort-v4-rel' filesep];   % change this path accordingly

%=== Unix
%basepathCode=[ '/home/urut/svnwork/neuro1/code/osort-v4-rel' filesep];   % change this path accordingly


%=================== Below, no changes necessary
path(path,[basepathCode ]);
path(path,[basepathCode '/code/continuous/']);
path(path,[basepathCode '/code/continuous/neuralynx']);
path(path,[basepathCode '/code/continuous/txt']);
path(path,[basepathCode '/code/sortingNew/']);
path(path,[basepathCode '/code/sortingNew/projectionTest']);
path(path,[basepathCode '/code/sortingNew/model']);
path(path,[basepathCode '/code/sortingNew/model/detection']);
path(path,[basepathCode '/code/sortingNew/evaluation']);
path(path,[basepathCode '/code/osortTextUI']); % text user interface
path(path,[basepathCode '/code/osortGUI']); %graphical user interface
path(path,[basepathCode '/code/helpers']);
path(path,[basepathCode '/code/plotting']);

path(path,[basepathCode '/code/3rdParty/gabbiani']);
path(path,[basepathCode '/code/3rdParty/MClust-3.5/ClusterQuality/']);

if ispc
    %Windows version
    path(path,[basepathCode '/code/3rdParty/neuralynxWindows']);
else
    %Unix version
    path(path,[basepathCode '/code/3rdParty/neuralynxUnixAll/binaries']);    
end

path(path,[basepathCode '/code/3rdParty/cwtDetection']);

function [neuronId] = querySomaVoxelLoc(locStr,varargin)
%queryVoxelLocation. Returns list of neurons in specific anatomical
%location based on voxel search. Can include dilation of search area.

%% Read settings file.
[cFolder,~,~] = fileparts(which('writeNeuronSwc'));
jsonText = fileread(fullfile(cFolder,'settings.json'));
settings = jsondecode(jsonText);

%% Parse input.
p = inputParser;
p.addRequired('locStr',@(x) ischar(x) | iscell(x));
p.addOptional('Dilation',0,@(x) isnumeric(x) & length(x)==1);
p.addParameter('UrlTracings',settings.Database.TracingsUrl,@(x) ischar(x));
p.addParameter('MeshFile',fullfile('//nrs/mouselight/Shared Files/Mesh Info/allenMesh.mat'),@(x) ischar(x));
p.parse(locStr,varargin{:});
Inputs = p.Results;

if ischar(Inputs.locStr), Inputs.locStr = {Inputs.locStr}; end

%% Load mesh info.
fprintf('\nLoading mesh file..');
load(Inputs.MeshFile);

%% get firstNode info.
query ='{ tracings(queryInput:{}){ tracings{ firstNode{ x y z } tracingStructure{ name } swcTracing { neuron{ idString } } } } }';
[ response ] = callgraphql( Inputs.UrlTracings, query);
neuronId = [response.tracings.tracings.swcTracing];
neuronId = [neuronId.neuron];
neuronId = {neuronId.idString};
somaLoc = [response.tracings.tracings.firstNode];
somaLoc = [[somaLoc.x]' [somaLoc.y]' [somaLoc.z]'];
% select only axon.
type = [response.tracings.tracings.tracingStructure];
type = {type.name};
indAxon = cellfun(@(x) strcmpi(x,'axon'),type);
somaLoc = somaLoc(indAxon,:);
neuronId = neuronId(indAxon);

%% Make voxel search area.
[ ontIm ] = VoxelizedBrainArea( [100,100,100], Inputs.locStr, allenMesh );
if Inputs.Dilation>0
    se = strel('sphere',Inputs.Dilation);
    ontIm = imdilate(ontIm,se);
end

%% Check Soma Locations.
pixLoc = ceil(somaLoc./100);
indPos = sub2ind(size(ontIm),pixLoc(:,1),pixLoc(:,2),pixLoc(:,3));
neuronId = neuronId(ontIm(indPos));

end



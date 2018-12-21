function [neuronStr,areaName] = queryNeuronLoc( locStr, url )
%% Parameters.
if nargin==1
    [cFolder,~,~] = fileparts(which('getNeuronfromIdString'));
    jsonText = fileread(fullfile(cFolder,'settings.json'));
    settings = jsondecode(jsonText); 
    url = settings.Database.SampleUrl;
end

%% brain area to id.
query = '{ brainAreas{ name structureId } }';
[ areaResponse ] = callgraphql( url, query);
areas = {areaResponse.brainAreas.name};
ind = find(cellfun(@(x) strcmpi(x,locStr),areas));
if isempty(ind)
    error('Could not find location: %s',locStr);
end
id = areaResponse.brainAreas(ind).structureId;

%% get neuron ID list.
query = '{ neurons { items{ idString brainArea{ safeName structureIdPath } } } }';
[ neuronResponse ] = callgraphql( url, query);
somaLocs = [neuronResponse.neurons.items.brainArea];
somaLocs = {somaLocs.structureIdPath};
ind = find(cellfun(@(x) ~isempty(regexp(x,sprintf('/%i/',id),'start')),somaLocs));

%% 
neuronStr = {neuronResponse.neurons.items(ind).idString}';
%@12/25/18: fixed missing entity
areaName = {neuronResponse.neurons.items(ind).brainArea}';
% areaName = {areaName.safeName}';

end


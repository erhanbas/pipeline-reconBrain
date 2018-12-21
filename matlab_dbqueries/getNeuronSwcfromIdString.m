function [ neuron ] = getNeuronSwcfromIdString( idString, tracingDbUrl, swcDbUrl )
if nargin~=3
    [cFolder,~,~] = fileparts(which('getNeuronfromIdString'));
    jsonText = fileread(fullfile(cFolder,'settings.json'));
    settings = jsondecode(jsonText);
end
if nargin<2
    tracingDbUrl = settings.Database.TracingsUrl;
end
if nargin<3 
    swcDbUrl = settings.Database.SWCUrl;
end
query = '{ swcTracings{ id tracingStructure{ name } neuron{ idString } } }';
[ response ] = callgraphql( tracingDbUrl, query);
id = [response.swcTracings(:).neuron];
ind = find(cellfun(@(x) strcmp(x,idString), {id.idString}));
selected = response.swcTracings(ind);
% Find axon/dendrite.
structNames = [selected.tracingStructure];
structNames = {structNames.name};
axonInd = find(strcmp(structNames,'axon'));
dendInd = find(strcmp(structNames,'dendrite'));
% Store axon info.
if isempty(axonInd)
    warning('Found no axon tracing for: %s',idString);
    neuron.axon = [];
else
    if length(axonInd)>1
        warning('Found multiple axon tracings for %s',idString);
    end
    neuron.axon = getSwcfromId( selected(axonInd(1)).id,swcDbUrl);
end
% Store dendrite.
if isempty(dendInd)
    warning('Found no dendrite tracing for: %s',idString);
    neuron.dendrite = [];
else
    if length(dendInd)>1
        warning('Found multiple axon tracings for %s',idString);
    end
    neuron.dendrite = getSwcfromId( selected(dendInd(1)).id,swcDbUrl);
end

end


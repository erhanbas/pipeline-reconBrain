function [ data ] = getSomafromIdString( idString,varargin)
%% Read settings file.
[cFolder,~,~] = fileparts(which('getSomafromIdString'));
jsonText = fileread(fullfile(cFolder,'settings.json'));
settings = jsondecode(jsonText);
    
%% Parse input.
p = inputParser;
p.addRequired('idString',@(x) ischar(x) && length(x)==6);
p.addParameter('Url',settings.Database.TracingsUrl,@(x) ischar(x));
p.addParameter('ForceHemi','no',@(x) ischar(x));
p.parse(idString,varargin{:});
Inputs = p.Results;

% Parameters.
halfPoint = 5695;
%getNeuronfromIdString
query = '{ swcTracings{ id tracingStructure{ name } neuron{ idString } } }';
[ response ] = callgraphql( Inputs.Url, query);
id = [response.swcTracings(:).neuron];
ind = find(cellfun(@(x) strcmp(x,idString), {id.idString}));
if isempty(ind), error('Could not find neuron: %s in database',idString); end
selected = response.swcTracings(ind);
% Find axon
structNames = [selected.tracingStructure];
structNames = {structNames.name};
axonInd = find(strcmp(structNames,'axon'));

% get firstNode info.
query = sprintf('{ tracings(queryInput: {swcTracingIds: "%s"}) { tracings { firstNode { sampleNumber x y z parentNumber structureIdValue brainArea { structureId atlasId safeName acronym structureIdPath } } } } }',...
    selected(axonInd(1)).id);
[ response ] = callgraphql( Inputs.Url, query);
data = response.tracings.tracings.firstNode;

%% Clean up output to single structure.
fn2 = fieldnames(data(1).brainArea);
ana = squeeze(struct2cell([data.brainArea]));
data = rmfield(data,'brainArea');
fn1 = fieldnames(data(1));
fn = [fn1',fn2'];
data = struct2cell(data);
data = [data',ana'];
data = cell2struct(data',fn);

% Mirror if requested
if (strcmpi(Inputs.ForceHemi,'Left') && data.x<halfPoint)
    data.x = -data.x+halfPoint*2;
end
if (strcmpi(Inputs.ForceHemi,'Right') && data.x>halfPoint)
    data.x = data.x-halfPoint;
end

end


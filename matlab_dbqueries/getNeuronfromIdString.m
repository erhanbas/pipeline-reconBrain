function [ neuron ] = getNeuronfromIdString( idString,varargin)
%% Read settings file.
[cFolder,~,~] = fileparts(which('writeNeuronSwc'));
jsonText = fileread(fullfile(cFolder,'settings.json'));
settings = jsondecode(jsonText);
    
%% Parse input.
p = inputParser;
p.addRequired('idString',@(x) ischar(x) && length(x)==6);
p.addParameter('Url',settings.Database.TracingsUrl,@(x) ischar(x));
p.addParameter('ForceHemi','no',@(x) ischar(x));
p.addParameter('Type',{'axon','dendrite'},@(x) iscell(x) || ischar(x));
p.addParameter('Registered',true,@(x) islogical(x));
p.parse(idString,varargin{:});
Inputs = p.Results;
if ischar(Inputs.Type), Inputs.Type={Inputs.Type}; end
% Parameters.
halfPoint = 5695;
%getNeuronfromIdString
query = '{ swcTracings{ id tracingStructure{ name } neuron{ idString } } }';
[ response ] = callgraphql( Inputs.Url, query);
id = [response.swcTracings(:).neuron];
ind = find(cellfun(@(x) strcmp(x,idString), {id.idString}));
if isempty(ind), error('Could not find neuron: %s in database',idString); end
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
elseif any(strcmpi(Inputs.Type,'axon'))==0
    neuron.axon = [];
else
    if length(axonInd)>1
        warning('Found multiple axon tracings for %s',idString);
    end
    neuron.axon = getTracingfromId( selected(axonInd(1)).id,Inputs.Url);
            % Force hemisphere.
        if (strcmpi(Inputs.ForceHemi,'Left') && neuron.axon(1).x<halfPoint)...
                || ( strcmpi(Inputs.ForceHemi,'Right') && neuron.axon(1).x>halfPoint)
            tMat = eye(4,4);
            tMat(1,1) = -1;
            for i = 1:size(neuron.axon,1)
                newCoord = tMat*[neuron.axon(i).x,neuron.axon(i).y,neuron.axon(i).z,0]';
                neuron.axon(i).x = newCoord(1)+halfPoint*2;
                neuron.axon(i).y = newCoord(2);
                neuron.axon(i).z = newCoord(3);
            end
        end
end
% Store dendrite.
if isempty(dendInd)
    warning('Found no dendrite tracing for: %s',idString);
    neuron.dendrite = [];
elseif any(strcmpi(Inputs.Type,'dendrite'))==0
    neuron.dendrite = [];
else
    if length(dendInd)>1
        warning('Found multiple axon tracings for %s',idString);
    end
    neuron.dendrite = getTracingfromId( selected(dendInd(1)).id,Inputs.Url);
             % Force hemisphere.
        if (strcmpi(Inputs.ForceHemi,'Left') && neuron.dendrite(1).x<halfPoint)...
                || ( strcmpi(Inputs.ForceHemi,'Right') && neuron.dendrite(1).x>halfPoint)
            tMat = eye(4,4);
            tMat(1,1) = -1;
            for i = 1:size(neuron.dendrite,1)
                newCoord = tMat*[neuron.dendrite(i).x,neuron.dendrite(i).y,neuron.dendrite(i).z,0]';
                neuron.dendrite(i).x = newCoord(1)+halfPoint*2;
                neuron.dendrite(i).y = newCoord(2);
                neuron.dendrite(i).z = newCoord(3);
            end
        end
end
end


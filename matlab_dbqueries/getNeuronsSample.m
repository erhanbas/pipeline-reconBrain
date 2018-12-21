function [ cellNames, coords ] = getNeuronsSample( sampleStr, varargin )
%% Read settings file.
[cFolder,~,~] = fileparts(which('getNeuronsSample'));
jsonText = fileread(fullfile(cFolder,'settings.json'));
settings = jsondecode(jsonText);

%% Parse input.
p = inputParser;
p.addRequired('sampleStr',@(x) ischar(x) && length(x)==10);
p.addParameter('SWCUrl',settings.Database.SWCUrl,@(x) ischar(x));
p.parse(sampleStr,varargin{:});
Inputs = p.Results;

%% Get samples.
query = '{ samples{ id sampleDate } }';
response = callgraphql(Inputs.SWCUrl,query);  
samples = [response.samples.sampleDate]./1000;
samples = datetime(samples,'ConvertFrom', 'posixtime');
samples = cellstr(datestr(samples,'yyyy-mm-dd'));
ind = find(cellfun( @(x) strcmp(x,sampleStr),samples));
sampleId = response.samples(ind).id;

%% Get neurons.
query = sprintf('{ neurons(sampleId:"%s"){ id idString } }',sampleId);
cells = callgraphql(Inputs.SWCUrl,query);  
cellNames = {cells.neurons.idString}';

%% get soma coords.
coords = [];
for iNeuron=1:size(cellNames,1)
    % get tracings.
    query = sprintf('{ tracings(pageInput:{ neuronIds: "%s"}){ tracings{ id tracingStructure { name } } } }',cells.neurons(iNeuron).id);
    response = callgraphql(Inputs.SWCUrl,query);  
    type = [response.tracings.tracings.tracingStructure];
    type = {type.name};
    ind = find(cellfun(@(x) strcmp(x,'axon'),type));
    swcId = response.tracings.tracings(ind).id;
    % % get swc.
    swc=getSwcfromId(swcId,Inputs.SWCUrl);
    coords = [coords;swc(1).x,swc(1).y,swc(1).z];
end
end
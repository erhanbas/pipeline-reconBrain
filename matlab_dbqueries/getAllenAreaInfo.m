function info = getAllenAreaInfo(structureId,varargin)
% structureId can be an array.
%% Read settings file.
[cFolder,~,~] = fileparts(which('getBrainAreaInfo'));
jsonText = fileread(fullfile(cFolder,'settings.json'));
settings = jsondecode(jsonText);

%% Parse input.
p = inputParser;
p.addRequired('structureId',@(x) isnumeric(x));
p.addParameter('Url',settings.Database.SampleUrl,@(x) ischar(x));
p.parse(structureId,varargin{:});
Inputs = p.Results;

%% Call database.
query = sprintf('{ brainAreas(input: {sortField: "structureId"}){ structureId structureIdPath acronym name safeName geometryColor atlasId graphOrder} }');
[ data ] = callgraphql( Inputs.Url, query);

%% Find info.
[~,ind] = ismember(structureId,[data.brainAreas.structureId]);
if sum(ind==0)>0
    warning('Could not find location info of certain nodes (structureIDs not found)');
    ind(ind==0)=987; %brain.
end
info = data.brainAreas(ind);

end
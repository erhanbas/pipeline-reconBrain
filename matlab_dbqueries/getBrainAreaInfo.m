function getBrainAreaInfo(structureId,varargin)
%% Read settings file.
[cFolder,~,~] = fileparts(which('getBrainAreaInfo'));
jsonText = fileread(fullfile(cFolder,'settings.json'));
settings = jsondecode(jsonText);

%% Parse input.
p = inputParser;
p.addRequired('structureId',@(x) ischar(x) && length(x)==6);
p.addParameter('Url',settings.Database.TracingsUrl,@(x) ischar(x));
p.parse(idString,varargin{:});
Inputs = p.Results;

%% Call database.
query = sprintf('{ tracings(queryInput: {swcTracingIds: "%s"}) { tracings { nodes { sampleNumber x y z parentNumber structureIdValue brainArea { structureId atlasId safeName acronym structureIdPath hexColor } } } } }',...
    id);
[ data ] = callgraphql( url, query);
end
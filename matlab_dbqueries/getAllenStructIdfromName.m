function structIds = getAllenStructIdfromName(names,varargin)
%% Read settings file.
[cFolder,~,~] = fileparts(which('getAllenStructIdfromName'));
jsonText = fileread(fullfile(cFolder,'settings.json'));
settings = jsondecode(jsonText);

%% Parse input.
p = inputParser;
p.addRequired('names',@(x) ischar(x) || iscell(x));
p.addParameter('Url',settings.Database.TracingsUrl,@(x) ischar(x));
p.parse(names,varargin{:});
Inputs = p.Results;
if ischar(Inputs.names), Inputs.names={Inputs.names}; end

%% Call database.
query = sprintf('{ brainAreas{ structureId name structureIdPath} }');
[ data ] = callgraphql( Inputs.Url, query);
% find name.
names = {data.brainAreas.name};
%% Lookup areas
structIds =[];
for iName = 1:size(Inputs.names,2)
    cName = Inputs.names{iName};
    ind = find(strcmpi(names,cName));
    if isempty(ind)
        warning('Could not find location %s in database',cName);
        ind = 0;
    end
    structIds = [structIds;data.brainAreas(ind).structureId];
end
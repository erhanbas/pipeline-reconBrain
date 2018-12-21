function writeNeuronSwc( idString, folder, varargin )
%% Read settings file.
[cFolder,~,~] = fileparts(which('writeNeuronSwc'));
jsonText = fileread(fullfile(cFolder,'settings.json'));
settings = jsondecode(jsonText);
    
%% Parse input.
p = inputParser;
p.addRequired('idString',@(x) ischar(x) && length(x)==6);
p.addRequired('folder',@(x) ischar(x));
p.addOptional('type',{'axon','dendrite'},@(x) (iscell(x) || ischar(x)) && any(strcmp(x,{'axon','dendrite'})));
p.addParameter('Color',rand(1,3),@(x) min(x)>=0 && max(x)<=1 && isequal(size(x),[1,3]));
p.addParameter('Radius',1,@(x) isnumeric(x));
p.addParameter('Url',settings.Database.TracingsUrl,@(x) ischar(x));
p.addParameter('UrlSwc',settings.Database.SWCUrl,@(x) ischar(x));
p.addParameter('ForceHemi','no',@(x) ischar(x));
p.addParameter('Registered',true,@(x) islogical(x));

p.parse(idString,folder,varargin{:});
Inputs = p.Results;
if ischar(Inputs.type), Inputs.type = {Inputs.type}; end

%% Check if output folder exist otherwise crease
if ~isdir(Inputs.folder), mkdir(Inputs.folder); end

%% Get neuron from database
if Inputs.Registered
    neuron = getNeuronfromIdString(Inputs.idString,'ForceHemi',Inputs.ForceHemi,'Url', Inputs.Url);
else
    neuron = getNeuronSwcfromIdString(Inputs.idString,Inputs.Url,Inputs.UrlSwc);
end

%% Write Neurite Type.
for field = Inputs.type
    if isempty(neuron.(field{:}))
        warning('%s not present for %s',field{:},Inputs.idString);
    else
        % Assemble swc.
        cNeuron = neuron.(field{:});
        swc = [[cNeuron.sampleNumber]',[cNeuron.structureIdValue]',...
            [cNeuron.x]',[cNeuron.y]',[cNeuron.z]',...
            repmat(Inputs.Radius,size(cNeuron,1),1), [cNeuron.parentNumber]'];
        % Write swc.
        outputFile = fullfile(Inputs.folder,sprintf('%s_%s.swc',Inputs.idString,field{:}));
        fid = fopen(outputFile,'w');
        % Header.
        fprintf(fid,'# ORIGINAL_SOURCE MouseLight Database');
        fprintf(fid,'\n# OFFSET 0 0 0');
        fprintf(fid,'\n# NEURON %s %s',Inputs.Color);
        fprintf(fid,'\n# GENERATED ON %s',datestr(now,'yy/mm/dd HH:MM'));
        fprintf(fid,'\n%i %i %.6f  %.6f  %.6f  %.6f %i',swc');
        fclose(fid);
    end
end


end


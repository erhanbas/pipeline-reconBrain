function [ data ] = getTracingfromId( id,url)
%getTracingfromId

%% Call database.
query = sprintf('{ tracings(queryInput: {swcTracingIds: "%s"}) { tracings { nodes { sampleNumber x y z parentNumber structureIdValue brainArea { structureId atlasId safeName acronym structureIdPath } } } } }',...
    id);
[ data ] = callgraphql( url, query);
%% Sort according to node ID.
data = data.tracings.tracings.nodes;
[~,ind] = sort([data.sampleNumber]','ascend');
data = data(ind);
%% check for empty anatomy fields.
empty_arr = cellfun(@(x) isempty(x), {data.brainArea});
sample_id = find(~empty_arr,1);
emptyInd = find(empty_arr);
%% Clean up output to single structure.
fn2 = fieldnames(data(sample_id).brainArea);
for i=emptyInd
data(i).brainArea = struct('structureId',[],'atlasId',[],'safeName',[],'acronym',[],'structureIdPath',[]); 
end
ana = squeeze(struct2cell([data.brainArea]));
data = rmfield(data,'brainArea');
fn1 = fieldnames(data(1));
fn = [fn1',fn2'];
data = struct2cell(data);
data = [data',ana'];
data = cell2struct(data',fn);

end


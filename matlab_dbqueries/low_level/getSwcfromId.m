function [ data ] = getSwcfromId( id,url)
%getTracingfromId

%% Call database.
query = sprintf('{ tracingNodes(id:"%s"){ sampleNumber structureIdValue x y z radius parentNumber } }',...
    id);
[ data ] = callgraphql( url, query);
%% Sort according to node ID.
data = data.tracingNodes;
[~,ind] = sort([data.sampleNumber]','ascend');
data = data(ind);

%% Get offset.
query = sprintf('{ tracing(id: "%s") { offsetX offsetY offsetZ } }',...
    id);
[ offset ] = callgraphql( url, query);
offset = offset.tracing;

%% Add offset.
% This is ugly and I dotn care.
for i =1:size(data,1)
    data(i).x = data(i).x + offset.offsetX;
    data(i).y = data(i).y + offset.offsetY;
    data(i).z = data(i).z + offset.offsetZ;
end
end

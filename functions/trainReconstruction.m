function trainReconstruction(opt)
% loads reconstructions from repo and trains a desicion classifier
%%
neurons = repo.crawler([]); % will return all axons in the repo.
save(fullfile(tempfold,sprintf('neurons-%s.mat',datestr(now,'yymmdd'))),'neurons')

numneurons = length(neurons);
[neuron_list,areaName] = deal([]);
for in = 1:numneurons
    neuron_list{in} = [neurons{in}.name(1:6),'-',neurons{in}.acronym];
    areaName{in}.safeName = neurons{in}.soma.safeName;
    areaName{in}.structureIdPath = neurons{in}.soma.structureIdPath;
end
[allen_neuron_color,structureId] = funcs.areaName2Color(areaName);

% envelope
sliceRange = [6500,7000]; dimSelection = [1 2];
[ env, color, varargout] = funcs.getMask( region, dimSelection, sliceRange );

rundate = datestr(now,'yymmdd');

end
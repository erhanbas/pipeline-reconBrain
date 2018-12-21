function feats = featsFromRecons(neurons_repo,atlas)
%create feature for groundtruth reconstruction
% neurons_repo = neuronrepo.neurons;
numGTneuron = length(neurons_repo);
dims_atlas = size(atlas);
feats = cell(1,numGTneuron);

try parfor_progress(0);catch;end
parfor_progress(numGTneuron)

parfor ineuron = 1:numGTneuron
    %%
    neuro_in = neurons_repo{ineuron};
    Gneuron = graph(max(neuro_in.recon.A',neuro_in.recon.A)); %#ok<UDIM>
    subsneuron = neuro_in.recon.subs;
    
    [branches_neuron,nodeBrid_neuron,branch_pair_distance_neuron,branch_pair_connectivity_neuron,branch_pair_dissimilarity_neuron] = ...
        feats4recon(Gneuron,subsneuron);
    %%
    % get sub location for junction
    [ia,ib,ic] = find(branch_pair_connectivity_neuron);
    junctionlocation = zeros(length(ia),3);
    for ix = 1:length(ia)
        pd = pdist2(branches_neuron(ia(ix)).tipsubs, branches_neuron(ib(ix)).tipsubs);
        [minval,imy] = min(pd(:));
        [imx,imy] = ind2sub([2,2],imy);
        junctionlocation(ix,:) = (branches_neuron(ia(ix)).tipsubs(imx,:) + branches_neuron(ib(ix)).tipsubs(imy,:))/2;
    end
    
    subsneuron_ = ceil(junctionlocation/10);
    
    indPoints_junctions = sub2ind_withnan(dims_atlas,subsneuron_);
    % indPoints_junctions = sub2ind(dims_atlas,subsneuron_(:,1),subsneuron_(:,2),subsneuron_(:,3));
    inds_junctions = indPoints_junctions;
    inds_junctions(isfinite(indPoints_junctions)) = atlas(indPoints_junctions(isfinite(indPoints_junctions)));
    
    % feat set for junctions [allen compartment, D]
    % D(1): Euclidean
    % D(2): theta
    % D(3): PCA
    % D(4): KL
    % D(5 [end]): hit
    vals = zeros(nnz(branch_pair_connectivity_neuron),1);
    for ik=1:size(branch_pair_dissimilarity_neuron,2)
        vals(:,ik) = branch_pair_dissimilarity_neuron{ik}(branch_pair_connectivity_neuron);
    end
    feats{ineuron} = [inds_junctions(:), vals];
    parfor_progress;
end
parfor_progress(0)


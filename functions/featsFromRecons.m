function feats = featsFromRecons(branches,branch_pair_dissimilarity,branch_pair_connectivity,atlas)
dims_atlas = size(atlas);
[ia,ib,ic] = find(branch_pair_connectivity);
junctionlocation = zeros(length(ia),3);
parfor ix = 1:length(ia)
    pd = pdist2(branches(ia(ix)).tipsubs, branches(ib(ix)).tipsubs);
    [minval,imy] = min(pd(:));
    [imx,imy] = ind2sub([2,2],imy);
    junctionlocation(ix,:) = (branches(ia(ix)).tipsubs(imx,:) + branches(ib(ix)).tipsubs(imy,:))/2;
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
for ik=1:size(branch_pair_dissimilarity,2)
    vals(:,ik) = branch_pair_dissimilarity{ik}(branch_pair_connectivity_neuron);
end
feats = [inds_junctions(:), vals];

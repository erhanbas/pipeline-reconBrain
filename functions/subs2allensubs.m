function [subsallen,ids] = subs2allensubs(subs,params,vectorField,tFormVecField,allentForm,atlas)

%%
dims_atlas = size(atlas);

%% voxel -> um -> allen
A_vox2um = [[diag([params.voxres]) [params.ox params.oy params.oz]'/1e3];[0 0 0 1]];
subs_um = [subs-1 ones(size(subs,1),1)]*A_vox2um'; % added -1 to match to neuronbrowser transformation. not sure which one is correct....
subs_vecfield = subs_um*tFormVecField';

%%
pixPosVecField = ceil(subs_vecfield);
nPoints = size(pixPosVecField,1);

% ravel multiindex
dims = size(vectorField);
inds = sub2ind(dims(1:3),pixPosVecField(:,1),pixPosVecField(:,2),pixPosVecField(:,3));
rvectorField = reshape(vectorField,[],3);
vector = rvectorField(inds,:);
subs_um_updated = ceil(subs_um(:,1:3) + vector);
pixPos = [subs_um_updated,ones(nPoints,1)]*allentForm';
pixPos_ceil = ceil(pixPos);

%%
ids = nan(size(subs,1),1);
subsallen = pixPos(:,1:3);
val_outofbounds = (pixPos_ceil(:,1)>dims_atlas(:,1) | pixPos_ceil(:,2)>dims_atlas(:,2) | pixPos_ceil(:,3)>dims_atlas(:,3) | any(pixPos_ceil<1,2));
subsallen(val_outofbounds,:) = nan;

inds_inofbounds = find(~val_outofbounds);
% subsallen(inds_inofbounds,:) = pixPos(inds_inofbounds,:);

%%
indPoints = sub2ind(dims_atlas,pixPos_ceil(inds_inofbounds,1),pixPos_ceil(inds_inofbounds,2),pixPos_ceil(inds_inofbounds,3));
ids(inds_inofbounds) = atlas(indPoints);
end

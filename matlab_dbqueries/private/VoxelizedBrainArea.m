function [ ontIm ] = VoxelizedBrainArea( voxelSize, locStr, allenMesh )
names = {allenMesh.name};

%% Creat result variable.
ontIm = zeros(round(11400/voxelSize(1)),round(8000/voxelSize(2)),round(13200/voxelSize(3)),'logical');
for iAna = 1:size(locStr,2)
    %% find correct mesh.
    meshId = find(strcmpi(names,locStr{iAna}));
    if isempty(meshId), error('Could not area: %s',locStr{iAna}); end
    %% Prepare Result image.
    gridX = voxelSize(1)/2:voxelSize(1):(1140*10);
    gridY = voxelSize(2)/2:voxelSize(2):(800*10);
    gridZ = voxelSize(3)/2:voxelSize(3):(1320*10);
    meshFV = struct('faces',allenMesh(meshId).f,'vertices',allenMesh(meshId).v);
    BW = VOXELISE(gridX,gridY,gridZ,meshFV) ;
    ontIm(BW) = true;
end

end


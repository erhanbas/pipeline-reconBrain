function [ids,subsallen,vectorField,tFormVecField,allentForm,atlas] = getAllenAtlasIds(params,subs,transformFile,allenH5)

if nargin<3
    allenH5 = '/groups/mousebrainmicro/mousebrainmicro/Software/Matlab/Registration tools/Dependencies/OntologyAtlas.h5';
    transformFile = '/groups/mousebrainmicro/mousebrainmicro/registration/Database/2018-08-01/Transform.2018-08-01.h5';
end
try
    tic
    vectorField = h5read(transformFile,'/DisplacementField');
    toc
    tFormVecField = h5readatt(transformFile,'/DisplacementField','Transformation_Matrix');
catch
    error('Could not read Displacement Field file: %s',transformFile);
end

allentForm = h5readatt(allenH5,'/OntologyAtlas','Transformation_Matrix');
atlas = h5read(allenH5,'/OntologyAtlas');

[subsallen,ids] = subs2allensubs(subs,params,vectorField,tFormVecField,allentForm,atlas);

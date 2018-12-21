function [outputArgs] = reconBrain(configfile)
%reconstruction of whole mouse brain

% $Author: base $	$Date: 2018/12/15 00:47:01 $
% Copyright: HHMI 2016
if nargin<1
    configfile = './config_files/config_reconBrain_20180801_prob0.cfg';
end
%%

addpath(genpath('./common'))
addpath(genpath('./functions'))
opt = configparser(configfile);
if ~isfield(opt,'sampling')
    opt.sampling = 'uni';
end

[aa,bb,cc] = fileparts(opt.inputh5);
tempfold = fullfile(pwd,'temp',opt.tag,bb);
mkdir(tempfold)

myh5 = opt.inputh5;
myh5prob = opt.h5prob;
[~,name] = fileparts(myh5);
h5infofile = fullfile('./h5infos',['h5inf_',name,'.mat']);
mkdir(fileparts(h5infofile))

[brainSize,RR,chunk_dims,rank] = h5parser(myh5,myh5prob);
%%
origin = h5read(opt.inputh5,[opt.h5prob,'_props/origin']);
spacing = h5read(opt.inputh5,[opt.h5prob,'_props/spacing']);
level = h5read(opt.inputh5,[opt.h5prob,'_props/level']);
params.outsiz = brainSize;
params.ox = origin(1);
params.oy = origin(2);
params.oz = origin(3);
params.sx = spacing(1);
params.sy = spacing(2);
params.sz = spacing(3);
params.level = level;

params.voxres = [params.sx params.sy params.sz]/2^(params.level)/1e3; % in um
opt.params = params;
%%
% create output folders
mkdir(fullfile(opt.outfolder,'full'))
mkdir(fullfile(opt.outfolder,'frags'))
%%
[subs,edges,A,weights] = skel2graph(opt);
subs_ori = subs;
edges_ori = edges;
weights_ori = weights;
A_ori = A;
%% apply allen mask if provided
% 0: outside
% 507: olfactory
% 
maskids = [0 507 73 81 89 98 108 116 124 129 140 145 153 164]; % outside and ventrical system

% maskids = [56 672 1022 1031]; %ACB/Caudoputamen/Globus pallidus, external segment/Globus pallidus, internal segment
%maskids = [956 844 882 686 56 1022 1031 1021 1085 719 882 583 182305705 182305709 182305713]
% gt_swcfolder = '/nrs/mouselight/seggui/swcfiles/GT/2017-09-25_striatum_neurons_temp'
if ~isempty(maskids)
    %%
    transformFile = '/groups/mousebrainmicro/mousebrainmicro/registration/Database/2018-08-01/Transform.2018-08-01.h5';

    [hits_allen_brain] = maskWithAllenAtlas(params,subs,maskids,1);
    if exist('gt_swcfolder','var')
        addpath(genpath('./scripts'))
        [hits_gt,hits_delete,swcout,gtfile] = cropSectionBasedOnGT(params,gt_swcfolder,subs);
        % keep allen_brain and gt then substract delete
        keepthese = union(setdiff(hits_allen_brain,hits_delete),hits_gt);
        figure, myplot3(subs(keepthese,:),'.')
    else
        keepthese = hits_allen_brain;
    end
    %%
    subs = subs(keepthese,:);
    A = A(keepthese,:);
    A = A(:,keepthese);
end


%%
subs = subs_ori;
edges = edges_ori;
weights = weights_ori;
A = A_ori;
Gin = graph(max(A,A'));
%% workflow 1, fast but not accurate
if 0
    workflow1(Gin,subs,opt)
end
%% TODO
workflow3(Gin,subs,opt)
%%
%
% %
% tstart = tic;
% affinityBuilder(opt,A,subs)
% sprintf('FINISHED IN: %d', round(toc(tstart)))

end

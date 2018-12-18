function [outputArgs] = reconBrain(configfile)
%reconstruction of whole mouse brain

% $Author: base $	$Date: 2018/12/15 00:47:01 $
% Copyright: HHMI 2016

if nargin<1
    configfile = './config_files/config_reconBrain_20180801_prob0.cfg';
end

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



%% PREPROCESSING: prune & filter graph
tr=tic;
if opt.prune == 1
    size_pruning_Thr = [100 50];
    [Gpruned,subspruned] = pruneGraph(Gin,subs,size_pruning_Thr,[]);
    save(fullfile(tempfold,'Gpruned'),'Gpruned','subspruned')
elseif opt.prune == 2
    load(fullfile(tempfold,'Gpruned'),'Gpruned','subspruned')
else
end

%% CONVERSION to branch
G = Gpruned;
subs = subspruned;
save(fullfile(tempfold,'G_and_subs.mat'),'G','subs','-v7.3')

if opt.graph2branch == 1
    % [branches,nodeBrid] = graph2branch(Gfilt,subsfilt);
    [branches,nodeBrid] = graphfuncs.graph2branch(G,subs);
    save(fullfile(tempfold,['branchlist',tag]),'branches','nodeBrid')
elseif opt.graph2branch == 2
    load(fullfile(tempfold,['branchlist',tag]),'branches','nodeBrid')
else
end
sprintf('TOTAL TIME: %d secs',round(toc(tr)))

%% TODO
workflow3(G,subs,opt)
%%
%
% %
% tstart = tic;
% affinityBuilder(opt,A,subs)
% sprintf('FINISHED IN: %d', round(toc(tstart)))

end
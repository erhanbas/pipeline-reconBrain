function workflow3(Gin,subs,opt)

addpath(genpath('functions'))
addpath(genpath('common'))
tr=tic;


%% PREPROCESSING: prune & filter graph
if opt.prune == 1
    if isfield(opt,'size_pruning_Thr')
        size_pruning_Thr = opt.size_pruning_Thr;
    else
        size_pruning_Thr = [100 50];
    end
        
    [Gpruned,subspruned] = pruneGraph(Gin,subs,size_pruning_Thr,[]);
    save(fullfile(tempfold,'Gpruned'),'Gpruned','subspruned')
elseif opt.prune == 2
    load(fullfile(tempfold,'Gpruned'),'Gpruned','subspruned')
else
end
sprintf('TOTAL TIME: %d secs',round(toc(tr)))

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

%%
if opt.branch2conn == 1
    if isfield(opt,'querdist')
        querdist = opt.querdist; % masking distance to calculate pair scores
    else
        querdist = [50];
    end
    tr=tic;
    branch_pair_distance = graphfuncs.branchConn(branches,subs,nodeBrid,querdist); % pwdist of branches
    sprintf('branchConn done in: %d secs',round(toc(tr)))
    [aa,bb,cc] = find(branch_pair_distance);
    tr=tic;
    branch_pair_connectivity = sparse(aa,bb,cc<=querdist,size(branch_pair_distance,1),size(branch_pair_distance,1)); %pwcost of branches
    branch_pair_dissimilarity = graphfuncs.calcDists(branches,branch_pair_connectivity);
    sprintf('calcDists done in: %d secs',round(toc(tr)))
    save(fullfile(tempfold,['branchConn',tag]),'branch_pair_distance','branch_pair_dissimilarity','opt')
elseif  opt.branch2conn == 2
else
end

%% apply allen transformation
allenH5 = '/groups/mousebrainmicro/mousebrainmicro/Software/Matlab/Registration tools/Dependencies/OntologyAtlas.h5';
transformFile = '/groups/mousebrainmicro/mousebrainmicro/registration/Database/2018-08-01/Transform.2018-08-01.h5';
[ids,subsallen,vectorField,tFormVecField,allentForm,atlas] = getAllenAtlasIds(opt.params,subs,transformFile,allenH5);
[subsallen,ids] = subs2allensubs(subs,opt.params,vectorField,tFormVecField,allentForm,atlas);

%% sample feature set
branch_pair_dissimilarity

%% train classifier based on repo
neuronrepo = load('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/pipeline-reconBrain/temp/wrapJRC/20180801_prob0_lev-6_chunk-111_111_masked-0/neurons-181220.mat')
feats = featsFromGTRecons(neuronrepo.neurons,atlas);
feats_juncs = featsFromRecons(branches,branch_pair_dissimilarity,branch_pair_connectivity,atlas)

% feat set for junctions [allen compartment, D]
% D(1): Euclidean
% D(2): theta
% D(3): PCA
% D(4): KL
% D(5 [end]): hit
% save(fullfile(tempfold,['GT_feats_qd_1']),'feats')

%%
feats_array = cat(1,feats{1:917}); % upto 2018-08-01 sample
feats_array = feats_array(isfinite(feats_array(:,1)),:);
compartments = unique(feats_array(:,1))';
hists = histcounts(feats_array(:,1),compartments);

%%
for icomp = compartments
    % skip;
    if icomp==0;continue;end
    inds_in_compartment = feats_array(:,1)==icomp;
    feat_subset = feats_array(inds_in_compartment,:);
    
    %% 
    inds_in_sub = ids==icomp;
    feat_sample = 
    
    %%
    figure(13), 
    subplot(3,1,1); 
    histogram(feat_subset(feat_subset(:,3)<10,3))
    subplot(3,1,2); 
    histogram(feat_subset(:,4))
    subplot(3,1,3); 
    histogram(feat_subset(:,5))
    
    
end
%% filter out based on Allen mask and context (edge density or alignment). 
% Implements a binary classifier to decide if an edge should be kept or
% deleted as a function of Allen compartment ids, and
% branch_pair_dissimilarity 



function workflow3(Gin,subs,opt)

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
%% filter out based on Allen mask and context (edge density or alignment). 
% Implements a binary classifier to decide if an edge should be kept or
% deleted as a function of Allen compartment ids, and
% branch_pair_dissimilarity measures.





function [Gfilt,subsfilt] = pruneGraph(Gin,subsin,size_pruning_Thr,opt)
%PRUNEGRAPH Summary of this function goes here
%
% [OUTPUTARGS] = PRUNEGRAPH(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2017/03/29 16:19:40 $	$Revision: 0.1 $
% Copyright: HHMI 2017
% distThr = 10;
% sizeThr = 10;
if issparse(Gin)
    Gin = graph(max(Gin,Gin'));
elseif class(Gin)=='graph'
else
    error('check Ain. it can be sparse matrix or graph structure')
end
if length(size_pruning_Thr)==1
    sizeThr=size_pruning_Thr(1);
    pruningThr=size_pruning_Thr(1);
elseif length(size_pruning_Thr)==2
    sizeThr=size_pruning_Thr(1);
    pruningThr=size_pruning_Thr(2);
else
    error('sizeThr can be a scalar (sizeThr==pruningThr) or 1x2 vertor [sizeThr, pruningThr]')
end

tic_prune = tic;
filt{1}{1} = 'graphfuncs.filtGsize';
filt{1}{2} = {sizeThr};
[Gfilt,subsfilt] = graphfuncs.filterGraph(Gin,subsin,filt);
sprintf('FILTER DONE in %d sec',round(toc))

tic
[Gfilt,subsfilt] = graphfuncs.deleteLoops(Gfilt,subsfilt);
sprintf('del loop DONE in %d sec',round(toc))

tic
iter=0;
G = Gfilt;
subs=subsfilt;
pruned = 1;
while pruned
    iter=iter+1;
    % [iter length(leafnodes)]
    [G,subs,pruned] = graphfuncs.pruneleafbranches(G,subs,pruningThr);
    [G,subs,totdel] = graphfuncs.deleteLoops(G,subs);
end
sprintf('PRUNE DONE in %d sec in %d iter',round(toc(tic_prune)),iter)
Gfilt = G;
subsfilt = subs;

end

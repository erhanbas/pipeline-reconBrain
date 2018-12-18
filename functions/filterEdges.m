function [A,subs] = filterEdges(A,subs,params)
% TODO: deletes edges based on some criteria: user mask or size or some
% heuristics
[A,subs] = maskEdges(A,subs);
end

function [A_,subs_] = maskEdges(A,subs,junk)
[S,Comps] = graphconncomp(A,'DIRECTED',false);
Y = histcounts(Comps,1:S+1);
% filter out results based on manual mask - or allen registration or some other criteria
[~,mY] = max(Y);
[ia,ib]=sort(Y,'descend');
%%
iter=1;
ic_ids=zeros(1,length(Comps));
for ic = junk(:)'
    ic_ids= ic_ids | Comps==ib(ic);
    iter=iter+1;
end
subs_ = subs;
subs_(ic_ids,:) = [];
A_ = A;
A_(ic_ids,:) = [];
A_(:,ic_ids) = [];

end

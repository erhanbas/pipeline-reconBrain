function [branches_neuron,nodeBrid_neuron,branch_pair_distance_neuron,branch_pair_connectivity_neuron,branch_pair_dissimilarity_neuron] = feats4recon(Gneuron,subsneuron,querdist)
%%
if nargin<3
    querdist = 50;
end

[branches_neuron,nodeBrid_neuron] = graphfuncs.graph2branch(Gneuron,subsneuron);
branch_pair_distance_neuron = graphfuncs.branchConn(branches_neuron,subsneuron,nodeBrid_neuron,querdist);
[aa,bb,cc] = find(branch_pair_distance_neuron);
branch_pair_connectivity_neuron = sparse(aa,bb,cc<=querdist,size(branch_pair_distance_neuron,1),size(branch_pair_distance_neuron,1)); %pwcost of branches
branch_pair_dissimilarity_neuron = graphfuncs.calcDists(branches_neuron,branch_pair_connectivity_neuron);

end
function [outputArgs] = affinityBuilder(opt,Ain,subsin)
%AFFINITYBUILDER Summary of this function goes here
%
% [OUTPUTARGS] = AFFINITYBUILDER(INPUTARGS) Explain usage here
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

[S,Comps] = graphconncomp(Ain,'DIRECTED',false);
Yin = histcounts(Comps,1:S+1);
Ain = tril(Ain,-1);

% $Author: base $	$Date: 2017/03/28 13:20:37 $	$Revision: 0.1 $
% Copyright: HHMI 2017
opt
if 0
    tempfold = './temp'
    if ~exist('opt','var')
        load opt
    end
else
    if 0
        % save to local
        tempfold = '/data3/temp/';
    else
        % save to dm11
        [aa,bb,cc] = fileparts(opt.inputh5);
        tempfold = fullfile(pwd,'temp',opt.tag,bb);
        mkdir(tempfold)
    end
    [ia,ib]=sort(Yin,'descend');
end
% ww=[100 100 100]
% h5im = h5read(opt.inputh5,opt.h5prob,[4286 6682 4272]-ww,2*ww);
% figure, imshow3D(permute(h5im,[2 1 3]))

params = opt.params;
if ~isfield(opt,'tag')
    opt.tag = ''
end
tag = opt.tag;

%% prune graph
tr=tic
if opt.prune == 1
    sizeThr = 50;
    
    [Gpruned,subspruned] = pruneGraph(Ain,subsin,sizeThr,opt);
    save(fullfile(tempfold,'Gpruned'),'Gpruned','subspruned')
elseif opt.prune == 2
    load(fullfile(tempfold,'Gpruned'),'Gpruned','subspruned')
else
end
%%
if opt.filterGraph == 1
    [Gfilt,subsfilt] = filterGraph(Gpruned,subspruned,{'filtGsize'},{100});
    save(fullfile(tempfold,'Gfilt'),'Gfilt','subsfilt')

elseif opt.filterGraph == 2
    load(fullfile(tempfold,'Gfilt'),'Gfilt','subsfilt')
else
end
%% %%%%%%%%%%%%%%%%%%%%%
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
% workflow-1: dump all branches after filtering junctions, downsampling
% and smoothing branches
if 1
    workflow1(G,subs,opt)
else
    %% workflow-2: affinity propagation
    workflow2(opt,tempfold,[opt.tag,datestr(now,'yymmdd')],G,subs,branches,nodeBrid)
end
%%
return
%%
% %%
% branchLength = arrayfun(@(x)(length(x.set)), branches);
% branchClustID = arrayfun(@(x)((x.idxC)), branches);
% % cluster with most branches
% Yclust = histcounts(branchClustID,1:max(branchClustID)+1);
% [~,maxclust] = max(Yclust);
% %%
% numbr = length(branches);
% querrylocs = zeros(numbr*2,3);
% for ibr = 1:numbr
%     querrylocs((ibr-1)*2+1:ibr*2,:) = branches(ibr).subs([1 end],:);
% end
%
% %%
% % build kdtree
% [Idx,dist] = rangesearch(querrylocs,querrylocs,2);
% len = cellfun(@length,Idx);
%%
if 0
    querdist = 50;
    connBr = branchConn(branches,subs,nodeBrid,querdist); % pwdist of branches
    [aa,bb,cc] = find(connBr);
    connbwBr = sparse(aa,bb,cc<=querdist,size(connBr,1),size(connBr,1)); %pwcost of branches
    distBr = calcDists(branches,connbwBr);
    save(fullfile(tempfold,['branchConn',tag]),'connBr','distBr')
else
    load(fullfile(tempfold,['branchConn',tag]),'connBr','distBr')
end
if 0
    %%
    swcfold = '/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2015-06-19/Consensus'
    swcfold = './Consensus'
    swcfile{1} = fullfile(swcfold,'2015-06-19_G-001_consensus.swc');
    swcfile{2} = fullfile(swcfold,'2015-06-19_G-003_AZ.swc');
    swcfile{3} = fullfile(swcfold,'2015-06-19_G-004_consensus.swc');
    swcfile{4} = fullfile(swcfold,'2015-06-19_G-005_CA.swc');
    swcfile{5} = fullfile(swcfold,'2015-06-19_G-006_consensus_CA-PB.swc');
    swcfile{6} = fullfile(swcfold,'2015-06-19_G-007-consensus.swc');
    swcfile{7} = fullfile(swcfold,'2015-06-19_G-008_consensus_new.swc');
    swcfile{8} = fullfile(swcfold,'2015-06-19_G-010_agreement.swc');
    swcfile{9} = fullfile(swcfold,'2015-06-19_G-011_consensus.swc');
    [somalocs,swcpixlocs,connMatrix] = loadGT(swcfile,params);
    GT.swcpixlocs = swcpixlocs;
    GT.connMatrix = connMatrix;
    
    %%
    ran = max(subs);
    branchcent=reshape([branches.cent]',3,[])';
    
    tipsubs=[branches.tipsubs];tipsubs = permute(reshape(tipsubs,2,3,[]),[2 1 3]);tipsubs = reshape(tipsubs,3,[]);
    
    [sources,dist] = knnsearch(tipsubs',somalocs,'k',1);
    sources = round(sources/2);
    nsource =length(sources) ;
    
    viz = 1;
    if viz
        % draw GT
        %     close all
        figure,
        %     myplot3(subs(1:10:end,:),{'k+','MarkerSize',2})
        hold on
        valinds=zeros(1,nsource);
        for ii=1:nsource
            if isempty(swcpixlocs{ii})
                continue
            end
            valinds(ii)=1;
            % myplot3(swcpixlocs{ii},{'o','MarkerSize',2})
            gplot3(connMatrix{ii},swcpixlocs{ii},'-')
        end
        %     ran = max(subsin);
        xlim([0 ran(1)])
        ylim([0 ran(2)])
        zlim([0 ran(3)])
        legend(num2str([0 find(valinds)]'))
        leg= {'SK','1','3','4','5','6','7','8','10','11'}
        legend(leg,'location','NorthEast')
        view([0 90])
        set(gca,'Ydir','reverse')
        set(gcf,'UserData',params)
        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'UpdateFcn',@mycallback)
    end
    
    %%
    % class labels of adjacent tiles
    % D(1): Euclidean
    % D(2): theta
    % D(3): PCA
    % D(4): KL
    % D(5 [end]): hit
    
    [w1,w2,w3]= find(distBr{1}); % based on euclidean tip distance
    [w1,w2,w3]= find(connBr); % based on tip2shorthest distance
    mm = w3<50;
    mask = sparse(w1(mm),w2(mm),ones(sum(mm),1),size(distBr{1},1),size(distBr{1},1));
    %
    distin = connBr/1e3 .* distBr{3} .* mask;
    nA = size(distin,1); % number of branches
    % make is symetric
    distin = max(distin,distin');
    
    Yinit = zeros(nA,nsource);
    for ii=1:nsource
        Yinit(sources(ii),ii)=1;
    end
    %%%%%%%%%%%%%%%%%%%%%% QUERIES @@@@@@@@@@@@@@@@
    % distin([38979],38826 ) = eps;
    % distin(38826,38979) = eps;
    % %%%%%%%%%%%%%%%%%%%%%% QUERIES - GLOBAL @@@@@@@@@@@@@@@@
    % Yinit([132375, 71798,117809,118728,120008,155836,11342],2) = 1;
    % Yinit([104562,153337],3) = 1;
    % Yinit([132922,126964,131488,132585],1) = 1;
    %%%%%%%%%%%%%%%%%%%%%% QUERIES - LOCAL @@@@@@@@@@@@@@@@
    % Yinit([132375, 71798,117809,118728,120008,155836,11342],2) = 1;
    % Yinit([104562,153337],3) = 1;
    % Yinit([132922,126964,131488,132585],1) = 1;
    
    % one pass
    alpha = 1000;
    conf = .5;
    labY = Yinit;
    sig = .01;
    
    %%%%%%%%%%%%%%%%%%%%%% ------ @@@@@@@@@@@@@@@@@
    L = getlaplacian(distin,sig);
    [P,probs,clust] = clusterGraph(L,alpha,labY,nsource);
    clust(isnan(P)) = 0;
    userdatin = {params,subs,nodeBrid,branchcent,opt}
    fignum=8
    vizClust(fignum,userdatin,clust,branches,GT,ran)
    %%
    [w1,w2,w3]= find(distBr{1}); % based on euclidean tip distance
    [w1,w2,w3]= find(connBr); % based on tip2shorthest distance
    mm = w3<50;
    mask = sparse(w1(mm),w2(mm),ones(sum(mm),1),size(distBr{1},1),size(distBr{1},1));
    %
    distin = distBr{1}/1e3 .* distBr{3} .* mask;
    fignum=90
    vizDemo(fignum,userdatin,clust,branches,GT,ran)
    
    %%
    fignum=8
    vizClust(fignum,userdatin,clust,branches,GT,ran)
    %%
    traverseNeuron(distin,userdatin,sources(2),connBr, distBr,sig,labY,branches,GT,ran)
    
    %%
    fignum=10
    vizComp(fignum,2,userdatin,clust,branches,GT,ran)
    
    %%
    Y = histcounts(clust(isfinite(P)),1:nsource+1);
    [~,mY] = max(Y);
    %% draw
    figure(7),
    % clf,
    % cla
    userdatin = {params,subs,nodeBrid}
    set(gcf,'UserData',userdatin)
    vizClust(clust,branches,swcpixlocs,ran)
    % vizGraph(clust,branches,swcpixlocs,ran)
    %%
    que1 = [15524 7018 2818];
    que2 = [15361 6995 2824];
    [aa,bb]=min(pdist2(subs,que2));
    ii=nodeBrid{bb};round([aa,bb,ii])
    hold on, myplot3(branches(130993).subs,'go')
    
    %%
    figure(100), clf
    hold on
    myplot3(subsin,'.')
    myplot3(subspruned,'d')
    myplot3(subsfilt,'o')
    rr=[15957 6985 2814]
    
    kk=500
    xlim([rr(1)+[-2 1]*kk])
    ylim([rr(2)+[-1 1]*kk])
    zlim([rr(3)+[-1 1]*kk])
    view([0 90])
    set(gca,'Ydir','reverse')
    set(gcf,'UserData',params)
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@mycallback)
    %%
    que1 = [15524 7018 2818];
    que2 = [15361 6995 2824];
    [aa,bb]=min(pdist2(subsfilt,que2));
    ii=nodeBrid{bb};round([aa,bb,ii])
    hold on, myplot3(branches(ii).subs,'go')
    
    
    %%
    kk=[100 100 100];
    sr=[11866 6067 1289];
    h5im = permute(h5read(opt.inputh5,opt.h5prob,sr-kk,2*kk),[2 1 3]);
    figure, imshow3D(h5im)
    % myplot3(subs(roots,:),{'ms','MarkerFaceColor','g','MarkerSize',20})
    %%
    figure, myplot3
    
    
    
    %%
    Comps = conncomp(graph(max(Ain,Ain')));
    S = max(Comps);
    Y = histcounts(Comps,1:S+1);
    %%
    % Gfilt
    GG = G;
    subb = subs;
    
    AA = GG.adjacency;
    Comps = conncomp(GG);
    S = max(Comps);
    Y = histcounts(Comps,1:S+1);
    %
    [~,mY] = max(Y);
    % AA = max(Ain,Ain');
    A_S=AA(Comps==mY,Comps==mY);
    su = subb(Comps==mY,:);
    
    figure, gplot3(A_S,su)
    
    
    
    %%
    if 1
        algorithmkeys = {'spb','dij','bel','spa'};
        algorithm = 2;
        debug_level = 0;
        directed = 0;
        W = [];
        colscomp = jet(length(validC));
        colscomp = colscomp(randperm(size(colscomp,1)),:);
        %
        runtic = tic;
        Eout = [];
        iter = 0;
        N = max(validC);
        %%
        skipthese = zeros(1,N);
        parfor mC=validC%1:size(Y,2)%ib(2:500)%
            if Y(mC)>opt.sizethreshold %& ismember(mC,validC) %& ~any(mC==skipThese)
                %%
                %             swcoutfolder = fullfile(opt.outfolder,'full');
                %             if 1
                %                 [aa,bb,cc] = fileparts(opt.inputh5);
                %                 pre = ['-',bb(end-6:end-4)];
                %             else
                %                 pre = '';
                %             end
                %             swcoutname = sprintf('auto%s',pre);
                %             fragname = sprintf('%s_cc-%04d.swc',swcoutname,mC);
                %             outname = fullfile(swcoutfolder,fragname);
                %             if exist(outname,'file')
                %                 skipthese(mC) = 1;
                %             end
            else
                skipthese(mC) = 1;
            end
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        skipthese(ib(1))=1 %% $$ FOR FIX MEMORY ERROR< NEED TO RECODE THIS PART
        %%
        for mC=ib(1)%1:validC(end)
            % for each cluster run reconstruction
            if skipthese(mC)
                continue
            end
            %%
            subidx = find(Comps==mC);
            subs_ = subs(subidx,:);
            nidx = length(subidx);
            Asub = A_(subidx,subidx); % get lower portion to make it directed
            [Asub,subs_] = A2G(Asub,subs_);
            %%
            graph2branch(Asub,subs_)
            %%
            % figure,
            % gplot3(A,subs)
            % hold on
            if 1
                pruned =1;
                iter = 0;
                while pruned
                    iter = iter+1;
                    [Asub,subs_,pruned,branchlist] = simplifiedGraph(Asub,subs_);
                    %     gplot3(A,subs)
                end
                if isempty(branchlist) | ~all(isfinite([branchlist.length]))
                    continue
                end
                pD = branch2dist(branchlist,subs_);
                nA=size(pD,1);
                pD(1:nA+1:nA^2)=inf;
                connG = pD<20;
                gC = graph(connG);
                dists = calcWeights(branchlist,connG,Asub,subs_(:,:));
                distin = dists(:,:,2);
                % make is symetric
                weights = exp(-max(distin,distin'));
                filename = sprintf('./affinitygraphs/clust_%05d',mC);
                locs = pix2um(params,subs_); % center anisotropy to compansate imresize
                dumper(filename,branchlist,locs,weights);
            end
        end
        toc(runtic)
        
    end
    %     hits = cellfun(@length,nodeBrid);
    %     depth = max(hits);
    %     ed = cell(1,size(nodeBrid,2));
    %     spnodeBrid = zeros(size(nodeBrid,2),depth);
    %     for ii=1:size(nodeBrid,2)
    %         tar = nodeBrid{ii};
    %         spnodeBrid(ii,1:length(tar))=tar;
    %     end
    
end
end
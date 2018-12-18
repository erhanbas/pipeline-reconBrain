function workflow1(G,subs,opt)
params = opt.params;
if ~isfield(opt,'medianFiltSize')
    opt.medianFiltSize = 10;
end
CompsC = conncomp(G,'OutputForm','cell');
A = G.adjacency;
A_ = tril(A,-1);
S = length(CompsC);
Y = cellfun(@length,CompsC);
validC = 1:size(Y,2);
[ia,ib]=sort(Y,'descend');
algorithmkeys = {'spb','dij','bel','spa'};
algorithm = 2;
debug_level = 0;
directed = 0;
W = [];
colscomp = jet(length(validC));
colscomp = colscomp(randperm(size(colscomp,1)),:);
SX = params.sx;
SY = params.sy;
SZ = params.sz;
voxres = [SX SY SZ]/2^(params.level)/1e3; % in um
params.voxres = voxres;
runtic = tic;
Eout = [];
iter = 0;
N = max(validC);
%%
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
    parpool(feature('numcores'))
else
    poolsize = poolobj.NumWorkers
end
skipthese = zeros(1,N);
skipthese(Y<=opt.sizethreshold) = 1;
%%
initialpass = 1;
acc=0;
if ~initialpass
    parfor mC=validC%1:size(Y,2)%ib(2:500)%
        if Y(mC)>opt.sizethreshold %& ismember(mC,validC) %& ~any(mC==skipThese)
            %%
            [aa,bb,cc] = fileparts(opt.inputh5);
            bb_ = strsplit(bb,'_');
            pre = bb_{end-1};%['-',bb(end-6:end-4)];
            swcoutname = sprintf('auto%s',pre);
            fragname = sprintf('%s_cc-%04d.swc',swcoutname,mC);
            swcoutfolder = fullfile(opt.outfolder,'full');
            outname = fullfile(swcoutfolder,fragname);
            if exist(outname,'file')
                skipthese(mC) = 1;
                acc=acc+1;
            end
        else
            skipthese(mC) = 1;
        end
    end
end
disp('done skip')
%% 
try parfor_progress(0);catch;end
parfor_progress(max(validC))
parfor mC=validC
    parfor_progress
    % for each cluster run reconstruction
    if skipthese(mC)
        continue
    end
    %%
    iter = mC;%iter+1;
    subidx = CompsC{mC};%find(Comps==mC);
    subs_ = subs(subidx,:); % get back to matlab image coordinates
    nidx = length(subidx);
    % get lower portion to make it directed
    Asub = A_(subidx,subidx); % faster
    %Gsub = G.subgraph(subidx);
    
    leafs = find(sum(Asub,2)==0);%find(sum(max(Asub,Asub'))==1,1);
    
    [eout] = graphfuncs.buildgraph(Asub,leafs(1));
    inupdate.dA = sparse(eout(:,1),eout(:,2),1,nidx,nidx);
    inupdate.D = ones(nidx,1);
    inupdate.R = ones(nidx,1);
    inupdate.X = subs_(:,1);
    inupdate.Y = subs_(:,2);
    inupdate.Z = subs_(:,3);
    %%
    deleteThese = NaN;
    while length(deleteThese)
        [inupdate, deleteThese] = prunTree(inupdate,opt.lengthThr,voxres);
        if opt.viz
            hold on
            gplot3(inupdate.dA,[inupdate.X,inupdate.Y,inupdate.Z]);
            drawnow
        end
    end
    %%
    % shuffle root to one of the leafs for efficiency and not
    % splitting long stretches into half
    if size(inupdate.dA,1)>1
        [eoutprun] = graphfuncs.buildgraph(inupdate.dA);
        nidx = max(eoutprun(:));
        inupdate.dA = sparse(eoutprun(:,1),eoutprun(:,2),1,nidx,nidx);
    else
    end
    if opt.viz
        hold on
        gplot3(inupdate.dA,[inupdate.X,inupdate.Y,inupdate.Z],'LineWidth',3);
        drawnow
    end
    %%
    if length(inupdate.dA)<opt.sizethreshold
        continue
    end
    %%
    [inupdate] = smoothtree(inupdate,opt);
    %%
    % [L,list] = getBranches(inupdate.dA);
    if 0
        outtree = inupdate;
    else
        % outtree_old = downSampleTree(inupdate,opt);
        outtree = sampleTree(inupdate,opt);
    end
    if opt.viz
        cla
        gplot3(inupdate.dA,[inupdate.X,inupdate.Y,inupdate.Z],'LineWidth',3);
        hold on
        gplot3(outtree.dA,[outtree.X,outtree.Y,outtree.Z],'--','LineWidth',3);
        drawnow
    end
    %%
    Aout = outtree.dA;
    nout = size(Aout,1);
    XYZout = [outtree.X,outtree.Y,outtree.Z]-1;
    Rout = outtree.R;
    Dout = outtree.D;
    % transform location
    XYZout = pix2um(params,XYZout); % center anisotropy to compansate imresize
    %
    
    At = Aout+Aout';
    [DISC,PRED,CLOSE] = graphtraverse(At,1,'Method','DFS');
    At(1:end,:) = At(DISC,:);
    At(:,1:end) = At(:,DISC);
    XYZout = XYZout(DISC,:);
    Rout = Rout(DISC,:);
    Dout = Dout(DISC,:);
    %%
    if strcmp(version('-release'),'2015b')
        [dist,pred] = graphalgs(algorithmkeys{algorithm},debug_level,directed,At,1);
    else
        [dist,path,pred] = graphshortestpath(At,1,'DIRECTED',false);
    end
    swcData = [[1:nout]' Dout XYZout Rout pred(:)];
    swcData(1,7) = -1;
    offset = min(swcData(:,3:5),[],1);
    swcData(:,3:5) = swcData(:,3:5)-ones(size(swcData,1),1)*offset;
    % check quadrant
    quadid = checkQuant(opt.params.outsiz,median([outtree.X,outtree.Y,outtree.Z]));
    quadid = 0;
    %%
    %% WRITE components
    if opt.writefull
        %%
        if quadid
            swcoutfolder = fullfile(opt.outfolder,'full',num2str(quadid));
        else
            swcoutfolder = fullfile(opt.outfolder,'full');
        end
        
        if 1
            [aa,bb,cc] = fileparts(opt.inputh5);
            bb_ = strsplit(bb,'_');
            pre = bb_{end-1};%['-',bb(end-6:end-4)];
            %                 [aa,bb,cc] = fileparts(opt.inputh5);
            %                 pre = ['-',bb(end-6:end-4)];
        else
            pre = '';
        end
        swcoutname = sprintf('auto%s',pre);
        fragname = sprintf('%s_cc-%04d.swc',swcoutname,mC);
        outname = fullfile(swcoutfolder,fragname);
        if ~exist(fileparts(outname),'dir')
            mkdir(fileparts(outname))
        end
        % write swc
        fid = fopen(outname,'w');
        mytxt = sprintf('# Generated by pw skel algorithm\n');
        fprintf(fid,'%s',mytxt);
        mytxt = sprintf('# OFFSET %.6f %.6f %.6f\n',offset);
        fprintf(fid,'%s',mytxt);
        mytxt = sprintf('# COLOR %f,%f,%f\n',colscomp(iter,1),colscomp(iter,2),colscomp(iter,3));
        fprintf(fid,'%s',mytxt);
        mytxt = sprintf('# NAME %s\n',fragname);
        fprintf(fid,'%s',mytxt);
        fprintf(fid,'%d %d %f %f %f %d %d\n',swcData');
        fclose(fid);
    end
    %% WRITE fragments
    if opt.writefrag
        %%
        %Aout = outtree.dA;
        %nout = size(Aout,1);
        XYZ = [outtree.X,outtree.Y,outtree.Z]-1; % -1 one is heuristic/bug
        R = outtree.R;
        D = outtree.D;
        [L,list] = getBranches(outtree.dA);
        % write to disk as swc file
        cols = jet(length(L));
        cols = cols(randperm(length(L)),:);
        for ii=1:length(L)
            %%
            out = L(ii).set;
            if isempty(out)
                continue
            else
                %out(end) = [];
                if isempty(out)
                    continue
                end
            end
            %%
            XYZout = XYZ(out,:);
            Rout = R(out);
            Dout = D(out);
            
            % transform location
            XYZout = pix2um(opt.params,XYZout); % center anisotropy to compansate imresize
            %XYZout = pix2um(opt,[XYZout(:,1)*opt.aniscale(1) XYZout(:,2)*opt.aniscale(2) (XYZout(:,3)-1)*opt.aniscale(3)+1]); % center anisotropy to compansate imresize
            nout = size(XYZout,1);
            swcData = [[1:nout]' Dout XYZout Rout [0:nout-1]'];
            swcData(1,7) = -1;
            offset = mean(swcData(:,[3:5]),1);
            %%
            % write swc file
            if 1
                if 1
                    [aa,bb,cc] = fileparts(opt.inputh5);
                    bb_ = strsplit(bb,'_');
                    pre = bb_{end-1};%['-',bb(end-6:end-4)];
                else
                    pre = '';
                end
                
                if quadid
                    swcoutfolder = fullfile(opt.outfolder,'frags',num2str(quadid));
                else
                    swcoutfolder = fullfile(opt.outfolder,'frags');
                end
                
                swcoutname = sprintf('auto%s',pre);
                %swcoutname = 'myskel';
                fragname = sprintf('%s_cc-%04d_branch-%04d.swc',swcoutname,mC,ii);
                outname = fullfile(swcoutfolder,fragname);
                if ~exist(fileparts(outname),'dir')
                    mkdir(fileparts(outname))
                end
                
                swcData(:,3:5) = swcData(:,3:5)-ones(size(swcData,1),1)*offset;
                % write swc
                fid = fopen(outname,'w');
                mytxt = sprintf('# Generated by pw skel algorithm\n');
                fprintf(fid,'%s',mytxt);
                mytxt = sprintf('# OFFSET %.6f %.6f %.6f\n',offset);
                fprintf(fid,'%s',mytxt);
                mytxt = sprintf('# COLOR %f,%f,%f\n',cols(ii,1),cols(ii,2),cols(ii,3));
                fprintf(fid,'%s',mytxt);
                mytxt = sprintf('# NAME %s\n',fragname);
                fprintf(fid,'%s',mytxt);
                fprintf(fid,'%d %d %f %f %f %d %d\n',swcData');
                fclose(fid);
            end
        end
    end
end
parfor_progress(0)
toc(runtic)
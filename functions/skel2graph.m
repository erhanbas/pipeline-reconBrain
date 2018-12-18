function [subs,edges,A,weights_] = skel2graph(opt)

mydir = dir([opt.skelfolder,'*.txt']);
% check the format for a non zero file
mysiz = [mydir.bytes];
fid = fopen(fullfile(opt.skelfolder,mydir(find(mysiz>0,1)).name));
tline = fgetl(fid);
fclose(fid);
tlines = strsplit(tline,' ');
numcloumns = size(tlines,2); % first 2 are indicies, rest are weights
format = ['%f %f', sprintf('%s',repmat(' %f',ones(1,numcloumns-2)))];

fclose all;
clear pairs edges
pairs = cell(1,length(mydir));
%%
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
    parpool(feature('numcores'))
else
    poolsize = poolobj.NumWorkers
end
%%
parfor idx = 1:length(mydir)
    if mysiz(idx)==0 %| mysiz>20e6
        continue
    end
    % read text file line by line
    fid = fopen(fullfile(opt.skelfolder,mydir(idx).name));
    tmp = textscan(fid,format);
    fclose(fid);
    pairs{idx} = cat(2,tmp{:})';
end
edges = [pairs{:}]';clear pairs;

%%
clc
clear subs
[keepthese,ia,ic] = unique(edges(:,[1 2]));
[subs(:,1),subs(:,2),subs(:,3)] = ind2sub(opt.params.outsiz([1 2 3]),keepthese);
edges_ = reshape(ic,[],2);
weights_ = edges(ia,3:end);
if isempty(edges_);return;end
if isempty(weights_);weights_ = ones(size(edges_,1),1);end

selfinds = find((edges_(:,1)==edges_(:,2)));
if ~isempty(selfinds);edges_(selfinds,:)=[];end
A = sparse(edges_(:,1),edges_(:,2),1,max(edges_(:)),max(edges_(:)));
A = max(A',A);

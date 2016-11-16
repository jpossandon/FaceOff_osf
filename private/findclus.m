function  [cluster] = findclus(data,neighboursmat,tipo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  = findclus(data,neighboursmat,type)
% find cluster in time and channel dimension
% data          = logical array time x channels , 1s are the event to cluster
% neighboursmat = channels neighbour matrix
% tipo          = 'id'  , gives to every cluster element the same corelative number
%                 'sum' , gives to every cluster element the number of the
%                 sum of element that are part of the cluster (used for TCFE)
%
% JPO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data                = double(data);
[r,c]               = size(data);
clusaux             = zeros(r+2,c);    % two more rows to find cluster at the begining and end
clusaux(2:end-1,:)  = data;
clusaux             = diff(clusaux);   % find start and end of continuous segments of 1s
[~,j]               = ind2sub(size(clusaux),find(clusaux==1));
[~,jj]              = ind2sub(size(clusaux),find(clusaux==-1));
indx_chclus         = [find(clusaux==1)-j+1,find(clusaux==-1)-jj+1-1];
for e = 1:size(indx_chclus,1)          % here we give to every cluster (only in the time[rows] dimension)
    data(indx_chclus(e,1):indx_chclus(e,2)) = e;
end
data = data';
% this is adapted from fieldtrip findcluster, it extends the cluster in the channel dimension 
replaceby=1:max(data(:));
for chan = 1:c
    neighbours    = find(neighboursmat(chan,:));
    for nbindx=neighbours
        indx = find((data(chan,:)~=0) & (data(nbindx,:)~=0));
        for i=1:length(indx)
          a = data(chan, indx(i));
          b = data(nbindx, indx(i));
          if replaceby(a)==replaceby(b)
            % do nothing
            continue;
          elseif replaceby(a)<replaceby(b)
            % replace all entries with content replaceby(b) by replaceby(a).
            replaceby(find(replaceby==replaceby(b))) = replaceby(a); 
          elseif replaceby(b)<replaceby(a)
            % replace all entries with content replaceby(a) by replaceby(b).
            replaceby(find(replaceby==replaceby(a))) = replaceby(b); 
          end
        end
    end
end

% renumber the cluster to the number of element withint th cluster(the support)
switch tipo
    case ('sum')
        cluster = zeros(size(data));
        for uniquelabel=unique(replaceby(:))'
            clusauxindx = ismember(data(:),find(replaceby==uniquelabel)); 
            cluster(clusauxindx) = sum(clusauxindx);
        end
    case ('id') 
        cluster = zeros(size(data));
        a = 1;
        for uniquelabel=unique(replaceby(:))'
            clusauxindx = ismember(data(:),find(replaceby==uniquelabel)); 
            cluster(clusauxindx) = a;
            a = a+1;
        end
end
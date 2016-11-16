function  [clusters] = clustereeg(st,H,elec,ch,times)
     HPos               = zeros(size(H));
     HPos(st>0 & H==1)   = 1;
     HNeg               = zeros(size(H));
     HNeg(st<0 & H==1)   = 1;

     if ch>1 && times>1
        [clusterp] = findclus(squeeze(HPos),elec.channeighbstructmat,'id');
        [clustern] = findclus(squeeze(HNeg),elec.channeighbstructmat,'id');
    elseif ch==1 && times>1
        [clusterp] = findclus(squeeze(HPos),elec.channeighbstructmat,'id');
        [clustern] = findclus(squeeze(HNeg),elec.channeighbstructmat,'id');
    elseif ch>1 && times==1
        [clusterp] = findclus(squeeze(HPos),elec.channeighbstructmat,'id');
        [clustern] = findclus(squeeze(HNeg),elec.channeighbstructmat,'id');
     end
     if any(clusterp(:))
        for cn = 1:max(clusterp(:))
%             auxclusp(cn) = sum(squeeze(st(:,:,find(clusterp==cn))));
        	auxclusp(cn) = sum(squeeze(st(find(clusterp'==cn))));
        end
     else
         auxclusp = [];
     end
     if any(clustern(:))
        for cn = 1:max(clustern(:))
%              auxclusn(cn) = sum(squeeze(st(:,:,find(clustern==cn))));
            auxclusn(cn) = sum(squeeze(st(find(clustern'==cn))));
        end
     else
         auxclusn = [];
     end

         if ~isempty(auxclusp)
            clusters.maxt_pos         = auxclusp;
            clusters.clus_pos         = clusterp;
        else
            clusters.maxt_pos         = [];
            clusters.clus_pos         = [];
        end
        if ~isempty(auxclusn)
            clusters.maxt_neg         = auxclusn;
            clusters.clus_neg         = clustern;
        else
            clusters.maxt_neg         = [];
            clusters.clus_neg         = [];
        end
%                     result.Bt(:,:,b)        = T;
        if isempty(auxclusp) && isempty(auxclusn)
            clusters.MAXst = 0;
            clusters.MAXst_noabs = [0 0];
        else

            clusters.MAXst = max(abs([auxclusp,auxclusn]));
            if isempty(auxclusp) 
                clusters.MAXst_noabs = [0 max([auxclusn])];
            elseif isempty(auxclusn) 
                clusters.MAXst_noabs = [max([auxclusp]) 0];
            else
               clusters.MAXst_noabs = [max([auxclusp]) min([auxclusn])]; 
            end
        end
 
    
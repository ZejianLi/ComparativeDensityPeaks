function ClusterIdx = assignClusterIdx(cluster_centers, ordrho, DeltaParent)

% this function is to propagate the labels from the cluster centers to the descendants 
% in the Delta tree

% the number of data
N = length(ordrho);
% the number of cluster
K = sum(DeltaParent(cluster_centers)==cluster_centers);

% cluster index
ClusterIdx = zeros(1,N);
ClusterIdx(cluster_centers) = 1:K; % ¥ÿÀ˘ Ù

% ordrho indicates the descending order of rho values
for i = 1:N
    thisIdx = ordrho(i);
    path = false(1,N);
    path(thisIdx) = true;
    while ClusterIdx(DeltaParent(thisIdx))==0
        if( thisIdx==DeltaParent(thisIdx))
            error('rho max may be not in cluster center');
        end
        if(path(thisIdx))
            error('cycle delta parent');
        end
        thisIdx = DeltaParent(thisIdx);
        path(thisIdx) = true;
    end
    ClusterIdx(path) = ClusterIdx(DeltaParent(thisIdx)); % follows the parent
end

end
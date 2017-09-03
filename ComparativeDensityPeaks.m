%% The Comparative Density Peaks algorithm
% input
% features  : the d*N input, where d is the dimension of the sample and N is the number of sample
% K         : the predefined number of clusters
% dc        : the cut-off parameter d_c
% DataSetName:the name of the dataset
% varargin{1}:if it is one, the algorithm is in oracle mode, which construct the \Delta-tree
%               according to the labels specified by labels
% output
% ClusterIdx: the cluster assignment on data
% rho       : the density
% delta     : the delta distance
% cluster_centers   : the indexes of cluster centers
% DeltaParent       : the indexes of \Delta
% theta     : the theta distance


function [ClusterIdx,rho,delta,cluster_centers,DeltaParent,gamma_,theta] = ComparativeDensityPeaks(features, K, dc, DataSetName, varargin)

    disp('ComparativeDensityPeaks called');

    rhoMeasure = @mutualKnnGaussianRho;
    distMeasure = @l2distMeasure;

    % whether in Oracle mode
    oracle=false;
    theclass=[];
    if(nargin>5)
        if(varargin{1}==1)
            disp('Oracle mode');
            oracle=true;
            theclass=varargin{2};
        end
    end


    N=size(features,2); % num of data points

    disp(['calculate rho with ',func2str(rhoMeasure)]);
    [rho,AdjaccentKnnMatrix] = rhoMeasure(features, dc, DataSetName);
    rho = rho + 1e-4 .* ( 1:numel(rho) );


    disp('calculate delta with mutualKnnMatrix for geodesic distance');
    [delta,DeltaParent,tau] = deal(zeros(1,N));

    knnDist = sparse([],[],[],N,N,nnz(AdjaccentKnnMatrix));
    for i = 1:N
        knnDist(AdjaccentKnnMatrix(:,i),i) = distMeasure(features(:,i),features(:,AdjaccentKnnMatrix(:,i)))+eps;
    end


    [~, maxRhoIdx] = max(rho);% the index of the maximal Rho
    [~, minRhoIdx] = min(rho);% the index of the smallest Rho
    DeltaParent(maxRhoIdx) = maxRhoIdx; % the density peak has the delta parent as itself


    for i = 1:N % calculate delta and Delta for every point
        if(i ~= maxRhoIdx)
            thisKnnIdx = AdjaccentKnnMatrix(:, i); % N*1
            kNNRho = rho(thisKnnIdx); % 1*k
            kNNDist = nonzeros(knnDist(:,i)); % k*1
            kNNRhoLargerIdx = kNNRho>rho(i); % 1*k
            
            if(oracle)
                kNNRhoLargerIdx = kNNRhoLargerIdx & theclass(thisKnnIdx)==theclass(i);
            end
                
            if(any(kNNRhoLargerIdx)) % one of the k nearest neighbors has a larger Rho
                [delta(i),idx] = min(kNNDist(kNNRhoLargerIdx));
                tmp = find(thisKnnIdx); % 1*k
                tmp = tmp(kNNRhoLargerIdx); % knn rho larger
                DeltaParent(i) = tmp(idx);            
            else
                terminal = @(j)rho(j)>rho(i); % termination condition of the Dijkstra 
                if(oracle)
                    terminal = @(j)rho(j)>rho(i) & theclass(j) == theclass(i);
                end
                
                [delta(i),DeltaParent(i)] = geodesic_delta(i,knnDist,rho,terminal,features,distMeasure);
                
                
            end
            if(isinf(delta(i)))
                error(['delta as inf in node ',num2str(i)]);
            end
            

            kNNRhoSmallerIdx = kNNRho<rho(i);
            
            if(oracle)
                kNNRhoSmallerIdx = kNNRhoSmallerIdx & theclass(thisKnnIdx) == theclass(i);
            end
            
            if(any(kNNRhoSmallerIdx)) % one of the KNN has a smaller Rho
                tau(i) = min(kNNDist(kNNRhoSmallerIdx));
            else % 
                if(i~=minRhoIdx) % not point with the smallest density
                    terminal=@(j)rho(j)<rho(i); % termination condition of Dijkstra
                    tau(i)=geodesic_delta(i,knnDist,rho,terminal,features,distMeasure);
                else % 
                    tau(i)=delta(i);
                end
            end        
        end
    end

    delta(maxRhoIdx) = max(delta); % to adjust the max delta value
    theta = delta - tau; % 

    gamma_ = rho.*theta+1e-8;
    [~,~,idxK] = minK(-gamma_, K); % the max K gamma_ values
    cluster_centers = idxK;

    if(oracle == 1) % in oracle mode, choose the maximal density in the group
        [~,cluster_centers] = arrayfun(@(i)max(rho.*(theclass==i)),1:K);
    end
        

    DeltaParent(cluster_centers) = cluster_centers; % cluster centers have themself as the Delta Parent
    [~,ordrho] = sort(rho,'descend');
    if(~any(ordrho(1) == cluster_centers))
        error('max rho not in clustering centers');
    end
    disp('clustering according to DeltaParent');
    ClusterIdx = assignClusterIdx(cluster_centers,ordrho,DeltaParent);

end





function [GeodesicDistance, DeltaParent] = geodesic_delta(i,knnDist,rho,ter,features,DistMeasure)
    GeodesicDistance = 0;
    DeltaParent = 0;
    N = size(knnDist,2);
    unset = true(N,1); 
    unset(i) = false; 
    dist = inf(N,1); % the final distance
    dist(i) = 0;
    thisIdx = i;
    possibleOnes = ter(1:N); % 
    anyPossible = any(possibleOnes); % 
    while(true)
        thisIdxKnnDist = knnDist(:,thisIdx); % 1*N
        thisIdxKnnDistNonNeg = thisIdxKnnDist > 0;
        dist(thisIdxKnnDistNonNeg) = min(dist(thisIdx) + thisIdxKnnDist(thisIdxKnnDistNonNeg),dist(thisIdxKnnDistNonNeg)); % about k elements

        [v,minone] = min(dist(unset)); % the smallest distance in the unset points

          if isinf(v) || ~anyPossible % v is inf, all the unset points are not reachable
            % compromise, follow the delta parent, but the distance is another
            % the original DeltaParent
            [~,minone] = min(DistMeasure(features(:,i),features(:,possibleOnes))); % the nearest one with a higher density
            tmp = find(possibleOnes);
            if isempty(tmp) % in orachle mode, the density peak
                [~,minone] = min( DistMeasure(features(:,i) , features(:, rho>rho(i) ) ) ); % the nearest one with a higher density
                tmp = find( rho>rho(i) );            
            end
            originDeltaParent = tmp(minone);
            % the distance from the reachable point to origin Delta Parent
            reachable2parent = DistMeasure( features(:,originDeltaParent) , features(:,~isinf(dist)) );
            % the smallest one
            [v,idx] = min(reachable2parent);
            reachableDist = dist(~isinf(dist));        
            % the sum as the compromist
            GeodesicDistance = v+reachableDist(idx);
            DeltaParent = originDeltaParent;
            return
        end

        ind = find(unset,minone,'first');
        setOne = ind(end);
        if( possibleOnes(1,setOne) )
            GeodesicDistance = v;
            DeltaParent = setOne;
            return
        end
        unset(setOne) = false;
        thisIdx = setOne;    
    end

end



function [Rho, MutualKnnMatrix] = mutualKnnGaussianRho(features, dc, varargin)
% This function is to compute Rho in mutual KNN graph

    N = size(features,2); 

    DataSetName = [];
    if nargin>2
        DataSetName = varargin{1};
    end

    P = 0.01;
    [~, knnMatrix, knnIdx, knnDist] = mutualKnnMatrix(features, P, DataSetName);


    MutualKnnMatrix = knnMatrix | knnMatrix';

    sb = sum(features.^2,1); 

    Rho = arrayfun(...
    @(i) sum( exp( ...
    bsxfun( @minus,...
    bsxfun( @minus,  2*full(features(:,i)'*features(:,MutualKnnMatrix(:,i)) ) , sb(i) ),...
    sb(MutualKnnMatrix(:,i) ) ) / (dc.^2) ...
    ) ) ,...
    1:N);

end
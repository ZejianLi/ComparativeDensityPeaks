function [thisKnn, KnnMatrix, KnnIdx, KnnDist] = mutualKnnMatrix(features, varargin)

    % this function is to get the KNN relation

    % p is k/N
    p = 0.01;
    if(nargin>1 && ~isempty(varargin{1}))
        p = varargin{1};
    end

    global P;
    if ~isempty(P)
        p = P;
    end

    % the name of dataset
    DataSetName = [];
    if(nargin>2 && ~isempty(varargin{2}))
        DataSetName = varargin{2};
        disp(DataSetName);
    end


    N = size(features,2); % number of data points
    thisKnn = min((max(ceil(N*p),5)),N); % the k in knn

    KnnMatrixFileName = ['KNNMatrix/KnnMatrix_',DataSetName,'.mat'];

    if ~exist('KNNMatrix','dir')
        mkdir('KNNMatrix');
    end

    KnnIdx = [];
    KnnDist = [];
    KnnMatrix = [];
    
    % to load the calculated result given the dataset name
    if(~isempty(DataSetName) && exist(KnnMatrixFileName,'file')) % 
        
        load(KnnMatrixFileName);
        
        disp(['load ',KnnMatrixFileName]);
    end
        
    if isempty(KnnIdx) || isempty(KnnDist) %
        

        [KnnIdx,KnnDist] = knnsearch(features',features','K',N,'Distance','euclidean');% 最近的是自己，所以不应该考虑
        
        KnnIdx = KnnIdx(:,2:end)'; % Knn*N
        KnnDist = KnnDist(:,2:end)'; % Knn*N

        if(~isempty(DataSetName))
            save(KnnMatrixFileName,'KnnIdx','KnnDist','DataSetName','-v7.3');
            disp([KnnMatrixFileName,' saved!']);
        end
        
    end

    
    KnnIdx = KnnIdx(1:thisKnn,:); % 
    KnnDist = KnnDist(1:thisKnn,:);


    KnnMatrix = spalloc(N,N,2*thisKnn*N); % 
    for i1=1:N
        KnnMatrix(KnnIdx(:,i1),i1) = true;
    end
    KnnMatrix = (KnnMatrix == true);


end

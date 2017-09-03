function [result,tmpK,idxK] = minK(vector,K)
% to search the smallest K elements in the input vector
% output:
% result is the Kth smallest value
% tmpK is the K smallest values
% idxK is the index

N = length(vector);
if N <= K
    [result,tmpK,idxK] = deal(vector(K),vector,1:K);
    return
end

if N <= 1e5
    [sortedVal,sortedIdx] = sort(vector,'ascend');
    tmpK = sortedVal(1:K);
    idxK = sortedIdx(1:K);
else
    
    tmpK = vector(1:K);
    idxK = 1:K;
    for i = (K+1):N
        thisone = vector(i);
        [val,idx] = min(thisone - tmpK);
        if val<0
            tmpK(idx) = thisone;
            idxK(idx) = i;
        end
    end
end
result = tmpK(end);

end

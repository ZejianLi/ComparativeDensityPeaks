function [dist,sB]=l2distMeasure(A,B,varargin)
	if( isempty(A) || isempty(B))
	    dist=0;
	    return;
	end
	% dist=pdist2(A',B');
	sA = (sum(A.^2, 1)); 
	if nargin>2 && ~isempty(varargin{1}) && length(varargin{1})==size(B,2)
	    sB = varargin{1};
	else
	    sB = (sum(B.^2, 1));
	end
	% sB= sum( bsxfun(@power,B,2),1);
	dist = bsxfun(@plus,bsxfun(@minus, sA',2*A'*B), sB);
	dist = max(dist,0);
    n =  size(dist);
    dist=sqrt(dist);
end
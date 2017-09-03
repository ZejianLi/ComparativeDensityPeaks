function rho = gaussianRhoMeasure(features, dc, varargin)
% Gaussian rho 

N=size(features,2); % num of data points
rho=zeros(1,N);

sb=sum(features.^2,1);

rho=arrayfun(...
@(i) sum( exp( bsxfun( @minus, bsxfun(@minus, 2*full(features(:,i)'*features), sb(i)), sb)/(dc.^2) ) ) ,...
1:N);

rho=max(rho-1,eps);



end
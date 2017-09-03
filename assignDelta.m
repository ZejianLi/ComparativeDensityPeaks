function [delta, DeltaParent, ordrho] = assignDelta(features, rho, distMeasure)

N=size(features,2);

[delta,DeltaParent]=deal(zeros(1,N));

% sort rho in descending order
[~,ordrho]=sort(rho,'descend');

% the index of points with a larger rho
largerRhoIdx=[];

maxRhoIdx=ordrho(1); % the maximal rho index
DeltaParent(maxRhoIdx)=maxRhoIdx; 
for i=2:N %
    thisIdx=ordrho(i);
    largerRhoIdx =ordrho(1:(i-1)); % 
    % the distance to the points with larger rho
    largerRhoDist=distMeasure(features(:,thisIdx),features(:,largerRhoIdx)); 
    % the nearest one
    [val,minone]=min(largerRhoDist); % 

    delta(thisIdx)=val; % 
    DeltaParent(thisIdx)=largerRhoIdx(minone); % 
end

% the delta value of the point with the globally maximal rho
delta(maxRhoIdx)=0;
delta(maxRhoIdx)=max(delta); % 

end

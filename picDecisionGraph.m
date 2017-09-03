function PicGraph = picDecisionGraph(rho, delta, ClusterIdx, cluster_centers, DeltaParent, titletoset, xl, yl, varargin)

% this function is to plot the decision graph according to rho and delta
% ClusterIdx is the assignment
% cluster_center is the indexes of cluster centers
% the function returns the handle of the figure

    if(isempty(rho) || isempty(delta))
        return;
    end
    
    DotOpacity = 0.2;
    LineOpacity = 0.6;
    RaibowLineWidth = 4;
    
    [ClusterIdx,cluster_centers] = rearrgeClusterIdx(rho,delta,ClusterIdx,cluster_centers);
    delta = cropMax(delta);

    PicGraph = figure;
    hold on

    if(isempty(titletoset))
        titletoset = 'Decision Graph';
    end
    title(titletoset,'FontSize',40)
    if isempty(xl)
        xl = '\rho';
    end
    if isempty(yl)
        yl = '\delta';
    end
    xlabel(xl);
    ylabel(yl);

    % preparation for the rainbow
    minRho = min(rho);
    maxRho = 1.01*max(rho);
    maxDelta = max(delta)*1.05;
    seriess = linspace(minRho,maxRho);

    numOfCluster = ceil(max(unique(ClusterIdx))); % number of clusters

    cmap = colormap(jet);
    cmap = min(cmap + 0.16,1);
    colorss = cmap( int16( ( (1:numOfCluster).*56 )/( numOfCluster ) ) , :);

    DotSize = 12;
    CenterSize = 36;
    if nargin>8
        DotSize = varargin{1};
        CenterSize = varargin{2};
    end
    if nargin>10
        colorss = varargin{3};
    end


    % all the points
    for i =  1:numOfCluster
       TheColor = colorss(i,:);
       ClusterPointsHandle = plot(rho(ClusterIdx==i),delta(ClusterIdx==i),...
        'o','MarkerSize',DotSize,'MarkerFaceColor',TheColor,'MarkerEdgeColor',TheColor);
%        alpha(ClusterPointsHandle,DotOpacity);
    end
    

    % draw the rainbow
    for i = 1:numel(cluster_centers)
       % rainbow line
       gamma = rho(cluster_centers(i))*delta(cluster_centers(i));
       XXX = seriess;
       YYY = max(gamma./seriess,0);
       tmp = find(YYY<maxDelta,1,'first');
       
       TheColor = colorss(i,:);
       RainbowHandle = plot(XXX(tmp:end),YYY(tmp:end),'Color',TheColor,'LineWidth',RaibowLineWidth);
       RainbowHandle.Color(4)=LineOpacity;
       
       
    end
    
    % draw center points
    for i=1:numel(cluster_centers)
       % centers
       TheColor = colorss(i,:);
       CenterPointsHandle = plot(rho(cluster_centers(i)),delta(cluster_centers(i)),...
        'o','MarkerSize',CenterSize,'MarkerFaceColor',TheColor,'MarkerEdgeColor','w');
%        alpha(CenterPointsHandle,DotOpacity);

    end
    
    
    
    hold off
    
end

function [NewClusterIdx,NewCluster_centers] = rearrgeClusterIdx(rho,delta,ClusterIdx,cluster_centers)
    
    NewClusterIdx = zeros(size(ClusterIdx));
    gamma = rho.*delta;
    [~,GammaSortedIndex]=sort(gamma(cluster_centers),'descend');
    NewCluster_centers=cluster_centers(GammaSortedIndex);
    
    for i1=1:numel(cluster_centers)
        Idx = ( ClusterIdx == ClusterIdx(cluster_centers(i1)) );
        NewClusterIdx(Idx) = i1;
    end
end


function outputVector = cropMax(inputVecter)
% clip the maximal value of the vector into the second largest on

    outputVector = inputVecter;
    [~,maxIdx] = max(inputVecter,[],2);
    outputVector(maxIdx) = max( inputVecter( 1, [1:(maxIdx-1), (maxIdx+1):end] ) ) + eps; % 

end
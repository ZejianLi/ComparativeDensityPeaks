%% The Density Peaks algorithm
% input
% features  : the d*N input, where d is the dimension of the sample and N is the number of sample
% K         : the predefined number of clusters
% dc        : the cut-off parameter d_c
% rhoMeasure: the function handle to calculate rho
% distMeasure:the function handle to calculate the distance
% output
% ClusterIdx: the cluster assignment on data
% rho       : the density
% delta     : the delta distance
% cluster_centers   : the indexes of cluster centers
% DeltaParent       : the indexes of \Delta
% gamma     : the gamma value
% core_halo : the flag to indicate whether the point is in the core region


function [ClusterIdx, rho, delta, cluster_centers, DeltaParent, gamma, core_halo] = dp(features, K, dc, rhoMeasure, distMeasure)
% Density Peaks algorithm
	disp('Standard DP called');

	N=size(features, 2); % num of data


	% [rho,delta,DeltaParent,gamma,ClusterIdx] = deal(zeros(1,N)); 

	disp('calculating rho...');
	rho=rhoMeasure(features,dc);

	disp('assigning delta...');
	[delta,DeltaParent,ordrho]=assignDelta(features,rho,distMeasure);

	gamma=rho.*delta+1e-6;
	[~,~,idxK]=minK(-gamma,K); % maximal K gamma values
	cluster_centers=idxK;
	DeltaParent(cluster_centers)=cluster_centers; % the cluster centers are split

	disp('assigning cluster labels...');
	ClusterIdx=assignClusterIdx(cluster_centers,ordrho,DeltaParent);


	disp('assigning the core/halo flag...')
	core_halo=false(1,N);
	 rho_bs=zeros(1,K);
	for i=1:K
	    diffDist=distMeasure(features(:,ClusterIdx~=i),features(:,ClusterIdx==i));
	    border_distance=(diffDist<dc); % 
	    border_area = sum(border_distance,1)>0; % 
	    rho_b=max(rho(ClusterIdx==i).*border_area); % 
	    rho_bs(i)=rho_b;
	end
	rho_b=max(rho_bs); % 
	core_halo(rho>=rho_b)=true; % 

end
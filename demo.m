% This script load the Optdigits first
% then it gives the best performance evaluation with dc from 0.1 to 1
% finally the decision graph when dc=1 is displayed

DataSetName = 'Optdigits';
load(['data/',DataSetName]);

features = double(features);
features = normc(features);

[labels, num_labels] = utility_classAdjust(labels); 

K = num_labels;
dc_range = 0.1:0.1:1;

evaluation_DP = zeros(1, 2);
evaluation_CDP = zeros(1, 2);

for i1 = 1:numel(dc_range)
    
	dc = dc_range(i1);
    
    [ClusterIdx, rho, delta, cluster_centers, DeltaParent] = dp(features, K, dc, @gaussianRhoMeasure, @l2distMeasure);
    
	evaluations = clusteringEvaluation(ClusterIdx, labels);
    evaluation_DP = max(evaluation_DP, evaluations(1, [1,3]));
    
    
    if dc==1
	    PicGraph = picDecisionGraph(rho, delta, ClusterIdx, cluster_centers, DeltaParent, ['DP d_c=',num2str(dc)], '\rho', '\delta');
        set(gca,'XLim',[min(rho),max(rho)*1.05]);
        set(gca,'YLim',[0,max(delta*1.01)]);
    end

	[ClusterIdx, rho, delta, cluster_centers, DeltaParent, ~, theta] = ComparativeDensityPeaks(features, K, dc, DataSetName);

	evaluations = clusteringEvaluation(ClusterIdx, labels);
    evaluation_CDP = max(evaluation_CDP, evaluations(1, [1,3]));
    
    
    if dc==1
	    PicGraph = picDecisionGraph(rho, theta, ClusterIdx, cluster_centers, DeltaParent, ['CDP d_c=',num2str(dc)], '\rho', '\theta');
	    set(gca,'XLim',[min(rho),max(rho)*1.05]);
        set(gca,'YLim',[0,max(theta*1.01)]);
	end

    
end

disp()
disp(['DP ACC and NMI values: ' , num2str(evaluation_DP) ]);
disp(['CDP ACC and NMI values: ' , num2str(evaluation_CDP) ]);

    
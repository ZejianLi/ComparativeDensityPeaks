function [new_labels,num_labels,new_LabelValues] = utility_classAdjust(labels)
% adjust the class labels from arbitrary integer to 1:num_labels
% written by Zejian Li

	original_LabelValues = unique(labels);
	num_labels = length(original_LabelValues);

	if all(unique(labels)==[1:num_labels])
	    new_labels = labels;
	    new_LabelValues = original_LabelValues;
	else
	    new_labels = zeros(size(labels));    
	    for i = 1:num_labels
	        new_labels(labels==original_LabelValues(i))=i;
	    end    
	    new_LabelValues = 1:num_labels;
	end

end
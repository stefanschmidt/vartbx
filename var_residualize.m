function residualized_scans = var_residualize(motion_parameters, original_scans)
%
% variability toolbox: residualize scans
%__________________________________________________________________________
%
% residualizes a set of scans (original_scans) on a set of
% predictors (motion_parameters) using linear regression
%

  % define intercept term required for linear regression
	number_of_data_rows = size(original_scans, 1);
	intercept_terms = ones(number_of_data_rows, 1);

	% perform linear regression on motion parameters and scans
	regression_coefficients = [intercept_terms motion_parameters] \ original_scans;
	
	% predict scans with linear model thus generating
	% regression lines for all dimensions of the scans
	predictions = [intercept_terms motion_parameters] * regression_coefficients;

	% remove correlations between the motion parameters and
	% the scans by subtracting the regression lines from all
	% dimensions of the original scans
	residualized_scans = original_scans - predictions;

end

function output = tbx_run_variability(job)

	% Variability Toolbox
	% input  (job): job data structure from toolbox configuration dialog
	% output (out): Full paths to the computed variability volumes

	job_result_files = {};
	condition_scans = {};
	conditions = {job.run(1).condition.name};
	num_conditions = numel(conditions);
	mask_file = job.mask{1};
	subject_header = spm_vol(job.run(1).scans{1});
	image_dimensions = subject_header.dim;

	if isempty(mask_file)
		common_coordinates = 1:prod(image_dimensions);
	else
		mask_header = spm_vol(mask_file);
		headers = [subject_header, mask_header];
		if spm_check_orientations(headers) == 0
			error('Error: Scans differ in orientation, dimension or voxel size.');
		end
		mask_volume = spm_read_vols(mask_header);
		common_coordinates = find(mask_volume);
	end

	header = subject_header;
	header.descrip = 'Variability Toolbox 0.1b';	

	% scan_indices{condition}{run}{block} -> scan indices for a condition block
	% condition_data{condition}           -> numel(common_coordinates) x num_scans_condition
	% condition_data gets initialized with zero matrices
	[scan_indices, condition_data] = order_blocks_by_condition();

	for condition = 1:numel(conditions)
		count{condition} = 0;
	end

	% initialize window for graphical progress indicator
	% avoid pixel shifts happening with default renderer
	progress_window = spm_figure('GetWin','Interactive');
	set(progress_window, 'Renderer','OpenGL');

	for run = 1:numel(job.run)
		residualized_data = residualize_scans_by_motion(run);
		merge_blocks_by_condition(run, residualized_data);
	end

	% variability computation
	%
	% condition_scans{condition} -> voxels x scans
	% the variability computation is performed row-wise
	% and the scan matrix is collapsed to a single scan
	% e.g. 50000x20 => 50000x1
	var_data = zeros(numel(conditions), numel(common_coordinates));
	for condition = 1:numel(conditions)
		if strcmp(job.metric, 'detrended_sd')
			var_data(condition,:) = squeeze(std(condition_scans{condition}, 0, 2))';
		elseif strcmp(job.metric, 'mssd')
			var_data(condition,:) = mssd(condition_scans{condition});
		else
			error('Error: Invalid variability measure.');
		end
	end

	% save results
	%
	label = '\fontsize{16}Saving Results';
	spm_progress_bar('Init', 100, label, '', 't');
	result_image = zeros(image_dimensions);

	cd(job.resultdir{1})

	for condition=1:numel(conditions)
		result_image(common_coordinates) = var_data(condition,:);
		header.fname = [job.resultprefix '_' conditions{condition} '.nii'];
		header.dt = [spm_type('float32') spm_platform('bigend')];
		spm_write_vol(header, result_image);
		job_result_files{end+1} = header.fname;
		progress = round(100 * condition / numel(conditions));
		spm_progress_bar('Set', progress);
	end

	output.result_files = job_result_files;

	function [scan_indices, condition_data] = order_blocks_by_condition()

		for condition = 1:numel(conditions)
			num_scans_condition = 0;

			for run = 1:numel(job.run)
				onsets = job.run(run).condition(condition).onset;
				durations = job.run(run).condition(condition).duration;

				if strcmp(job.timing.units, 'seconds')
					[onsets, durations] = secs_to_scans(onsets, durations, job.timing.RT);
				else
					% Expect the user to specify the first scan as zero
					% regardless if the times unit is scans or seconds
					onsets = onsets + 1;
					if not(positive_ints(onsets) && positive_ints(durations))
						error('Error: Onsets and duration must be non-negative integers for unit "scans".')
					end
				end

				% if duration is defined only once extend it to a vector
				if numel(durations) == 1
					durations = ones(1, numel(onsets)) * durations;
				end

				for block = 1:numel(onsets)
					scan_indices{condition}{run}{block} = onsets(block) - 1 + [1:durations(block)];
					this_duration = numel(scan_indices{condition}{run}{block});
					num_scans_condition = num_scans_condition + this_duration;
				end
			end
			condition_data{condition} = zeros(numel(common_coordinates), num_scans_condition);
		end

	end

	function img = residualize_scans_by_motion(run)

		scan_files = job.run(run).scans;
		num_scans = numel(scan_files);
		motion_parameters_file = char(job.run(run).residualize);
		headers = cell2mat(spm_vol(scan_files));

		label_text = sprintf('Reading Run %i', run);
		label_format = '\fontsize{16}';
		label = [label_format label_text];
		spm_progress_bar('Init', 100, label, '', 't');

		volumes = zeros([headers(1).dim, num_scans]);

		for v = 1:num_scans
			volumes(:,:,:,v) = spm_read_vols(headers(v));
			spm_progress_bar('Set', round(100*v/num_scans));
		end

		img = double(reshape(volumes, [], size(volumes, 4)));
		clear volumes;
		img = img(common_coordinates, :);

		% residualize if parameter file is given
		if not(isempty(motion_parameters_file))
			label_text = sprintf('Processing Run %i', run);
			label_format = '\fontsize{16}';
			label = [label_format label_text];
			spm_progress_bar('Init', 100, label, '', 't');

			temporal_mean = mean(img, 2);
			mp_time_series = load(motion_parameters_file);

			batch_size = 5000;
			num_voxels = numel(common_coordinates);
			remaining_voxels = mod(num_voxels, batch_size) - 1;
			batching_end = num_voxels - remaining_voxels - batch_size;

			for voxel = 1:batch_size:batching_end
				voxel_range = voxel:(voxel + batch_size - 1);
				img(voxel_range,:) = var_residualize([mp_time_series], img(voxel_range,:)')' + repmat(temporal_mean(voxel_range), [1 size(img,2)]);
				spm_progress_bar('Set', round(99*voxel/num_voxels));
			end

			if remaining_voxels > 0
				voxel_range = (num_voxels - remaining_voxels):num_voxels;
				img(voxel_range,:) = var_residualize([mp_time_series], img(voxel_range,:)')' + repmat(temporal_mean(voxel_range), [1 size(img,2)]);
			end

			spm_progress_bar('Set', 100);
		end

	end

	function merge_blocks_by_condition(run, scan_data)

		% only do computations over the selected voxels
		voxel_indices = 1:numel(common_coordinates);
		
		for condition = 1:numel(conditions)
			for block = 1:numel(scan_indices{condition}{run})

				block_data = scan_data(:, scan_indices{condition}{run}{block});

				% normalize block_data to global block mean = 100
				block_data = 100 * block_data / mean(mean(block_data));

				if strcmp(job.metric, 'detrended_sd')

					% temporal mean of this block
					block_mean = mean(block_data, 2);

					% subtract temporal mean to block-wise detrend the data
					for scan = 1:size(block_data, 2)
						count{condition} = count{condition} + 1;
						adjusted_block_data = block_data(voxel_indices, scan) - block_mean(voxel_indices);
						condition_scans{condition}(voxel_indices, count{condition}) = adjusted_block_data;
					end
				elseif strcmp(job.metric, 'mssd')
					for scan = 1:size(block_data, 2)
						count{condition} = count{condition} + 1;
						condition_scans{condition}(voxel_indices, count{condition}) = block_data(voxel_indices, scan);
					end
				else
					error('Error: Invalid variability measure.');
				end
			end
		end
	end

	% compute the compute mean square successive difference (MSSD) for a single
	% condition over all scans assigned to that condition of all runs
	function mssd_scan = mssd(scans)
		mssd_scan = sum((diff(scans, 1, 2).^2) / (length(scans) - 1), 2);
	end

	% convert the unit of a condition from seconds to scans
	function [onsets, durations] = secs_to_scans(secs_onsets, secs_durations, RT)

		onsets_valid = sum(mod(secs_onsets, RT)) == 0;
		durations_valid = sum(mod(secs_durations, RT)) == 0;

		if ~(onsets_valid && durations_valid)
			message = 'Onsets and/or durations not a multiple of RT,\n';
			message = 'rounding result to nearest scan number.';
			warning('%s', message);
		end

		% convert definitions in seconds to scan index numbers
		onsets = round((secs_onsets / RT) + 1);
		durations = round(secs_durations / RT);

	end

	function valid = positive_ints(vector)
		valid = (sum(vector >= 1) == length(vector)) && (sum(mod(vector,1)) == 0);
	end

	%% scaffold for integration of a cli progress bar for the grid
	% function set_progress_mode(run_mode)

		% if strcmp(run_mode, 'gui')
			% @progress_bar = spm_progress_bar
		% else
			% @progress_bar = cli_progress_bar
		% end

end

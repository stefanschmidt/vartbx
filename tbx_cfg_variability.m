function variability = tbx_cfg_variability
% Variability Toolbox Configuration
%__________________________________________________________________________

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Variability')); end

%--------------------------------------------------------------------------
% metric Varibility Metric
%--------------------------------------------------------------------------

metric        = cfg_menu;
metric.tag    = 'metric';
metric.name   = 'Variability function';
metric.help   = {'Select a variability function to apply to the scans'};
metric.labels = {'Detrended SD' 'MSSD'}';
metric.values = {'detrended_sd' 'mssd'};
metric.val 		= {'detrended_sd'};

% ---------------------------------------------------------------------
% resultprefix Result prefix
% ---------------------------------------------------------------------

resultprefix         = cfg_entry;
resultprefix.tag     = 'resultprefix';
resultprefix.name    = 'Filename prefix';
resultprefix.help    = {'A filename prefix for the results'};
resultprefix.strtype = 's';
resultprefix.val     = {'var'};
resultprefix.num     = [0 Inf];

% ---------------------------------------------------------------------
% resultdir Result directory
% ---------------------------------------------------------------------

resultdir         = cfg_files;
resultdir.tag     = 'resultdir';
resultdir.name    = 'Result directory';
resultdir.help    = {'Select a directory where the variability images will be written.'
									 'By default they will be written to the current working directory.'};
resultdir.filter  = 'dir';
resultdir.ufilter = '.*';
resultdir.num     = [1 1];

% ---------------------------------------------------------------------
% units Units for design
% ---------------------------------------------------------------------

units         = cfg_menu;
units.tag     = 'units';
units.name    = 'Units for design';
units.help    = {'The onsets of events or blocks can be specified in either scans or seconds.'};
units.labels = {
                'Scans'
                'Seconds'
}';
units.values = {
                'scans'
                'seconds'
}';
% ---------------------------------------------------------------------
% RT Interscan interval
% ---------------------------------------------------------------------

RT         = cfg_entry;
RT.tag     = 'RT';
RT.name    = 'Interscan interval';
RT.help    = {'Interscan interval, TR, (specified in seconds).  This is the time between acquiring a plane of one volume and the same plane in the next volume.  It is assumed to be constant throughout.'};
RT.strtype = 'e';
RT.num     = [1 1];

% ---------------------------------------------------------------------
% timing Timing parameters
% ---------------------------------------------------------------------

timing         = cfg_branch;
timing.tag     = 'timing';
timing.name    = 'Timing parameters';
timing.val     = {units RT};
timing.help    = {
'Specify various timing parameters needed to construct the design matrix. This includes the units of the design specification and the interscan interval.'
''
'Also, with longs TRs you may want to shift the regressors so that they are aligned to a particular slice.  This is effected by changing the microtime resolution and onset. '
}';

% ---------------------------------------------------------------------
% run_name Run Name
% ---------------------------------------------------------------------

run_name         = cfg_entry;
run_name.tag     = 'run_name';
run_name.name    = 'Name';
run_name.help    = {'A distinctive name for the run (e.g. the run number)'};
run_name.strtype = 's';
run_name.val     = {''};
run_name.num     = [0 Inf];

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------

scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.check   = @scans_check;
scans.help    = {'Select the fMRI scans for this run.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];

% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------

name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Condition Name'};
name.strtype = 's';
name.num     = [1 Inf];

% ---------------------------------------------------------------------
% onset Onsets
% ---------------------------------------------------------------------

onset         = cfg_entry;
onset.tag     = 'onset';
onset.name    = 'Onsets';
onset.help    = {'Specify a vector of onset times for this condition type.'};
onset.strtype = 'e';
onset.num     = [Inf 1];

% ---------------------------------------------------------------------
% duration Durations
% ---------------------------------------------------------------------

duration         = cfg_entry;
duration.tag     = 'duration';
duration.name    = 'Durations';
duration.help    = {'Specify the event durations. Epoch and event-related responses are modeled in exactly the same way but by specifying their different durations.  Events are specified with a duration of 0.  If you enter a single number for the durations it will be assumed that all trials conform to this duration. If you have multiple different durations, then the number must match the number of onset times.'};
duration.strtype = 'e';
duration.num     = [Inf 1];

% ---------------------------------------------------------------------
% cond Condition
% ---------------------------------------------------------------------
condition         = cfg_branch;
condition.tag     = 'condition';
condition.name    = 'Condition';
condition.val     = {name onset duration};
condition.check   = @condition_check;
condition.help    = {'An array of input functions is contructed, specifying occurrence events or epochs (or both). These are convolved with a basis set at a later stage to give regressors that enter into the design matrix. Interactions of evoked responses with some parameter (time or a specified variate) enter at this stage as additional columns in the design matrix with each trial multiplied by the [expansion of the] trial-specific parameter. The 0th order expansion is simply the main effect in the first column.'};

% ---------------------------------------------------------------------
% generic Conditions
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Conditions';
generic1.help    = {'You are allowed to combine both event- and epoch-related responses in the same model and/or regressor. Any number of condition (event or epoch) types can be specified.  Epoch and event-related responses are modeled in exactly the same way by specifying their onsets [in terms of onset times] and their durations.  Events are specified with a duration of 0.  If you enter a single number for the durations it will be assumed that all trials conform to this duration.For factorial designs, one can later associate these experimental conditions with the appropriate levels of experimental factors. '};
generic1.values  = {condition };
generic1.num     = [1 Inf];

% ---------------------------------------------------------------------
% residualize Residualization
% ---------------------------------------------------------------------
residualize         = cfg_files;
residualize.tag     = 'residualize';
residualize.name    = 'Residualization';
residualize.help    = {'File with covariates of no interest that you wish to residualize from the data (e.g. motion parameters). The file is assumed to be column-based where the columns are what will be removed from the data.'};
residualize.filter  = 'any';
residualize.ufilter = '.*';
residualize.val     = {{''}};
residualize.num     = [0 1];

% ---------------------------------------------------------------------
% separator Separator
% ---------------------------------------------------------------------

separator      = cfg_const;
separator.name = '';
separator.tag  = 'separator';
separator.val  = {1};
separator.help = {['']};

% ---------------------------------------------------------------------
% fmri_run run
% ---------------------------------------------------------------------
fmri_run    		 = cfg_branch;
fmri_run.tag     = 'run';
fmri_run.name    = 'Run';
fmri_run.val     = {run_name scans generic1 residualize separator};
fmri_run.help    = {''};


% ---------------------------------------------------------------------
% generic Conditions & Onsets
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Data & Design';
generic.help    = {'Compute variability condition-wise across runs. Specifying multiple runs will result in one image per condition.'};
generic.values  = {fmri_run};
generic.num     = [1 Inf];

% ---------------------------------------------------------------------
% mask Mask
% ---------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Mask';
mask.val     = {{''}};
mask.help    = {'Specify an image for masking the analysis (e.g. gray matter mask).'};
mask.filter = 'image';
mask.ufilter = '.*';
mask.num     = [0 1];

%--------------------------------------------------------------------------
% variability Calculation of variability in series of volumes
%--------------------------------------------------------------------------
variability				= cfg_exbranch;
variability.tag		= 'variability';
variability.name	= 'Variability';
variability.val		= {metric timing mask resultprefix resultdir generic};
variability.help	= {'This toolbox is intended to be similar to first level analysis. It measures voxel-based variability and will output Nifti files. The assumption is to proceed with those to level two analysis to model effective interest, however the output could also be used with other statistics programs.'};
variability.prog	= @tbx_run_variability;
variability.vout	= @vout;

%--------------------------------------------------------------------------
% validation of input data
%--------------------------------------------------------------------------

function t = scans_check(scans)

	t = {};

	for idx = 1:numel(scans)

		[pth, nam, ext, ~] = spm_fileparts(deblank(scans{idx}));
		p = fullfile(pth,[nam ext]);

		if ~spm_existfile(p)
			t = {sprintf('File "%s" does not exist.', p)};
		end

		switch ext
	    case {'.nii','.NII'}
				% Do nothing

	    case {'.img','.IMG'}
				if ~spm_existfile(fullfile(pth,[nam '.hdr'])) && ...
					 ~spm_existfile(fullfile(pth,[nam '.HDR']))
					t = {sprintf('File "%s" does not exist.', fullfile(pth,[nam '.hdr']))};
				end
	        
	    case {'.hdr','.HDR'}
				ext = '.img';
				p = fullfile(pth,[nam ext]);
				if ~spm_existfile(p)
					t = {sprintf('File "%s" does not exist.', p)};
				end

			otherwise
				t = {sprintf('File "%s" is not of a recognised type.', p)};
		end % switch
	end % for
end % function

function t = condition_check(item)
	t = {};
	if (numel(item.onset) ~= numel(item.duration)) && (numel(item.duration)~=1),
		t = {sprintf('Number of event onsets (%d) does not match the number of durations (%d) for condition %s.', numel(item.onset),numel(item.duration),item.name)};
	end
end

function dep = vout(varargin)
	dep(1)            = cfg_dep;
	dep(1).sname      = 'Variability in series of volumes';
	dep(1).src_output = substruct('.','variability');
	dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

end
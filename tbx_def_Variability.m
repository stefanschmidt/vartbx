function tbx_def_Variability

variability.metric = 'detrended_sd';
variability.timing.units = 'scans';
variability.timing.RT = 1;

matlabbatch{1}.spm.tools.variability = variability;

% Set Modality
spm('defaults', 'FMRI')

% Run the job interactively (GUI)
spm_jobman('interactive', matlabbatch)

%end
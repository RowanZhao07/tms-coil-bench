function run_external_measurements
% Reanalyse eight published measurement campaigns after running prepare_inputs.py.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
startPath = path;
restorePath = onCleanup(@() path(startPath)); %#ok<NASGU>
addpath(genpath(here));
bench = fullfile(root,'placement_pilot','public_benchmark');
addpath(bench,fullfile(bench,'matched_figure8'));

validate_petrov_figure8;
plot_petrov_comparison;
rescore_petrov;
replay_pose_check;
replay_point_sensitivity;

run_mc_b70;
run_waveforms;
plot_lab_validation;
run_hall_distance;
run_lyon_distance;
rerun_abbasi;
run_italy;
run_finland;
run_roth_profiles;
fprintf('Completed all eight measurement campaigns.\n');
end

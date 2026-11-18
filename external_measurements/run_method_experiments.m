function run_method_experiments
% From the project root, run prepare_inputs.py once for the original scan.
% Reconstruct the two upstream physical inputs when starting a clean checkout.
here=fileparts(mfilename('fullpath'));
root=fileparts(here);
previous=path; cleanup=onCleanup(@() path(previous)); %#ok<NASGU>
addpath(genpath(here),fullfile(root,'placement_pilot','public_benchmark','matched_figure8'));
if ~isfile(fullfile(root,'placement_pilot','public_benchmark','matched_figure8','figure8_slice_results.mat'))
    validate_petrov_figure8;
end
if ~isfile(fullfile(here,'measured_coil_library','mc_b70_equivalent_sources.csv'))
    run_mc_b70;
end
run_pose_refinement;
run_sensor_design;
plot_sensor_design;
run_decay_calibration;
fprintf('Completed the three calibration and ablation experiments.\n');
end

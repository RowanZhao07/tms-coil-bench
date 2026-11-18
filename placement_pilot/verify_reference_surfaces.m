function checks = verify_reference_surfaces()
% Recompute spherical and tangent-plane results with the original MATLAB kernel.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'matlab'));
ref = readtable(fullfile(here,'reference_surface_comparison.csv'),'TextType','string');
checks = ref(:,{'coil','depth_mm'});
names = unique(ref.coil);
for k = 1:numel(names)
    coil = names(k);
    rows = ref.coil == coil;
    s = readmatrix(fullfile(here,'placement','target_frame_segments',coil+".csv"));
    delta = unique(ref.centerline_gap_mm(rows)) - min(s(:,[3,6]),[],'all');
    assert(max(abs(ref.shift_to_flat_mm(rows)-delta)) < 1e-9);
    s(:,1:6) = s(:,1:6)/1000;
    targets = [zeros(nnz(rows),2), -ref.depth_mm(rows)]/1000;
    sphere = field_segments(s,targets,1000)*1000;
    s(:,[3,6]) = s(:,[3,6])+delta/1000;
    flat = field_segments(s,targets,1000)*1000;
    checks.sphere_norm_mT(rows) = vecnorm(sphere,2,2);
    checks.flat_norm_mT(rows) = vecnorm(flat,2,2);
    sphere_ref = ref{rows,{'sphere_Bx_mT','sphere_By_mT','sphere_Bz_mT'}};
    flat_ref = ref{rows,{'flat_Bx_mT','flat_By_mT','flat_Bz_mT'}};
    checks.sphere_vector_difference_mT(rows) = vecnorm(sphere-sphere_ref,2,2);
    checks.flat_vector_difference_mT(rows) = vecnorm(flat-flat_ref,2,2);
end
assert(max(checks{:,5:6},[],'all') < 1e-7);
writetable(checks,fullfile(here,'matlab_reference_check.csv'));
disp(checks)
end

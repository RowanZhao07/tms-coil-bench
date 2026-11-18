function [placements,aggregate]=replay_point_sensitivity()
% Translate the entire six-point sampling grid within the same known MRI event.
% Both candidate predictions are frozen at the published current. Lower (unscaled)
% six-point L2 wins. No far-side observation enters this decision.
here=fileparts(mfilename('fullpath'));
source=fullfile(here,'..','..','placement_pilot','public_benchmark','matched_figure8','figure8_slice_results.mat');
d=load(source,'mask','conditions','prediction_T','measured_T');
assert(isequal(d.conditions(1:2),["v1.0","v1.1"]));
old=1e6*d.prediction_T(:,:,1); fixed=1e6*d.prediction_T(:,:,2);
measured=1e6*d.measured_T;
[base_r,base_c]=ndgrid([55 65 75],[60 105]);
offsets=-4:2:4;
far=d.mask;far(1:89,:)=false;
placements=table(); invalid=0;
for dr=offsets
    for dc=offsets
        id=sub2ind(size(d.mask),base_r(:)+dr,base_c(:)+dc);
        if ~all(d.mask(id)) || any(far(id))
            invalid=invalid+1;
            continue
        end
        y=measured(id);p0=old(id);p1=fixed(id);
        assert(all(isfinite([y;p0;p1])));
        r0=norm(p0-y);r1=norm(p1-y);
        winner="tie";
        if r1<r0, winner="v1.1"; elseif r1>r0, winner="v1.0"; end
        % Euclidean distance to the equal-squared-error decision hyperplane.
        % The reported RMS units allow direct comparison with per-probe noise.
        noise_to_tie=(r0^2-r1^2)/(2*norm(p1-p0)*sqrt(numel(id)));
        s=table(dr,dc,winner,100*r0/norm(y),100*r1/norm(y),...
            100*(r0-r1)/norm(y),(r0-r1)/sqrt(numel(id)),noise_to_tie, ...
            'VariableNames',{'row_shift_px','column_shift_px','selected','old_L2_percent',...
            'fixed_L2_percent','old_minus_fixed_L2_pp','old_minus_fixed_RMSE_uT','noise_to_tie_RMS_uT'});
        placements=[placements;s]; %#ok<AGROW>
    end
end
assert(height(placements)+invalid==numel(offsets)^2);
aggregate=table(numel(offsets)^2,height(placements),invalid,nnz(placements.selected=="v1.1"), ...
    nnz(far),min(placements.old_minus_fixed_L2_pp),min(placements.old_minus_fixed_RMSE_uT),...
    min(placements.noise_to_tie_RMS_uT), ...
    'VariableNames',{'requested_grids','valid_grids','invalid_grids','fixed_pose_selected','far_pixels',...
    'minimum_L2_margin_pp','minimum_RMSE_margin_uT','minimum_noise_to_tie_RMS_uT'});
writetable(placements,fullfile(here,'translated_grid_scores.csv'));
writetable(aggregate,fullfile(here,'translated_grid_summary.csv'));
disp(aggregate)
end

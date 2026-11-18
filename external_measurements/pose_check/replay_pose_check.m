function scores = replay_pose_check()
% Retrospective sparse-point registration check on the Petrov BSM-819 Bz slice.
here = fileparts(mfilename('fullpath'));
source = fullfile(here,'..','..','placement_pilot','public_benchmark','matched_figure8');
d = load(fullfile(source,'figure8_slice_results.mat'));
assert(isequal(d.conditions(1:2),["v1.0","v1.1"]));
[row,col] = ndgrid([55 65 75],[60 105]);
points = sub2ind(size(d.mask),row(:),col(:));
assert(all(d.mask(points)));
% Leave entire lower part of the mask out of both gain fitting and pose selection.
test = d.mask;
test(1:89,:) = false;
assert(~any(test(points)));
measured = 1e6*d.measured_T;
readings = table(row(:),col(:),measured(points),'VariableNames',{'image_row','image_column','MRI_uT'});
scores = table();
for v = 1:2
    field = 1e6*d.prediction_T(:,:,v);
    fitFixed = field(points); reference = measured(points);
    gain = dot(fitFixed,reference)/dot(fitFixed,fitFixed);
    readings.("predicted_"+replace(d.conditions(v),".","_")+"_uT") = fitFixed;
    for region=["six_points","held_out_rows_90_end"]
        if region=="six_points", sim=fitFixed; truth=reference;
        else, sim=field(test); truth=measured(test); end
        rowScore = table(d.conditions(v),region,numel(sim),gain, ...
            100*norm(sim-truth)/norm(truth),100*norm(gain*sim-truth)/norm(truth),sqrt(mean((sim-truth).^2)), ...
            'VariableNames',{'pose','region','pixels','gain_from_six_points','fixed_L2_percent','gain_corrected_L2_percent','fixed_RMSE_uT'});
        scores=[scores;rowScore]; %#ok<AGROW>
    end
end
writetable(readings,fullfile(here,'six_points.csv'));
writetable(scores,fullfile(here,'pose_holdout_scores.csv'));
fig=figure('Visible','off','Color','white','Position',[80 80 1020 445]);
t=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(t); h=imagesc(ax,measured); h.AlphaData=double(d.mask);
set(ax,'Color',[.94 .94 .94]); colormap(ax,parula); c=colorbar(ax); c.Label.String='B_z (µT)';
axis(ax,'image'); hold(ax,'on'); plot(ax,col(:),row(:),'wo','LineWidth',1.8,'MarkerSize',8);
yline(ax,89.5,'k--','LineWidth',1.3); xlabel(ax,'像素列'); ylabel(ax,'像素行');
title(ax,'六个检测点与留出的远侧区域');
ax=nexttile(t); far=scores(scores.region=="held_out_rows_90_end",:);
bar(ax,categorical(far.pose),[far.fixed_L2_percent,far.gain_corrected_L2_percent]); grid(ax,'on');
legend(ax,{'按公开电流','仅用六点估计增益'},'Location','northeast','Box','off');
ylabel(ax,'留出区域相对 L2 误差 (%)'); title(ax,'同一张实测切片的异地测点检验');
set(findall(fig,'-property','FontName'),'FontName','Microsoft YaHei');
exportgraphics(fig,fullfile(here,'pose_holdout.png'),'Resolution',160); close(fig);
disp(scores)
end

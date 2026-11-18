function plot_validation
% Plot off-axis directional errors and convergence from completed 3-D solves.
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results');
t=readtable(fullfile(out,'mesh_domain_convergence.csv'));
p=readtable(fullfile(out,'final_roi.csv'));
f=figure('Visible','off','Color','white','Position',[100 100 1080 440]);
layout=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout); hold(ax,'on');
for width=[256 320 480]
    rows=t(t.half_width_mm==width,:);
    plot(ax,rows.h_mm,rows.relative_vector_l2_pct,'o-','LineWidth',1.5,...
        'DisplayName',sprintf('Half-width %d mm',width));
end
set(ax,'XDir','reverse');xticks(ax,[1.5 2 4 8]);grid(ax,'on');
xlabel(ax,'Hexahedral mesh spacing (mm)');ylabel(ax,'Vector relative L2 difference (%)');
title(ax,'Mesh and air-domain sensitivity');legend(ax,'Location','northwest');
ax=nexttile(layout);
slice=p(abs(p.z_m-.02)<1e-10,:);
xy=unique(slice.x_m)*1000;
assert(height(slice)==numel(xy)^2);
slice=sortrows(slice,{'x_m','y_m'});
imagesc(ax,xy,xy,reshape(slice.relative_vector_error_pct,numel(xy),numel(xy)));
axis(ax,'xy','image'); cb=colorbar(ax);cb.Label.String='Pointwise vector difference (%)';
xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');title(ax,'Difference from MATLAB at z = 20 mm');
title(layout,'Double-cone coil in air | 3-D Galerkin FEM');
exportgraphics(f,fullfile(out,'validation.png'),'Resolution',160);
close(f);
end

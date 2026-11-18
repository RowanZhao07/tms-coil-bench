function plot_placement
% Reproduce the nominal registration overview from the exported source paths.
root=fileparts(mfilename('fullpath'));
ids={'h1','figure8','double_cone'};
names={'H1 | reconstructed','Figure-eight | idealized','Double cone | parameterized'};
fig=figure('Visible','off','Color','w','Position',[20 20 1320 480]);
lay=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for k=1:numel(ids)
    ax=nexttile(lay);
    s=readmatrix(fullfile(root,'target_frame_segments',[ids{k} '.csv']));
    [sx,sy,sz]=sphere(48);
    surf(ax,85*sx,85*sy,85*sz-85,'FaceColor',[.52 .63 .69], ...
        'FaceAlpha',.17,'EdgeColor','none');hold(ax,'on');
    v=linspace(0,2*pi,240);
    plot3(ax,85*cos(v),zeros(size(v)),85*sin(v)-85,'Color',[.64 .69 .72]);
    plot3(ax,85*cos(v),85*sin(v),zeros(size(v))-85,'Color',[.64 .69 .72]);
    xyz=cell(1,3);
    for a=1:3, b=[s(:,a),s(:,a+3),nan(size(s,1),1)]'; xyz{a}=b(:); end
    plot3(ax,xyz{1},xyz{2},xyz{3},'Color',[.03 .34 .47],'LineWidth',.65);
    plot3(ax,[0 0],[0 0],[0 -50],'Color',[.72 .19 .12],'LineWidth',1.5);
    scatter3(ax,[0 0 0],[0 0 0],[-10 -30 -50],26,[.72 .19 .12],'filled');
    scatter3(ax,0,0,0,37,'o','MarkerEdgeColor',[.72 .19 .12],'LineWidth',1.3);
    daspect(ax,[1 1 1]);view(ax,0,12);
    xlim(ax,[-135 135]);ylim(ax,[-135 135]);zlim(ax,[-175 95]);
    grid(ax,'on');ax.FontSize=9;
    xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');zlabel(ax,'z (mm)');title(ax,names{k});
end
title(lay,'Equal 12 mm winding clearance on an 85 mm sphere | apex and 10, 30, 50 mm targets in red');
exportgraphics(fig,fullfile(root,'nominal_registration.png'),'Resolution',160);
close(fig)
end

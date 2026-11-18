function plot_sensor_design
% Visualise the recorded six-probe layouts and their independent plane scores.
here=fileparts(mfilename('fullpath'));
sites=readtable(fullfile(here,'sensor_locations.csv'),'TextType','string');
s=readtable(fullfile(here,'layout_summary.csv'),'TextType','string');
s=s(s.fit=="Joint",:);
f=figure('Visible','off','Color','w','Position',[80 60 1120 480]);
t=tiledlayout(f,1,2,'Padding','compact','TileSpacing','loose');
ax=nexttile(t);
hold(ax,'on');
[xx,yy]=ndgrid(-56:7:56);
scatter(ax,xx(:),yy(:),7,[.80 .82 .83],'filled','HandleVisibility','off');
names=["Fixed6","Dopt6","Weighted6"]; labels=["Symmetric","D-optimal","Residual-weighted"];
colors=[.24 .28 .32;.83 .35 .10;0 .43 .42]; shapes=["s","^","d"];
for i=1:3
    sel=sites.layout==names(i);
    scatter(ax,sites.x_mm(sel),sites.y_mm(sel),72,colors(i,:),shapes(i),'LineWidth',1.5, ...
        'DisplayName',labels(i));
end
axis(ax,'equal'); xlim(ax,[-64 64]); ylim(ax,[-64 64]); xticks(ax,-56:28:56); yticks(ax,-56:28:56);
xlabel(ax,'Scanner x (mm)'); ylabel(ax,'Scanner y (mm)'); title(ax,'a   Six probe locations','FontWeight','normal');
legend(ax,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',9);
ax=nexttile(t); hold(ax,'on');
groups=["Random","Shuffled"]; cloudColors=[.66 .68 .72;.56 .46 .68];
cloudNames=["Random layouts","Shuffled weights"];
for i=1:2
    q=s(startsWith(s.layout,groups(i)),:);
    scatter(ax,q.mean_position_error_mm,q.mean_test_L2_pct,19,cloudColors(i,:),'filled', ...
        'MarkerFaceAlpha',.65,'DisplayName',cloudNames(i));
end
for i=1:3
    q=s(s.layout==names(i),:);
    scatter(ax,q.mean_position_error_mm,q.mean_test_L2_pct,76,colors(i,:),shapes(i),'filled', ...
        'MarkerEdgeColor',[.17 .20 .23],'HandleVisibility','off');
end
xlabel(ax,'Mean in-plane position error (mm)'); ylabel(ax,'Two test planes: relative field error (%)');
title(ax,'b   Position and field prediction','FontWeight','normal');
xlim(ax,[0 1.1]); ylim(ax,[2.4 3.45]); grid(ax,'on'); ax.GridAlpha=.13;
legend(ax,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',9);
set(findall(f,'-property','FontName'),'FontName','Arial');
exportgraphics(f,fullfile(here,'six_probe_design.png'),'Resolution',180);
close(f);
end

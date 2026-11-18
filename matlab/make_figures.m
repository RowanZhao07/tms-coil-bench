function make_figures(root)
%MAKE_FIGURES All scientific figures are generated from executed MATLAB data.
if nargin<1, root=fileparts(fileparts(mfilename('fullpath'))); end
d=load(fullfile(root,'results','study.mat'));
out=fullfile(root,'results','figures');
set(groot,'defaultAxesFontName','Arial','defaultTextFontName','Arial',...
    'defaultAxesFontSize',10,'defaultAxesLineWidth',.8,'defaultLineLineWidth',1.5);
for k=1:numel(d.coils)
    c=d.coils(k); fig=figure('Visible','off','Color','w','Position',[80 80 1180 850]);
    tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl); draw_coil(ax,c); title(ax,'a  Directed wire centerlines');
    ax=nexttile(tl); imagesc(ax,d.x_mm,d.y_mm,d.fieldXY{k}*1e3);axis(ax,'xy','image');
    cb=colorbar(ax); cb.Label.String='|B| (mT)';colormap(ax,parula);
    xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');title(ax,'b  Field at z = 20 mm');
    ax=nexttile(tl);plot(ax,d.depth_mm,vecnorm(d.allB{k},2,2)*1e3,'Color',[.05 .35 .55]);
    grid(ax,'on');xlabel(ax,'Depth z (mm)');ylabel(ax,'|B(0,0,z)| (mT)');title(ax,'c  Field on common reference axis');
    ax=nexttile(tl);imagesc(ax,d.x_mm,d.z_mm,d.fieldXZ{k}*1e3);axis(ax,'xy');
    cb=colorbar(ax);cb.Label.String='|B| (mT)';colormap(ax,parula);
    xlabel(ax,'x (mm)');ylabel(ax,'Depth z (mm)');title(ax,'d  Field at y = 0');
    title(tl,sprintf('%02d  %s | I = 1000 A | free space',k,c.name_en),'FontSize',15,'FontWeight','bold','Interpreter','none');
    exportgraphics(fig,fullfile(out,[c.id '.png']),'Resolution',160);
    exportgraphics(fig,fullfile(out,[c.id '.pdf']),'ContentType','vector'); close(fig);
end
fig=figure('Visible','off','Color','w','Position',[30 30 1500 1050]);
tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
for k=1:numel(d.coils)
    ax=nexttile(tl);draw_coil(ax,d.coils(k));title(ax,sprintf('%02d %s',k,d.coils(k).name_en),'Interpreter','none','FontSize',10);
end
title(tl,'Fourteen parameterized coils | coordinates in mm','FontSize',17);
exportgraphics(fig,fullfile(out,'geometry_overview.png'),'Resolution',170);close(fig);
fig=figure('Visible','off','Color','w','Position',[30 30 1180 580]);tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
colors=lines(7);
for g=1:2
    ax=nexttile(tl);hold(ax,'on');
    for j=1:7
        k=(g-1)*7+j; b=vecnorm(d.allB{k},2,2);
        plot(ax,d.depth_mm,b/b(1),'Color',colors(j,:),'DisplayName',d.coils(k).name_en);
    end
    xlabel(ax,'Depth z (mm)');ylabel(ax,'|B(z)| / |B(10 mm)|');grid(ax,'on');
    legend(ax,'Location','northeast','FontSize',9,'Interpreter','none');
    title(ax,sprintf('Coils %d-%d',7*g-6,7*g));
end
title(tl,'Decay on the common reference axis | normalization is per coil','FontSize',15);
exportgraphics(fig,fullfile(out,'normalized_depth.png'),'Resolution',170);close(fig);
fig=figure('Visible','off','Color','w','Position',[30 30 1200 650]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);barh(ax,d.metrics.B30_mT,'FaceColor',[.1 .4 .6]);
yticks(ax,1:14);yticklabels(ax,string({d.coils.name_en}));set(ax,'YDir','reverse','TickLabelInterpreter','none');
xlabel(ax,'|B(0,0,30 mm)| (mT)');title(ax,'a  Equal current: 1000 A');grid(ax,'on');
ax=nexttile(tl);barh(ax,d.metrics.B30_equal_10m_mT,'FaceColor',[.15 .5 .38]);
yticks(ax,1:14);yticklabels(ax,string({d.coils.name_en}));set(ax,'YDir','reverse','TickLabelInterpreter','none');
xlabel(ax,'|B(30 mm)| at equal I L_w (mT)');title(ax,'b  Equal I L_w = 10000 A m');grid(ax,'on');
title(tl,'Amplitude comparison depends on the excitation normalization','FontSize',15);
exportgraphics(fig,fullfile(out,'amplitude_comparison.png'),'Resolution',170);close(fig);
fprintf('Exported 14 individual figure sheets and 3 overview figures.\n');
end

function draw_coil(ax,c)
s=c.segments; x=[s(:,1),s(:,4),nan(size(s,1),1)]'*1000;
y=[s(:,2),s(:,5),nan(size(s,1),1)]'*1000;z=[s(:,3),s(:,6),nan(size(s,1),1)]'*1000;
plot3(ax,x(:),y(:),z(:),'Color',[.15 .4 .55],'LineWidth',.65);hold(ax,'on');
q=max(1,round(size(s,1)/14));r=s(1:q:end,:); mid=(r(:,1:3)+r(:,4:6))/2*1000;
v=(r(:,4:6)-r(:,1:3)).*sign(r(:,7));v=v./vecnorm(v,2,2)*7;
quiver3(ax,mid(:,1),mid(:,2),mid(:,3),v(:,1),v(:,2),v(:,3),0,'Color',[.75 .2 .1],'LineWidth',.8,'MaxHeadSize',.8);
axis(ax,'equal');grid(ax,'on');view(ax,35,28);xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');zlabel(ax,'z (mm)');
zEndpoints=[s(:,3);s(:,6)];
if max(zEndpoints)-min(zEndpoints)<.02, zticks(ax,0); end
if strcmp(c.id,'double_cone'), view(ax,0,35); end
end

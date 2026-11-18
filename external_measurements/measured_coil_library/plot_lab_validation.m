function plot_lab_validation
here=fileparts(mfilename('fullpath'));
cache=fullfile(here,'..','..','tmp','measured_coil_library');
t=readtable(fullfile(here,'mc_b70_heldout_points.csv'));
t=t(abs(t.z_mm-23.8)<.2,:);
fig=figure('Visible','off','Color','w','Position',[50 50 1100 390]);
set(fig,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei');
tiles=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
b=[t.measured_Bx_uT_A t.measured_By_uT_A t.measured_Bz_uT_A];
p=[t.predicted_Bx_uT_A t.predicted_By_uT_A t.predicted_Bz_uT_A];
scale=max(abs([b(:,3);p(:,3)]));
cmap=interp1([0,.5,1],[.15,.32,.61;.98,.98,.97;.68,.2,.16],linspace(0,1,255));
for k=1:3
    ax=nexttile(tiles);
    if k==1
        c=b(:,3); title_text='探头实测 B_z';
    elseif k==2
        c=p(:,3); title_text='留出层预测 B_z';
    else
        c=sqrt(sum((b-p).^2,2)); title_text='三分量残差模长';
    end
    scatter(ax,t.x_mm,t.y_mm,7,c,'filled'); axis(ax,'image'); axis(ax,'tight');
    xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); title(ax,title_text,'FontWeight','normal');
    cb=colorbar(ax); cb.Label.String='µT/A';
    if k<3
        clim(ax,[-scale scale]); colormap(ax,cmap);
    else
        clim(ax,[0 max(c)]); colormap(ax,parula(128));
    end
end
exportgraphics(fig,fullfile(here,'mc_b70_holdout.png'),'Resolution',150); close(fig)

fig=figure('Visible','off','Color','w','Position',[50 50 940 350]);
set(fig,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei');
tiles=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
names=["Magstim_DCC","Magstim_D70"];
for k=1:2
    d=readmatrix(fullfile(cache,names(k)+".tsv"),'FileType','text','NumHeaderLines',1,'Delimiter',char(9));
    ax=nexttile(tiles); hold(ax,'on');
    plot(ax,1e6*d(:,1),d(:,2),'Color',[.20 .45 .60],'LineWidth',1.4);
    plot(ax,1e6*d(:,1),d(:,3),'--','Color',[.75 .32 .16],'LineWidth',1.4);
    xlabel(ax,'时间 (µs)'); ylabel(ax,'发布波形的相对幅值');
    title(ax,strrep(names(k),'_',' '),'FontWeight','normal');
    legend(ax,{'拾取线圈波形','作者拟合'},'Box','off'); grid(ax,'on');
end
exportgraphics(fig,fullfile(here,'pickup_waveforms.png'),'Resolution',150); close(fig)
end

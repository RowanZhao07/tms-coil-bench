function plot_petrov_comparison()
% Visualize fixed-pose BSM-819 maps from validate_petrov_figure8, in microtesla.
here = fileparts(mfilename('fullpath'));
d = load(fullfile(here,'figure8_slice_results.mat'));
assert(isequal(d.conditions(1:2),["v1.0","v1.1"]));
outdir = fullfile(here,'figures');
if ~isfolder(outdir), mkdir(outdir); end
old = 1e6*d.prediction_T(:,:,1);
fixed = 1e6*d.prediction_T(:,:,2);
mri = 1e6*d.measured_T;
author = 1e6*d.author_T;
fieldMaps = cat(3,old,fixed,mri);
titles = ["旧摆位  v1.0","更正摆位  v1.1","体模 MRI 实测"];
blue = [.11 .32 .57]; red = [.68 .18 .17]; white = [.99 .99 .98];
colors = interp1([0 .5 1],[blue;white;red],linspace(0,1,257));
fig = figure('Visible','off','Color','w','Position',[40 40 1120 370]);
font = 'Microsoft YaHei';
set(fig,'DefaultAxesFontName',font,'DefaultTextFontName',font,'DefaultAxesFontSize',12);
t = tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');
for n = 1:3
    ax = nexttile(t); data = fieldMaps(:,:,n);
    h = imagesc(ax,data,[-180,180]); h.AlphaData = double(d.mask);
    axis(ax,'image'); ax.Color = [.93 .94 .94]; ax.Layer = 'top';
    xlim(ax,[10 150]); ylim(ax,[15 155]);
    ax.XTick = [30 80 130]; ax.YTick = [30 80 130];
    xlabel(ax,'像素列'); if n==1, ylabel(ax,'像素行'); end
    title(ax,titles(n),'FontWeight','normal'); colormap(ax,colors);
end
c = colorbar(ax); c.Layout.Tile = 'east'; c.Label.String='B_z  (µT)';
exportgraphics(fig,fullfile(outdir,'registered_fields.png'),'Resolution',175);
close(fig)

% Display coordinates after rot90: columns run toward decreasing scanner x,
% rows toward decreasing scanner z. The slice is 2 mm thick in scanner y.
fig = figure('Visible','off','Color','w','Position',[40 40 840 530]);
set(fig,'DefaultAxesFontName',font,'DefaultTextFontName',font,'DefaultAxesFontSize',12);
ax=axes(fig); hold(ax,'on'); ax.Color='w';
xmm=85.033-((1:160)-.5); zmm=-85.401+160.5-(1:160);
[~,phantom]=contourf(ax,xmm,zmm,double(d.mask),[.5 .5],'LineColor',[.45 .49 .47], ...
    'FaceColor',[.89 .93 .90],'LineWidth',1.1);
h=gobjects(1,2);
for n=1:2
    w=d.registered_wires_mm{n}; sep=nan(size(w,1),1);
    left=reshape([w(:,1),w(:,4),sep].',[],1); up=reshape([w(:,3),w(:,6),sep].',[],1);
    if n==1, col=red; else, col=blue; end
    h(n)=plot(ax,left,up,'Color',col,'LineWidth',.8);
    p=d.poses(1:3,4,n);
    plot(ax,p(1),p(3),'o','Color',col,'MarkerFaceColor','w','LineWidth',1.3);
end
axis(ax,'equal'); xlim(ax,[-105 105]); ylim(ax,[-92 108]); grid(ax,'on');
xlabel(ax,'扫描仪 x  (mm)'); ylabel(ax,'扫描仪 z  (mm)');
title(ax,'九匝八字绕组与实测体模切片的正面投影','FontWeight','normal');
legend(ax,[h,phantom],{'旧摆位','更正摆位','体模有效区'},'Location','southoutside','Orientation','horizontal','Box','off');
exportgraphics(fig,fullfile(outdir,'pose_projection.png'),'Resolution',175);
close(fig)

% Report RMS in pre-specified strips at the top/bottom of the published mask.
rowBands=[31 40;136 145];
strips = table(rowBands(:,1),rowBands(:,2),'VariableNames',{'first_row','last_row'});
for n=1:size(rowBands,1)
    valid = d.mask;
    valid(1:rowBands(n,1)-1,:)=false;
    valid(rowBands(n,2)+1:end,:)=false;
    strips.pixels(n)=nnz(valid);
    strips.old_to_author_rms(n)=norm(old(valid))/norm(author(valid));
    strips.fixed_to_author_rms(n)=norm(fixed(valid))/norm(author(valid));
end
writetable(strips,fullfile(here,'row_band_comparison.csv'));

% Row 50 is specified in the authors' original metrics.m.
row = 50; col = find(d.mask(row,:));
data = table(col',mri(row,col)',old(row,col)',fixed(row,col)',author(row,col)', ...
    'VariableNames',{'column_one_based','MRI_uT','v1_0_uT','v1_1_uT','author_uT'});
writetable(data,fullfile(here,'profile_row50.csv'));
fig = figure('Visible','off','Color','w','Position',[40 40 940 525]);
set(fig,'DefaultAxesFontName',font,'DefaultTextFontName',font,'DefaultAxesFontSize',12);
t = tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');
ax = nexttile(t); hold(ax,'on');
plot(ax,col,old(row,col),'--','Color',red,'LineWidth',1.7);
plot(ax,col,fixed(row,col),'-','Color',blue,'LineWidth',2);
plot(ax,col,mri(row,col),'k.','MarkerSize',8);
legend(ax,{'旧摆位','更正摆位','MRI 实测'},'Location','northeast','Box','off','NumColumns',3);
ylabel(ax,'B_z  (µT)'); title(ax,'同一条剖线上的有符号磁场','FontWeight','normal');
grid(ax,'on'); xlim(ax,[min(col),max(col)]); ax.XTickLabel = [];
ax = nexttile(t); hold(ax,'on'); yline(ax,0,':','Color',[.55 .55 .55]);
plot(ax,col,fixed(row,col)-mri(row,col),'-','Color',blue,'LineWidth',1.8);
plot(ax,col,author(row,col)-mri(row,col),'--','Color',[.2 .2 .2],'LineWidth',1.4);
legend(ax,{'','本项目更正摆位','作者发布的预测'},'Location','northeast','Box','off','NumColumns',2);
title(ax,'修正后仍保留的逐像素差值','FontWeight','normal');
ylabel(ax,'预测 − 实测  (µT)'); xlabel(ax,'第 50 行的像素列（从 1 起计）');
grid(ax,'on'); xlim(ax,[min(col),max(col)]);
exportgraphics(fig,fullfile(outdir,'profile_row50.png'),'Resolution',175);
close(fig)
end

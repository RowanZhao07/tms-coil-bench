function make_validation_figure
%MAKE_VALIDATION_FIGURE Plot independently executed solver verification.
root=fileparts(fileparts(mfilename('fullpath')));
v=readtable(fullfile(root,'results','validation','circle_analytic_convergence.csv'));
g=readtable(fullfile(root,'results','validation','geometry_convergence.csv'));
m=readtable(fullfile(root,'results','magpylib','comparison.csv'));
f=readtable(fullfile(root,'fem','results','axis_comparison.csv'),'VariableNamingRule','preserve');
c=readtable(fullfile(root,'fem','results','convergence_summary.csv'),'VariableNamingRule','preserve');
fig=figure('Visible','off','Color','w','Position',[30 30 1200 850]);
t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t);loglog(ax,v.segments,v.max_axis_relative_error*100,'o-','LineWidth',1.7);
xlabel(ax,'Segments per circle');ylabel(ax,'Max axis error (%)');grid(ax,'on');title(ax,'a  Polygonal circle vs closed-form circle');
ax=nexttile(t);ids=unique(m.id,'stable');values=zeros(numel(ids),1);
for k=1:numel(ids),values(k)=max(m.relative_L2(strcmp(m.id,ids{k})));end
barh(ax,values,'FaceColor',[.1 .4 .6]);set(ax,'XScale','log');yticks(ax,1:numel(ids));yticklabels(ax,ids);
set(ax,'TickLabelInterpreter','none','FontSize',8);xlabel(ax,'Relative vector L2 error');title(ax,'b  Magpylib vs MATLAB: all saved points');grid(ax,'on');
ax=nexttile(t);sel=strcmp(f.('case'),'boundary_080');p=f(sel,:);
plot(ax,p.z_m*1000,p.Bz_finite_cross_T*1000,'-','LineWidth',1.8,'DisplayName','Finite-section integral');hold(ax,'on');
plot(ax,p.z_m*1000,p.Bz_FEM_T*1000,'o','MarkerSize',7,'DisplayName','GetDP FEM');
xlabel(ax,'Depth z (mm)');ylabel(ax,'B_z (mT)');title(ax,'c  FEM benchmark: R=50 mm, 1000 A, one turn');grid(ax,'on');legend(ax,'Location','northeast');
ax=nexttile(t);ii=[4 3 5];plot(ax,c.domain_radius_m(ii),c.max_abs_finite_error_pct(ii),'s-','LineWidth',1.7);
xlabel(ax,'Air-domain radius (m)');ylabel(ax,'Max finite-section error (%)');title(ax,'d  FEM outer-boundary convergence');grid(ax,'on');
title(t,'Verification separates source discretization, solver implementation and FEM truncation','FontSize',13);
exportgraphics(fig,fullfile(root,'results','figures','validation_summary.png'),'Resolution',170);
exportgraphics(fig,fullfile(root,'results','figures','validation_summary.pdf'),'ContentType','vector');close(fig);
fprintf('Validation figure exported. Max 192/384 sampled geometry change %.6g%%\n',max(g.relative192_vs384)*100);
multiPath=fullfile(root,'fem','results','multiturn_main_circular','main_model_comparison.csv');
if exist(multiPath,'file')
    q=readtable(multiPath); fig=figure('Visible','off','Color','w','Position',[30 30 1100 480]);
    t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(t);plot(ax,q.z_m*1000,q.Bz_MATLAB_T*1000,'-','LineWidth',1.7,'DisplayName','MATLAB filament');hold(ax,'on');
    plot(ax,q.z_m*1000,q.Bz_FEM_T*1000,'o','MarkerSize',8,'DisplayName','GetDP finite section');
    plot(ax,q.z_m*1000,q.Bz_finite_cross_T*1000,'x','MarkerSize',8,'DisplayName','Finite-section integral');
    xlabel(ax,'Depth z (mm)');ylabel(ax,'B_z (mT)');grid(ax,'on');legend(ax,'Location','northeast');title(ax,'a  The actual 18-turn circular model');
    ax=nexttile(t);plot(ax,q.z_m*1000,q.FEM_vs_MATLAB_pct,'o-','DisplayName','FEM vs MATLAB');hold(ax,'on');
    plot(ax,q.z_m*1000,q.FEM_vs_finite_cross_pct,'s-','DisplayName','FEM vs finite-section integral');
    yline(ax,0,':','HandleVisibility','off');xlabel(ax,'Depth z (mm)');ylabel(ax,'Signed relative difference (%)');
    grid(ax,'on');legend(ax,'Location','southeast');title(ax,'b  Solver and source-model differences');
    title(t,'Main circular coil | 18 turns, 1000 A per turn | wire radius 2 mm in FEM','FontSize',14);
    exportgraphics(fig,fullfile(root,'results','figures','main_circular_fem_comparison.png'),'Resolution',170);
    exportgraphics(fig,fullfile(root,'results','figures','main_circular_fem_comparison.pdf'),'ContentType','vector');close(fig);
end
end

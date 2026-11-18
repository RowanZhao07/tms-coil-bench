function run_sensitivity(rootDir)
%RUN_SENSITIVITY Reproduce three supplementary numerical experiments.
% Usage: addpath('matlab'); run_sensitivity
% Writes only results/sensitivity. Uses MATLAB CPU double precision.
% All fields are free-space B, not tissue E or stimulation efficacy.
if nargin<1 || isempty(rootDir), rootDir=fileparts(fileparts(mfilename('fullpath'))); end
out=fullfile(rootDir,'results','sensitivity');
if ~exist(out,'dir'),mkdir(out);end
figdir=fullfile(out,'figures'); if ~exist(figdir,'dir'),mkdir(figdir);end
diary(fullfile(out,'matlab_sensitivity.log')); clean=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('SENSITIVITY START %s\n',char(datetime('now')));
tstart=tic; I=1000; n=192; coils=build_coils(n,rootDir,false);

%% 1. Conflicting pediatric figure-eight dimensions, equal NI comparison.
fprintf('1/3 Pediatric size conflict, equal current and equal effective turns\n');
compact=getcoil(coils,'pediatric_figure8');
large=larger_pediatric_variant(n);
dep=(10:1:50).'; obs=[zeros(numel(dep),2),dep/1000];
bc=field_segments(compact,obs,I,'exact',false);
bl=field_segments(large,obs,I,'exact',false);
nc=vecnorm(bc,2,2); nl=vecnorm(bl,2,2);
child=table(dep,bc(:,1),bc(:,2),bc(:,3),nc,bl(:,1),bl(:,2),bl(:,3),nl,nl./nc, ...
    'VariableNames',{'depth_mm','compact_Bx_T','compact_By_T','compact_Bz_T','compact_Bnorm_T', ...
    'larger_Bx_T','larger_By_T','larger_Bz_T','larger_Bnorm_T','larger_over_compact'});
writetable(child,fullfile(out,'pediatric_size_equal_NI.csv'));
write_segments(fullfile(out,'pediatric_larger_variant_geometry.csv'),large);
child_params=table(["compact_p9";"larger_p3_p4_controlled_variant"],[70;92],[17;14.5],[30.5;43], ...
    [4;4],[4;4],[4.5;4.5],[4;4],[32;32],[I;I],[I*32;I*32], ...
    'VariableNames',{'case_id','center_distance_mm','min_centerline_radius_mm','max_centerline_radius_mm', ...
    'turns_per_lobe_per_layer','layers','layer_pitch_mm','assumed_wire_OD_mm', ...
    'total_absolute_effective_turns','current_A','total_absolute_ampere_turns'});
writetable(child_params,fullfile(out,'pediatric_parameters.csv'));
f=figure('Visible','off','Color','w','Position',[100 100 920 520]);
plot(dep,nc*1e3,'LineWidth',2,'Color',[.1 .35 .65]); hold on;
plot(dep,nl*1e3,'LineWidth',2,'Color',[.8 .32 .12]);
xlabel('Depth from common wire plane (mm)'); ylabel('|B| on z axis (mT)');
title('Pediatric figure-of-eight: size uncertainty at equal NI');
legend({'Compact: source p.9','Larger: p.3-4 dimensions, controlled winding count'},'Location','northeast');
grid on; set(gca,'FontSize',11); exportplot(f,figdir,'pediatric_size_equal_NI');

%% 2. Same circular winding, only outer fold angle changes.
fprintf('2/3 Fold angle scan: 0,15,30,45,60 degrees\n');
base=getcoil(coils,'circular'); angles=[0 15 30 45 60];
xmm=linspace(-180,180,81); [X,Y]=meshgrid(xmm/1000,xmm/1000);
p=[X(:),Y(:),.02*ones(numel(X),1)];
fold=table(); fold_fields=cell(numel(angles),1);
for k=1:numel(angles)
    ang=angles(k); s=fold_outer_wings(base,.020,ang);
    b30=field_segments(s,[0 0 .03],I,'exact',false);
    bxy=field_segments(s,p,I,'exact',false);
    mag=reshape(vecnorm(bxy,2,2),size(X)); fold_fields{k}=mag;
    [area,clip]=half_area(mag,xmm);
    fold=[fold; table(ang,norm(b30)*1e3,max(mag(:))*1e3,area,clip, ...
        'VariableNames',{'fold_angle_deg','axis_B30_mT','xy20_peak_mT','xy20_half_area_cm2','half_area_clipped'})]; %#ok<AGROW>
    writetable(array2table([p,bxy,vecnorm(bxy,2,2)],'VariableNames', ...
        {'x_m','y_m','z_m','Bx_T','By_T','Bz_T','Bnorm_T'}), ...
        fullfile(out,sprintf('fold_angle_%02d_xy20.csv',ang)));
    fprintf('  angle%2d: B30 %.6f mT; A_half %.6f cm2\n',ang,norm(b30)*1e3,area);
end
writetable(fold,fullfile(out,'fold_angle_scan.csv'));
% Independent spatial-grid refinement for the actual main angle (30 deg).
xfine=linspace(-180,180,161); [XF,YF]=meshgrid(xfine/1000,xfine/1000);
pf=[XF(:),YF(:),.02*ones(numel(XF),1)];
bf=field_segments(fold_outer_wings(base,.020,30),pf,I,'exact',false);
mf=reshape(vecnorm(bf,2,2),size(XF)); [af,cf]=half_area(mf,xfine);
coarse=fold(fold.fold_angle_deg==30,:);
gridcheck=table([81;161],[4.5;2.25],[coarse.xy20_half_area_cm2;af], ...
    [coarse.xy20_peak_mT;max(mf(:))*1e3],[coarse.half_area_clipped;cf], ...
    'VariableNames',{'points_per_axis','grid_spacing_mm','xy20_half_area_cm2','xy20_peak_mT','half_area_clipped'});
writetable(gridcheck,fullfile(out,'fold30_area_grid_refinement.csv'));
f=figure('Visible','off','Color','w','Position',[100 100 1080 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(angles,fold.axis_B30_mT,'o-','LineWidth',2,'Color',[.1 .35 .65]);
xlabel('Outer fold angle (deg)'); ylabel('|B(0,0,30 mm)| (mT)'); title('Axial field'); grid on;set(gca,'FontSize',11);
nexttile; plot(angles,fold.xy20_half_area_cm2,'s-','LineWidth',2,'Color',[.8 .32 .12]);
xlabel('Outer fold angle (deg)'); ylabel('Half-peak area at z=20 mm (cm^2)'); title('Magnetic-field spatial extent'); grid on;set(gca,'FontSize',11);
sgtitle('Same 18-turn circular winding; hinge x=+/-20 mm; I=1000 A');
exportplot(f,figdir,'fold_angle_scan');

%% 3. Midpoint quadrature on long straight segments versus exact integral.
% Subdivision leaves every straight wire's geometry exactly unchanged.
fprintf('3/3 Long-segment midpoint quadrature refinement\n');
ids={'elliptical','h1_simplified','h1_alternative'};
names={'Racetrack','Simplified H1','Alternative H1'};
subdiv=[1 2 4 8 16 32];
testp=[0 0 .01;0 0 .02;0 0 .03;0 0 .04;0 0 .05;-.03 .02 .02;.03 -.02 .05];
writetable(array2table(testp,'VariableNames',{'x_m','y_m','z_m'}),fullfile(out,'quadrature_observers.csv'));
qt=table(); pointwise=table(); refcheck=table();
for c=1:numel(ids)
    s=getcoil(coils,ids{c}); exact=field_segments(s,testp,I,'exact',false);
    den=vecnorm(exact,2,2);
    for div=subdiv
        ss=split_segments(s,div);
        mid=field_segments(ss,testp,I,'midpoint',false);
        err=mid-exact;
        rel=norm(err,'fro')/norm(exact,'fro');
        pe=vecnorm(err,2,2)./max(den,1e-12);
        maxseg=max(vecnorm(ss(:,4:6)-ss(:,1:3),2,2))*1000;
        qt=[qt;table(string(ids{c}),div,size(ss,1),maxseg,rel,max(pe),max(vecnorm(err,2,2))*1e6, ...
            'VariableNames',{'id','subdivisions_per_original_segment','segments','max_segment_length_mm', ...
            'midpoint_relative_L2','max_pointwise_vector_relative_error','max_abs_vector_error_uT'})]; %#ok<AGROW>
        pointwise=[pointwise;table(repmat(string(ids{c}),7,1),repmat(div,7,1),(1:7)', ...
            exact(:,1),exact(:,2),exact(:,3),mid(:,1),mid(:,2),mid(:,3),pe, ...
            'VariableNames',{'id','subdivisions','observer_index','exact_Bx_T','exact_By_T','exact_Bz_T', ...
            'midpoint_Bx_T','midpoint_By_T','midpoint_Bz_T','vector_relative_error'})]; %#ok<AGROW>
    end
    % Verify that dividing the exact-integrated segments changes only rounding.
    exact32=field_segments(split_segments(s,32),testp,I,'exact',false);
    inv_err=norm(exact32-exact,'fro')/norm(exact,'fro');
    assert(inv_err<1e-10,'Straight-segment subdivision changed exact B.');
    refcheck=[refcheck;table(string(ids{c}),inv_err,'VariableNames',{'id','exact_subdivision32_invariance_relative_L2'})]; %#ok<AGROW>
    row1=qt(qt.id==ids{c}&qt.subdivisions_per_original_segment==1,:);
    row32=qt(qt.id==ids{c}&qt.subdivisions_per_original_segment==32,:);
    fprintf('  %s midpoint error: %.6g -> %.6g; exact invariance %.3g\n',ids{c},row1.midpoint_relative_L2,row32.midpoint_relative_L2,inv_err);
end
writetable(qt,fullfile(out,'midpoint_subdivision_convergence.csv'));
writetable(pointwise,fullfile(out,'midpoint_pointwise_fields.csv'));
writetable(refcheck,fullfile(out,'exact_subdivision_invariance.csv'));
f=figure('Visible','off','Color','w','Position',[100 100 920 520]);
clr=[.1 .35 .65;.8 .32 .12;.25 .55 .35];
for c=1:numel(ids)
    rows=qt(qt.id==ids{c},:);
    loglog(rows.subdivisions_per_original_segment,rows.midpoint_relative_L2,'o-','LineWidth',2,'Color',clr(c,:));hold on;
end
xlabel('Equal subdivisions per original straight segment'); ylabel('Relative vector L2 error versus exact segment integral');
title('Midpoint integration converges as long straight wires are subdivided');
legend(names,'Location','southwest');grid on;set(gca,'FontSize',11,'XTick',subdiv);
exportplot(f,figdir,'midpoint_subdivision_convergence');

%% Compact, experiment-oriented provenance and findings.
grid_rel=abs(af-coarse.xy20_half_area_cm2)/af;
meta=struct('matlab_version',version,'current_A',I,'segments_per_turn',n,'use_gpu',false, ...
    'source_document','Combination magnetic coils drawings, pages 3-4, 9, 16', ...
    'elapsed_s',toc(tstart),'completed_local',char(datetime('now')), ...
    'area_grid_relative_change',grid_rel, ...
    'quadrature_observers_count',size(testp,1), ...
    'interpretation','free-space magnetic B field; no tissue E, dose, or clinical inference');
fid=fopen(fullfile(out,'metadata.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true));fclose(fid);
fid=fopen(fullfile(out,'README_敏感性实验.md'),'w','n','UTF-8');
fprintf(fid,'# 补充参数敏感性实验\n\n运行：MATLAB中添加`matlab`目录后执行`run_sensitivity`。本目录所有CSV和图均来自实际MATLAB运行，I=1000 A，每匝192段，CPU双精度。\n\n');
fprintf(fid,'## 1. 儿童8字同名尺寸冲突\n\n原文p.9采用4层、4圈/层、中心距70 mm；p.3–4同名图则给出中心距92±2 mm、外径90±2 mm、内径25±1 mm。较大图未明确匝数/层数，故本实验保持两者每翼每层4圈、4层、总绝对32匝以及相同1000 A电流，只改变半径和中心距离；这是等NI的受控尺寸变体，不声称重建了较大版本的真实绕组。较大变体取中心线半径14.5…43 mm、假设导线OD4 mm、层距4.5 mm，满足内孔25、外径90和总宽182<185 mm。\n\n');
for d=[10 20 30 40 50]
    ix=find(dep==d); fprintf(fid,'- z=%d mm：较大/紧凑版本 |B| 比值 %.6f。\n',d,nl(ix)/nc(ix));
end
fprintf(fid,'\n这些比值同时受孔径、半径和中心距变化影响，不是单独半径的因果效应，也不能解释为儿童治疗效果。精确三分量和1 mm间隔深度数据见`pediatric_size_equal_NI.csv`。\n\n');
fprintf(fid,'## 2. 折叠角扫描\n\n保持与主圆形模型相同的18个闭合环、半径17…55 mm、2层、层距4.5 mm和相同电流，只将x=±20 mm以外部分向远离观测区方向折0/15/30/45/60°。这组角度与20 mm铰链距离是显式实验参数，原文p.16没有给出数值。该实验只回答所选几何在空气中的B场变化，不声称存在通用最优角。\n\n');
fprintf(fid,'| 折角(°) | B(0,0,30mm) (mT) | z20半值面积(cm²) |\n|---:|---:|---:|\n');
for k=1:height(fold), fprintf(fid,'| %d | %.6f | %.6f |\n',fold.fold_angle_deg(k),fold.axis_B30_mT(k),fold.xy20_half_area_cm2(k)); end
fprintf(fid,'\n半值面积定义为z=20 mm平面内 |B| ≥该平面最大值/2 的面积，正方形窗口x/y=±180 mm，主扫描81×81（4.5 mm）网格。30°加密到161×161（2.25 mm）后面积为%.6f cm²，相对变化%.3f%%。边界截断标志均写入CSV；二值阈值面积存在栅格误差，不宜把微小差异视为稳定结论。该面积不是脑组织电场聚焦面积。\n\n',af,100*grid_rel);
fprintf(fid,'## 3. 长直线的中点积分误差\n\n将每一条原有直线元等分1/2/4/8/16/32段；线元端点和直线路径不变，没有改变物理几何。对跑道形和两个示意H1，使用5个轴点(z=10…50 mm)加两个离轴点[-30,20,20]、[30,-20,50]mm比较中点积分与主程序的直线段解析积分。\n\n');
for c=1:numel(ids)
    r1=qt(qt.id==ids{c}&qt.subdivisions_per_original_segment==1,:);
    r32=qt(qt.id==ids{c}&qt.subdivisions_per_original_segment==32,:);
    fprintf(fid,'- %s：中点相对向量L2误差 %.4f%% → %.6f%%（1段→32段）。\n',ids{c},100*r1.midpoint_relative_L2,100*r32.midpoint_relative_L2);
end
fprintf(fid,'\n大误差来自把长直线只在一个中点采样。主实验从始至终采用每条直线的解析积分，不存在这类长直线中点求积误差。对同一路径将解析积分细分32次，结果相对变化均<1e-10，具体值见`exact_subdivision_invariance.csv`。这里仍不消除圆弧折线逼近、源几何尺寸、引线省略或H1源重建不确定性；这些与求积误差是不同问题。\n');
fclose(fid);
save(fullfile(out,'sensitivity.mat'),'child','child_params','fold','gridcheck','qt','pointwise','refcheck','fold_fields','xmm','meta');
fprintf('SENSITIVITY COMPLETE %.3f seconds\n',toc(tstart));
end

function s=getcoil(coils,id)
ix=find(strcmp({coils.id},id),1); assert(~isempty(ix));s=coils(ix).segments;
end

function s=larger_pediatric_variant(n)
t=linspace(0,2*pi,n+1).';s=zeros(0,7);
for side=[-1 1]
    for z=-(0:3)*4.5
        for r=linspace(14.5,43,4)
            p=[side*46+r*cos(t),r*sin(t),z+zeros(size(t))]*1e-3;p(end,:)=p(1,:);
            s=[s;p(1:end-1,:),p(2:end,:),side*ones(n,1)]; %#ok<AGROW>
        end
    end
end
end

function s=fold_outer_wings(s,h,deg)
for cols={1:3,4:6}
    q=cols{1};p=s(:,q);x=p(:,1);d=max(abs(x)-h,0);a=deg*pi/180;
    p(:,1)=sign(x).*(min(abs(x),h)+d*cos(a));p(:,3)=p(:,3)-d*sin(a);s(:,q)=p;
end
zmax=max([s(:,3);s(:,6)]);s(:,[3 6])=s(:,[3 6])-zmax;
end

function ss=split_segments(s,n)
a=s(:,1:3);d=s(:,4:6)-a;ss=zeros(size(s,1)*n,7);
for j=0:n-1
    ix=j*size(s,1)+(1:size(s,1));ss(ix,:)=[a+d*(j/n),a+d*((j+1)/n),s(:,7)];
end
end

function [a,clipped]=half_area(mag,xmm)
mask=mag>=.5*max(mag(:));a=trapz(xmm,trapz(xmm,double(mask),2))/100;
clipped=any(mask(1,:))||any(mask(end,:))||any(mask(:,1))||any(mask(:,end));
end

function write_segments(path,s)
writetable(array2table(s,'VariableNames',{'ax_m','ay_m','az_m','bx_m','by_m','bz_m','weight'}),path);
end

function exportplot(f,folder,name)
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',180);
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');close(f);
end

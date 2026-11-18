function run_study(useGPU)
%RUN_STUDY Reproduce all fourteen magnetic-field experiments in MATLAB.
% Run: addpath('matlab'); run_study(false)
% Air-core magnetostatic snapshot; current=1000 A. No tissue E-field claim.
if nargin<1, useGPU=false; end
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'results'); if ~exist(out,'dir'), mkdir(out); end
for d={'geometry','fields','figures','validation'}, if ~exist(fullfile(out,d{1}),'dir'), mkdir(fullfile(out,d{1})); end, end
diary(fullfile(out,'matlab_run.log')); cleanup=onCleanup(@()diary('off'));
fprintf('START %s\nMATLAB %s\n',char(datetime('now')),version);
timer=tic; I=1000; n=192;
coils=build_coils(n,root,true);
assert(numel(coils)==14,'Expected fourteen unique coil types.');
depth_mm=(10:1:100)'; axisPts=[zeros(numel(depth_mm),2),depth_mm/1000];
x_mm=linspace(-180,180,81); y_mm=x_mm;
[X,Y]=meshgrid(x_mm/1000,y_mm/1000);
xyPts=[X(:),Y(:),0.02*ones(numel(X),1)];
z_mm=linspace(10,100,46); [XS,ZS]=meshgrid(x_mm/1000,z_mm/1000);
xzPts=[XS(:),zeros(numel(XS),1),ZS(:)];
metrics=table(); depthAll=table(); conv=table();
allB=cell(14,1); fieldXY=cell(14,1); fieldXZ=cell(14,1);
for k=1:numel(coils)
    c=coils(k); tk=tic;
    fprintf('[%02d/14] %s: %d segments\n',k,c.id,size(c.segments,1));
    s=c.segments;
    writetable(array2table(s,'VariableNames',{'ax_m','ay_m','az_m','bx_m','by_m','bz_m','weight'}),fullfile(out,'geometry',[c.id '.csv']));
    Baxis=field_segments(s,axisPts,I,'exact',useGPU);
    Bxy=field_segments(s,xyPts,I,'exact',useGPU);
    Bxz=field_segments(s,xzPts,I,'exact',useGPU);
    allB{k}=Baxis; fieldXY{k}=reshape(vecnorm(Bxy,2,2),size(X)); fieldXZ{k}=reshape(vecnorm(Bxz,2,2),size(XS));
    write_field(fullfile(out,'fields',[c.id '_axis.csv']),axisPts,Baxis);
    write_field(fullfile(out,'fields',[c.id '_xy20.csv']),xyPts,Bxy);
    write_field(fullfile(out,'fields',[c.id '_xz.csv']),xzPts,Bxz);
    bn=vecnorm(Baxis,2,2); xy=fieldXY{k}; peak=max(xy(:));
    mask=xy>=0.5*peak; halfArea=trapz(y_mm,trapz(x_mm,double(mask),2))/100;
    clipped=any(mask(1,:))||any(mask(end,:))||any(mask(:,1))||any(mask(:,end));
    halfDepth=NaN; isMonotonic=all(diff(bn)<=max(bn)*1e-9);
    j=find(bn<=bn(1)/2,1);
    if isMonotonic && ~isempty(j) && j>1
        halfDepth=interp1(bn(j-1:j),depth_mm(j-1:j),bn(1)/2);
    end
    length_m=sum(vecnorm(s(:,4:6)-s(:,1:3),2,2).*abs(s(:,7)));
    % Equal wire-length weighted current: common 10 m reference.
    scale10m=10/length_m;
    vals=bn(ismember(depth_mm,[10,20,30,40,50]));
    metrics=[metrics; table(string(c.id),string(c.name_cn),size(s,1),length_m,vals(1)*1e3,vals(2)*1e3,vals(3)*1e3,vals(4)*1e3,vals(5)*1e3,peak*1e3,halfArea,halfDepth,isMonotonic,clipped,scale10m,vals(3)*1e3*scale10m,toc(tk),'VariableNames',{'id','name_cn','segments','weighted_length_m','B10_mT','B20_mT','B30_mT','B40_mT','B50_mT','xy20_peak_mT','xy20_half_area_cm2','axis_half_depth_mm','axis_monotonic','half_area_clipped','equal_10m_scale','B30_equal_10m_mT','runtime_s'})]; %#ok<AGROW>
    depthAll=[depthAll; table(repmat(string(c.id),numel(depth_mm),1),depth_mm,Baxis(:,1),Baxis(:,2),Baxis(:,3),bn,'VariableNames',{'id','depth_mm','Bx_T','By_T','Bz_T','Bnorm_T'})]; %#ok<AGROW>
    % Midpoint versus exact segments: quadrature only, not geometry uncertainty.
    pts=[axisPts(1:10:41,:); -0.03 0.02 0.02; 0.03 -0.02 0.05];
    ex=field_segments(s,pts,I,'exact',false); mid=field_segments(s,pts,I,'midpoint',false);
    err=norm(mid-ex,'fro')/norm(ex,'fro');
    conv=[conv; table(string(c.id),err,'VariableNames',{'id','midpoint_relative_L2'})]; %#ok<AGROW>
    save(fullfile(out,'fields',[c.id '.mat']),'c','I','axisPts','Baxis','xyPts','Bxy','xzPts','Bxz');
end
writetable(metrics,fullfile(out,'metrics.csv'));
writetable(depthAll,fullfile(out,'depth_profiles.csv'));
writetable(conv,fullfile(out,'validation','quadrature_comparison.csv'));
save(fullfile(out,'study.mat'),'coils','metrics','allB','fieldXY','fieldXZ','x_mm','y_mm','z_mm','depth_mm','I');
run_validation(coils,root,useGPU);
make_figures(root);
meta=struct('matlab_version',version,'current_A',I,'mu0_H_per_m',4*pi*1e-7,'segments_per_turn',n,'use_gpu',useGPU,'elapsed_s',toc(timer),'completed_local',char(datetime('now')),'coordinate_convention','all wire endpoints z<=0; observation depth z>0 from wire support plane','medium','homogeneous free space; no tissue, eddy current or ferromagnetic material');
fid=fopen(fullfile(out,'run_metadata.json'),'w'); fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true)); fclose(fid);
fprintf('COMPLETE %.3f seconds\n',toc(timer));
end

function write_field(path,p,B)
writetable(array2table([p B vecnorm(B,2,2)],'VariableNames',{'x_m','y_m','z_m','Bx_T','By_T','Bz_T','Bnorm_T'}),path);
end

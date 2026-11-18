function run_spatial_convergence(useGPU)
%RUN_SPATIAL_CONVERGENCE Recompute half-value area on a finer observation grid.
% This tests sampling error, independently of wire discretization error.
if nargin<1, useGPU=false; end
root=fileparts(fileparts(mfilename('fullpath'))); d=load(fullfile(root,'results','study.mat'));
out=fullfile(root,'results','spatial_convergence');if ~exist(out,'dir'),mkdir(out);end
xy_mm=linspace(-180,180,161);[X,Y]=meshgrid(xy_mm/1000);
p=[X(:),Y(:),.02*ones(numel(X),1)];t=table();
for k=1:14
    clock=tic;B=field_segments(d.coils(k).segments,p,d.I,'exact',useGPU);
    Bmag=reshape(vecnorm(B,2,2),size(X));peak=max(Bmag(:));mask=Bmag>=.5*peak;
    area=trapz(xy_mm,trapz(xy_mm,double(mask),2))/100;
    clipping=any(mask(1,:))||any(mask(end,:))||any(mask(:,1))||any(mask(:,end));
    a81=d.metrics.xy20_half_area_cm2(k);pk81=d.metrics.xy20_peak_mT(k);
    t=[t;table(string(d.coils(k).id),a81,area,(a81-area)/area,pk81,peak*1e3,clipping,toc(clock),...
        'VariableNames',{'id','area81_cm2','area161_cm2','area_relative_change','peak81_mT','peak161_mT','clipped161','runtime_s'})]; %#ok<AGROW>
    save(fullfile(out,[d.coils(k).id '.mat']),'xy_mm','Bmag','area','peak');
    fprintf('SPATIAL %s area81=%.4f area161=%.4f cm2\n',d.coils(k).id,a81,area);
end
writetable(t,fullfile(out,'plane_grid_convergence.csv'));
meta=struct('use_gpu',useGPU,'coarse_grid',81,'fine_grid',161,'fine_spacing_mm',2.25,'plane_depth_mm',20,'domain_mm',[-180 180],'threshold','half of the peak on each respective grid','note','Two observation grids estimate sampling sensitivity; not a rigorous continuum error bound. Main metrics.csv remains the recorded81x81 baseline.');
fid=fopen(fullfile(out,'metadata.json'),'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);
end

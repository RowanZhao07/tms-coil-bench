function scores = run_pose_refinement()
% Estimate scanner-z translation from sparse physical Bz measurements.
here=fileparts(mfilename('fullpath'));
root=fileparts(fileparts(here));
source=fullfile(root,'placement_pilot','public_benchmark','matched_figure8');
addpath(fullfile(root,'matlab'),source)
input=fullfile(source,'figure8_slice_results.mat');
if ~isfile(input), validate_petrov_figure8; end
data=load(input,'mask','measured_T','registered_wires_mm');
mask=logical(data.mask);
measured=1e6*data.measured_T;
% v1.0 is the only navigation file used in the inverse calculation.
segments=signed_segments(data.registered_wires_mm{1});

[r,c]=ndgrid(45:5:85,46:8:118);
pool=sub2ind(size(mask),r(:),c(:));
pool=pool(mask(pool));
[r,c]=ndgrid([55 65 75],[60 105]);
fixed=sub2ind(size(mask),r(:),c(:));
[row,col]=ind2sub(size(mask),pool);
% Choose six positions without looking at any MRI field values.
designOffsets=[-20 0 20 40];
fisher=zeros(numel(pool),2,numel(designOffsets));
poolNodes=pixel_nodes(size(mask),pool);
for k=1:numel(designOffsets)
    z=designOffsets(k);
    atPose=predict(segments,poolNodes,z);
    slope=(predict(segments,poolNodes,z+0.5)-predict(segments,poolNodes,z-0.5));
    fisher(:,:,k)=[slope/norm(slope),atPose/norm(atPose)];
end
pick=choose_sensors(row,col,fisher,6);
designed=pool(pick);
designs={"fixed6",fixed; "physics6",designed; "physics3",designed(1:3); "physics2",designed(1:2)};
hold=mask;
hold(1:89,:)=false;
assert(nnz(hold)==5325 && all(mask(fixed)) && ~any(hold([fixed;pool])));
siteTable=table();
parameters=table();
for i=1:size(designs,1)
    name=designs{i,1}; indices=designs{i,2};
    [pr,pc]=ind2sub(size(mask),indices);
    siteTable=[siteTable;table(repmat(name,numel(indices),1),(1:numel(indices))',pr,pc, ...
        measured(indices),'VariableNames',{'design','order','row','column','measured_uT'})]; %#ok<AGROW>
    nodes=pixel_nodes(size(mask),indices);
    yy=measured(indices);
    for model=["unchanged","gain","translation","translation_gain"]
        [z,g,objective]=fit_model(segments,nodes,yy,model);
        pp=predict(segments,nodes,z);
        dp=predict(segments,nodes,z+0.5)-predict(segments,nodes,z-0.5);
        similarity=abs(dp'*pp)/(norm(dp)*norm(pp));
        pars=table(name,model,numel(indices),z,g,sqrt(objective/numel(yy)), ...
            100*sqrt(objective)/norm(yy),similarity,'VariableNames', ...
            {'design','model','n_train','dz_mm','gain','train_RMSE_uT','train_L2_percent','derivative_cosine'});
        parameters=[parameters;pars]; %#ok<AGROW>
    end
end
writetable(siteTable,fullfile(here,'calibration_points.csv'));
% Random-location controls use the identical candidate pool and spacing rule.
% All designs use the same 0.5 mm grid of directly calculated fields.
[randomControls,locationContext]=random_layouts(segments,poolNodes,pool,fixed,designed, ...
    row,col,measured,mask);
% Open the test values and the v1.1 reference only after all fits are complete.
ref=load(input,'poses','prediction_T');
knownDz=ref.poses(3,4,2)-ref.poses(3,4,1);
randomControls.abs_dz_difference_from_author_mm=abs(randomControls.dz_mm-knownDz);
writetable(randomControls,fullfile(here,'random_layouts.csv'));
writetable(locationContext,fullfile(here,'random_layout_sites.csv'));
testNodes=pixel_nodes(size(mask),find(hold));
truth=measured(hold);
oldMap=1e6*ref.prediction_T(:,:,1);
correctMap=1e6*ref.prediction_T(:,:,2);
assert(max(abs(predict(segments,pixel_nodes(size(mask),fixed),0)-oldMap(fixed)))<1e-8);
scores=parameters;
scores.test_L2_percent=zeros(height(scores),1);
scores.test_RMSE_uT=zeros(height(scores),1);
% Cache repeated estimates, including the common unshifted geometry.
[offsets,~,idx]=unique(scores.dz_mm);
fields=zeros(numel(truth),numel(offsets));
for k=1:numel(offsets), fields(:,k)=predict(segments,testNodes,offsets(k)); end
for i=1:height(scores)
    residual=scores.gain(i)*fields(:,idx(i))-truth;
    scores.test_L2_percent(i)=100*norm(residual)/norm(truth);
    scores.test_RMSE_uT(i)=sqrt(mean(residual.^2));
end
scores.abs_dz_difference_from_author_mm=abs(scores.dz_mm-knownDz);
refLabels=["old v1.0";"corrected v1.1";"old plus author z only"];
v=[oldMap(hold),correctMap(hold),predict(segments,testNodes,knownDz)];
assert(max(abs(v(:,1)-fields(:,offsets==0)))<1e-8);
refScores=table(refLabels,[0;knownDz;knownDz],100*(vecnorm(v-truth)./norm(truth))', ...
    sqrt(mean((v-truth).^2))','VariableNames',{'reference','dz_mm','test_L2_percent','test_RMSE_uT'});
writetable(scores,fullfile(here,'pose_ablation.csv'));
writetable(refScores,fullfile(here,'reference_only.csv'));
negative_control(here,signed_segments(data.registered_wires_mm{2}),measured,mask, ...
    fixed,designed,testNodes,truth);
figure_summary(here,mask,measured,fixed,designed,scores,refScores);
disp(scores);disp(refScores)
end

function s=signed_segments(w)
wingBreak=find(vecnorm(w(1:end-1,4:6)-w(2:end,1:3),2,2)>1e-5);
assert(isscalar(wingBreak) && size(w,1)==2700);
s=[w*1e-3,[ones(wingBreak,1);-ones(size(w,1)-wingBreak,1)]];
end

function negative_control(here,segments,y,mask,fixed,designed,testNodes,testTruth)
% Reuse the same sensors with the already-corrected navigation file.
out=table();
for name=["fixed6","physics6"]
    indices=fixed;
    if name=="physics6", indices=designed; end
    nodes=pixel_nodes(size(mask),indices);
    for model=["unchanged","gain","translation","translation_gain"]
        [dz,g,loss]=fit_model(segments,nodes,y(indices),model);
        err=g*predict(segments,testNodes,dz)-testTruth;
        out=[out;table(name,model,dz,g,100*sqrt(loss)/norm(y(indices)), ...
            100*norm(err)/norm(testTruth),sqrt(mean(err.^2)), ...
            'VariableNames',{'design','model','additional_dz_mm','gain','train_L2_percent', ...
            'test_L2_percent','test_RMSE_uT'})]; %#ok<AGROW>
    end
end
writetable(out,fullfile(here,'corrected_start_control.csv'));
disp(out);
end

function [summary,locations]=random_layouts(s,nodes,pool,fixed,selected,row,col,y,mask)
grid=-40:0.5:60;
dict=zeros(numel(pool),numel(grid));fixedDict=zeros(numel(fixed),numel(grid));
fixedNodes=pixel_nodes(size(mask),fixed);
for k=1:numel(grid)
    dict(:,k)=predict(s,nodes,grid(k));fixedDict(:,k)=predict(s,fixedNodes,grid(k));
end
stream=RandStream('mt19937ar','Seed',1806);
summary=table();locations=table();
for l=1:102
    if l==1
        name="fixed6";rows=fixed;pred=fixedDict;
    elseif l==2
        name="physics6";rows=selected;[~,ix]=ismember(selected,pool);pred=dict(ix,:);
    else
        name="random";active=true(numel(pool),1);indices=zeros(6,1);
        for j=1:6
            avail=find(active);i=avail(randi(stream,numel(avail)));indices(j)=i;
            active=active & hypot(row-row(i),col-col(i))>=8;
        end
        rows=pool(indices);pred=dict(indices,:);
    end
    target=y(rows);
    gain=max(0,(target'*pred)./sum(pred.^2,1));
    loss=sum((pred.*gain-target).^2,1);
    [~,best]=min(loss);
    summary=[summary;table(name,l,grid(best),gain(best), ...
        100*sqrt(loss(best))/norm(target),'VariableNames', ...
        {'design','layout_id','dz_mm','gain','train_L2_percent'})]; %#ok<AGROW>
    [rr,cc]=ind2sub(size(mask),rows);
    locations=[locations;table(repmat(l,6,1),rr,cc,'VariableNames', ...
        {'layout_id','row','column'})]; %#ok<AGROW>
end
end

function chosen=choose_sensors(row,col,F,n)
% Greedy robust D-optimal criterion over the four assumed translations.
chosen=zeros(n,1); avail=true(numel(row),1);
G=repmat(1e-6*eye(2),1,1,size(F,3));
for j=1:n
    score=-inf(size(row));
    for i=find(avail)'
        candidate=inf;
        for k=1:size(F,3)
            a=F(i,:,k);
            candidate=min(candidate,log(det(G(:,:,k)+a'*a)));
        end
        score(i)=candidate;
    end
    [~,best]=max(score);chosen(j)=best;
    for k=1:size(F,3)
        a=F(best,:,k);G(:,:,k)=G(:,:,k)+a'*a;
    end
    avail=avail & hypot(row-row(best),col-col(best))>=8;
end
end

function [shift,gain,objective]=fit_model(segments,nodes,y,model)
freeGain=any(model==["gain","translation_gain"]);
if any(model==["translation","translation_gain"])
    % Coarse global bracket followed by direct full-field objective refinement.
    grid=-40:5:60;
    vals=arrayfun(@(z) evaluate(z),grid);
    [~,j]=min(vals);
    lo=grid(max(1,j-1));hi=grid(min(numel(grid),j+1));
    [candidate,d]=fminbnd(@evaluate,lo,hi,optimset('TolX',0.005,'Display','off'));
    if vals(j)<=d, shift=grid(j);else,shift=candidate;end
else
    shift=0;
end
[objective,gain]=evaluate(shift);
    function [loss,g]=evaluate(z)
        p=predict(segments,nodes,z);
        g=1;
        if freeGain, g=max(0,(p'*y)/(p'*p)); end
        loss=sum((g*p-y).^2);
    end
end

function nodes=pixel_nodes(sz,ids)
% Inverse of rot90(pixel grid), then eight corners per voxel in scanner space.
[r,c]=ind2sub(sz,ids(:));
i=c-1;j=sz(1)-r;
[di,dj,dk]=ndgrid(0:1,0:1,0:1);
ni=i+di(:)';nj=j+dj(:)';nk=zeros(size(i))+dk(:)';
[linear,~,back]=unique(sub2ind([161,161,2],ni(:)+1,nj(:)+1,nk(:)+1));
[ii,jj,kk]=ind2sub([161,161,2],linear);
nodes.xyz=1e-3*[85.033-(ii-1),53.931+2*(kk-1),-85.401+(jj-1)];
nodes.corners=reshape(back,numel(ids),8);
end

function p=predict(segments,nodes,shift)
observers=nodes.xyz;
observers(:,3)=observers(:,3)-1e-3*shift;
b=field_segments(segments,observers,3.576,'exact',false);
z=b(:,3);
p=1e6*mean(z(nodes.corners),2);
end

function figure_summary(here,mask,measured,fixed,design,scores,reference)
fig=figure('Visible','off','Color','white','Position',[50 50 1040 430]);
t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t); im=imagesc(ax,measured); im.AlphaData=double(mask);
ax.Color=[.93 .93 .93];axis(ax,'image'); clim(ax,[-130 130]);colormap(ax,parula);
cb=colorbar(ax);cb.Label.String='B_z (µT)';hold(ax,'on');
[r,c]=ind2sub(size(mask),fixed);plot(ax,c,r,'ko','MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.2);
[r,c]=ind2sub(size(mask),design);plot(ax,c,r,'mx','MarkerSize',10,'LineWidth',1.8);
yline(ax,89.5,'k--'); xlabel(ax,'Image column');ylabel(ax,'Image row');
legend(ax,{'Fixed six','Sensitivity-selected six','Test boundary'},'Location','south','FontSize',8);
title(ax,'Calibration sites and untouched test region');
ax=nexttile(t);a=scores(scores.design=="fixed6",:);b=scores(scores.design=="physics6",:);
hold(ax,'on');
plot(ax,(1:4)-0.08,a.test_L2_percent,'o','MarkerFaceColor',[0 .37 .59], ...
    'MarkerEdgeColor',[0 .37 .59],'MarkerSize',8);
plot(ax,(1:4)+0.08,b.test_L2_percent,'s','MarkerFaceColor',[.85 .37 .03], ...
    'MarkerEdgeColor',[.85 .37 .03],'MarkerSize',8);
ax.YScale='log';ax.YLim=[18 320];ax.YTick=[20 25 40 80 150 250];
ax.XLim=[0.6 4.45];grid(ax,'on');
ax.XTick=1:4;ax.XTickLabel={'Frozen','Gain','Shift','Shift + gain'};
yline(ax,reference.test_L2_percent(2),'k--');
legend(ax,{'Fixed six','Sensitivity-selected six','Full author correction'},'Location','northeast');
ylabel(ax,'Held-out relative L2 error (%)');title(ax,'5,325 held-out MRI pixels');
exportgraphics(fig,fullfile(here,'calibration_ablation.png'),'Resolution',160);close(fig);
end

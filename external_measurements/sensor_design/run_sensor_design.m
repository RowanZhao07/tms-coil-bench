function run_sensor_design
% Sparse field measurements to recover an in-plane coordinate error and gain.
here = fileparts(mfilename('fullpath'));
parent = fileparts(here);
cache = fullfile(fileparts(parent),'tmp','measured_coil_library','mc_b70_scan.mat');
d = load(cache,'xyz_mm','B_T_per_A');
dipoles = readmatrix(fullfile(parent,'measured_coil_library','mc_b70_equivalent_sources.csv'));
near = coarse_scan(d,23.8);
validation = coarse_scan(d,33.7);
test = [coarse_scan(d,-34.5),coarse_scan(d,43.65)];

% Strength, derivatives and reliability are calculated before reading the near-plane
% or test-plane observations. Only validation readings enter reliability.
candidate = near.xyz(all(abs(near.xyz(:,1:2))<=56,2),:);
b = magnetic_field(candidate,dipoles);
deriv = field_derivatives(candidate,dipoles);
val_residual = magnetic_field(validation.xyz,dipoles)-validation.obs;
square_error = sum(val_residual.^2,2);
floor2 = median(square_error);
rel = zeros(size(candidate,1),1);
for i = 1:size(candidate,1)
    dxy2 = sum((validation.xyz(:,1:2)-candidate(i,1:2)).^2,2);
    [~,nine] = mink(dxy2,9);
    rel(i) = 1/(floor2+mean(square_error(nine)));
end
rel = rel/median(rel);
info = zeros(3,3,size(candidate,1));
for i = 1:size(candidate,1)
    J = [7*deriv(i,:,1)',7*deriv(i,:,2)',0.1*b(i,:)'];
    info(:,:,i) = J'*J;
end
unweighted = design(candidate,info,ones(size(rel)),6);
weighted = design(candidate,info,rel,12);
[~,strong_first] = sort(vecnorm(b,2,2),'descend');
peak = separated_points(candidate,strong_first,6);
fixed_xy = [-35 -28;35 28;-35 28;35 -28;-35 0;35 0];
[found,idx] = ismember(fixed_xy,candidate(:,1:2),'rows'); assert(all(found));
names = ["Fixed6","Peak6","Dopt6","Weighted6","Weighted3","Weighted12"];
layouts = {idx(:),peak,unweighted,weighted(1:6),weighted(1:3),weighted};
rng(20261118,'twister');
for k=1:24
    names(end+1)=compose("Random%02d",k); %#ok<AGROW>
    layouts{end+1}=separated_points(candidate,randperm(size(candidate,1)),6); %#ok<AGROW>
end
% Follow-up negative control: preserve the reliability histogram while breaking
% its spatial association. Use the same fitter, costs and six-probe budget.
rng(27,'twister');
for k=1:24
    names(end+1)=compose("Shuffled%02d",k); %#ok<AGROW>
    layouts{end+1}=design(candidate,info,rel(randperm(numel(rel))),6); %#ok<AGROW>
end
locations = cell(size(layouts));
for k=1:numel(layouts)
    sel=layouts{k}(:);
    locations{k}=table(repmat(names(k),numel(sel),1),(1:numel(sel))', ...
        candidate(sel,1),candidate(sel,2),rel(sel),vecnorm(b(sel,:),2,2), ...
        'VariableNames',{'layout','pick','x_mm','y_mm','relative_weight','model_norm_uT_A'});
end
writetable(vertcat(locations{:}),fullfile(here,'sensor_locations.csv'));
writetable(array2table([candidate,rel,vecnorm(b,2,2)],'VariableNames', ...
    {'x_mm','y_mm','z_mm','relative_weight','model_norm_uT_A'}),fullfile(here,'candidate_weights.csv'));

% All shifted references are actual neighbouring scan readings. Gain changes
% are explicitly synthetic. Test points are aligned by exact lattice lookup.
for k=1:2
    test(k).nominal=test(k).xyz(all(abs(test(k).xyz(:,1:2))<=126,2),:);
    assert(size(test(k).nominal,1)==37^2);
end
% Precompute the smooth test-plane model; off-grid exact checks bound the
% numerical interpolation error separately from experimental error.
check=zeros(2,1);
surfaces=cell(2,1);
for k=1:2
    [surfaces{k},check(k)] = test_interpolants(test(k).nominal,dipoles);
end
writetable(table([test.z]',check,'VariableNames',{'z_mm','offgrid_exact_check_relative_L2_pct'}), ...
    fullfile(here,'interpolation_check.csv'));
assert(all(check<0.05),'Test interpolation insufficiently resolved.');
[dx,dy]=ndgrid([-7 0 7],[-7 0 7]);
trials=cell(0,1);
for n=1:numel(dx)
    shift=[dx(n),dy(n)];
    true_gain=1+0.10*any(shift~=0);
    target=cell(2,1);
    for k=1:2
        target{k}=true_gain*lookup(test(k),test(k).nominal(:,1:2)+shift);
    end
    for j=1:numel(names)
        pos=candidate(layouts{j},:);
        y=true_gain*lookup(near,pos(:,1:2)+shift);
        modes="Joint";
        if names(j)=="Weighted6", modes=["None","Gain","Joint"]; end
        for mode=modes
            [fit,calib,hit]=fit_pose(pos,y,dipoles,mode);
            sse=zeros(2,1); energy=sse;
            for k=1:2
                xy=test(k).nominal(:,1:2)+fit(1:2);
                estimated=zeros(size(target{k}));
                for component=1:3
                    F=surfaces{k}{component};
                    estimated(:,component)=fit(3)*F(xy(:,1),xy(:,2));
                end
                residual=estimated-target{k};
                sse(k)=sum(residual.^2,'all'); energy(k)=sum(target{k}.^2,'all');
            end
            trials{end+1,1}=table(names(j),mode,size(pos,1),dx(n),dy(n),true_gain, ...
                fit(1),fit(2),fit(3),norm(fit(1:2)-shift),100*abs(fit(3)-true_gain), ...
                calib,100*sqrt(sum(sse)/sum(energy)),100*sqrt(sse(1)/energy(1)), ...
                100*sqrt(sse(2)/energy(2)),hit,'VariableNames',{'layout','fit','budget', ...
                'injected_dx_mm','injected_dy_mm','injected_gain','fit_dx_mm','fit_dy_mm', ...
                'fit_gain','position_error_mm','gain_error_pp','calibration_L2_pct', ...
                'test_L2_pct','test_negative_L2_pct','test_positive_L2_pct','at_bound'}); %#ok<AGROW>
        end
    end
end
t=vertcat(trials{:});
writetable(t,fullfile(here,'injection_trials.csv'));
[gr,lay,mode]=findgroups(t.layout,t.fit);
fault=t.injected_dx_mm~=0 | t.injected_dy_mm~=0;
measure=@(v) splitapply(@(a,f) mean(a(f)),v,fault,gr);
control=splitapply(@(a,f) a(~f),t.test_L2_pct,fault,gr);
summary=table(lay,mode,measure(t.position_error_mm),measure(t.gain_error_pp), ...
    measure(t.test_L2_pct),measure(t.test_negative_L2_pct),measure(t.test_positive_L2_pct),control, ...
    'VariableNames',{'layout','fit','mean_position_error_mm','mean_gain_error_pp','mean_test_L2_pct', ...
    'mean_negative_L2_pct','mean_positive_L2_pct','control_test_L2_pct'});
writetable(summary,fullfile(here,'layout_summary.csv'));
% Post-scoring diagnostic: test whether high-error development locations rank
% similarly on other planes. No measured error from these planes enters design.
planes=[near,rmfield(test,'nominal')];
drift=zeros(numel(planes),1);
for k=1:numel(planes)
    s=planes(k); e=sum((magnetic_field(s.xyz,dipoles)-s.obs).^2,2);
    local_mse=zeros(size(candidate,1),1);
    for i=1:numel(local_mse)
        [~,neighbours]=mink(sum((s.xyz(:,1:2)-candidate(i,1:2)).^2,2),9);
        local_mse(i)=mean(e(neighbours));
    end
    cor=corrcoef(average_rank(1./rel),average_rank(local_mse));
    drift(k)=cor(1,2);
end
writetable(table([planes.z]',drift,'VariableNames',{'target_z_mm','spatial_spearman_rho'}), ...
    fullfile(here,'residual_transfer.csv'));
disp(summary(~startsWith(summary.layout,'Random') & ~startsWith(summary.layout,'Shuffled'),:))
disp(median(summary.mean_test_L2_pct(startsWith(summary.layout,'Random'))))
end

function r=average_rank(v)
[~,~,group]=unique(v);
counts=accumarray(group,1); mid=cumsum(counts)-(counts-1)/2;
r=mid(group);
end

function s=coarse_scan(d,z)
% The two source raster scans have slightly different z, and only the coarse
% one has |x|>50. Keep one measurement per (x,y) on its complete 7 mm grid.
zmm=d.xyz_mm(:,3);
level=median(zmm(abs(zmm-z)<.2 & abs(d.xyz_mm(:,1))>50));
take=abs(zmm-level)<.001;
s.z=level; s.xyz=d.xyz_mm(take,:); s.obs=1e6*d.B_T_per_A(take,:);
assert(size(s.xyz,1)==41^2 && size(unique(s.xyz(:,1:2),'rows'),1)==41^2);
end

function y=lookup(s,xy)
[found,idx]=ismember(xy,s.xyz(:,1:2),'rows'); assert(all(found));
y=s.obs(idx,:);
end

function b=magnetic_field(p,src)
% SI dipole source positions (input mm), moments m^2, field output µT/A.
b=zeros(size(p));
for a=1:256:size(p,1)
    ix=a:min(a+255,size(p,1));
    rx=(p(ix,1)-src(:,1)')*.001; ry=(p(ix,2)-src(:,2)')*.001;
    rz=(p(ix,3)-src(:,3)')*.001;
    r2=rx.^2+ry.^2+rz.^2;
    h=.1*r2.^(-1.5);
    scale=3*(rx.*src(:,4)'+ry.*src(:,5)'+rz.*src(:,6)')./r2;
    b(ix,:)=[sum(h.*(rx.*scale-src(:,4)'),2),sum(h.*(ry.*scale-src(:,5)'),2), ...
        sum(h.*(rz.*scale-src(:,6)'),2)];
end
end

function j=field_derivatives(p,src)
h=.2; j=zeros(size(p,1),3,2);
for axis=1:2
    inc=[0 0 0]; inc(axis)=h;
    j(:,:,axis)=(magnetic_field(p+inc,src)-magnetic_field(p-inc,src))/(2*h);
end
end

function ids=design(candidate,info,w,n)
total=eye(3)*1e-6*max(squeeze(sum(sum(info.*eye(3),1),2)));
valid=true(size(w)); ids=zeros(n,1);
for k=1:n
    objective=-inf(size(w));
    for i=find(valid)'
        R=chol(total+w(i)*info(:,:,i));
        objective(i)=2*sum(log(diag(R)));
    end
    [~,pick]=max(objective); ids(k)=pick;
    total=total+w(pick)*info(:,:,pick);
    valid=valid & sum((candidate(:,1:2)-candidate(pick,1:2)).^2,2)>=14^2;
end
end

function ids=separated_points(x,order,n)
ids=zeros(n,1); count=0;
for i=order(:)'
    if all(sum((x(ids(1:count),1:2)-x(i,1:2)).^2,2)>=14^2)
        count=count+1; ids(count)=i;
        if count==n, return; end
    end
end
error('Insufficient separated candidates.')
end

function [fit,err,bound]=fit_pose(x,y,src,mode)
fit=[0 0 1]; lower=[-14 -14 .6]; upper=[14 14 1.4];
b=magnetic_field(x,src);
if mode=="Gain"
    fit(3)=max(lower(3),min(upper(3),sum(y.*b,'all')/sum(b.^2,'all')));
elseif mode=="Joint"
    for iter=1:24
        p=x+[fit(1:2),0];
        b=magnetic_field(p,src); der=field_derivatives(p,src);
        residual=fit(3)*b-y;
        J=[fit(3)*reshape(der(:,:,1),[],1),fit(3)*reshape(der(:,:,2),[],1),b(:)];
        step=-(J\residual(:))';
        prev=sum(residual.^2,'all');
        next=fit;
        for alpha=2.^-(0:12)
            proposal=max(lower,min(upper,fit+alpha*step));
            pred=proposal(3)*magnetic_field(x+[proposal(1:2),0],src);
            if sum((pred-y).^2,'all')<prev, next=proposal; break; end
        end
        done=norm(next-fit)<1e-6; fit=next;
        if done, break; end
    end
end
pred=fit(3)*magnetic_field(x+[fit(1:2),0],src);
err=100*norm(pred-y,'fro')/norm(y,'fro');
bound=any(abs(fit-lower)<1e-5 | abs(fit-upper)<1e-5);
end

function [fields,check]=test_interpolants(nominal,src)
axis=(-147:3.5:147)';
[xx,yy]=ndgrid(axis,axis);
surface=magnetic_field([xx(:),yy(:),repmat(nominal(1,3),numel(xx),1)],src);
fields=cell(3,1);
points=nominal(1:13:end,:)+[1.3 -2.3 0]; estimated=zeros(size(points));
for k=1:3
    fields{k}=griddedInterpolant({axis,axis},reshape(surface(:,k),size(xx)),'spline','none');
    estimated(:,k)=fields{k}(points(:,1),points(:,2));
end
actual=magnetic_field(points,src);
check=100*norm(actual-estimated,'fro')/norm(actual,'fro');
end

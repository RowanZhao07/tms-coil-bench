function scores = run_mc_b70
% Fit regularized equivalent dipoles to measured MC-B70 scans; hold out whole planes.
here=fileparts(mfilename('fullpath'));
cache=fullfile(here,'..','..','tmp','measured_coil_library');
scan=load(fullfile(cache,'mc_b70_scan.mat'));
xyz=scan.xyz_mm/1000;
B=scan.B_T_per_A;
assert(isequal(size(xyz),[95773 3]) && all(isfinite(B),'all'));
casing=stlread(fullfile(cache,'MagVenture_MC-B70_dipole-fit.stl'));
sites=source_sites(casing,xyz);
% Merge the slight offsets between the two grid spacings on each scanner layer.
z_mm=scan.xyz_mm(:,3);
test_z=[-34.5,23.8,43.65];
test=any(abs(z_mm-test_z)<0.2,2);
validation=abs(z_mm-33.7)<0.2;
train=find(~test & ~validation);
train=train(1:8:end); % thin only training to bound memory and runtime
idx_valid=find(validation);
npar=3*size(sites,1);
gram=zeros(npar); rhs=zeros(npar,1);
for k=1:256:numel(train)
    sel=train(k:min(k+255,end));
    G=dipole_matrix(xyz(sel,:),sites);
    target=B(sel,:); gram=gram+G'*G; rhs=rhs+G'*target(:);
end
[V,eigsG]=eig((gram+gram')/2,'vector');
top=max(eigsG);
weights=V'*rhs;
lambda_grid=10.^(-10:-2);
validation_error=zeros(numel(lambda_grid),1);
Gvalid=dipole_matrix(xyz(idx_valid,:),sites);
target=B(idx_valid,:); yvalid=target(:);
for j=1:numel(lambda_grid)
    moments=V*(weights./(max(eigsG,0)+top*lambda_grid(j)));
    validation_error(j)=norm(Gvalid*moments-yvalid)/norm(yvalid);
end
[~,best]=min(validation_error);
moments=V*(weights./(max(eigsG,0)+top*lambda_grid(best)));
pred=zeros(size(B));
for k=1:512:size(B,1)
    sel=k:min(k+511,size(B,1));
    pred(sel,:)=reshape(dipole_matrix(xyz(sel,:),sites)*moments,[],3);
end
group=["training_subsample";"validation_z33.7";"test_z-34.5";"test_z23.8";"test_z43.65";"all_test_planes"];
sets=[{train};{idx_valid};arrayfun(@(z)find(abs(z_mm-z)<0.2),test_z','UniformOutput',false);{find(test)}];
scores=table(group);
for j=1:numel(sets)
    sel=sets{j}; measured=B(sel,:); estimated=pred(sel,:); delta=estimated-measured;
    scores.n_points(j)=numel(sel);
    scores.relative_vector_L2_percent(j)=100*norm(delta,'fro')/norm(measured,'fro');
    scores.vector_RMSE_uT_per_A(j)=1e6*sqrt(sum(delta(:).^2)/numel(sel));
    scores.rms_field_ratio(j)=norm(estimated,'fro')/norm(measured,'fro');
    scores.bx_relative_L2_percent(j)=100*norm(delta(:,1))/norm(measured(:,1));
    scores.by_relative_L2_percent(j)=100*norm(delta(:,2))/norm(measured(:,2));
    scores.bz_relative_L2_percent(j)=100*norm(delta(:,3))/norm(measured(:,3));
end
writetable(scores,fullfile(here,'mc_b70_scores.csv'));
writetable(table(lambda_grid',validation_error,'VariableNames',{'lambda_relative','validation_relative_L2'}), ...
    fullfile(here,'mc_b70_regularization.csv'));
dipoles=[sites*1000 reshape(moments,[],3)];
writetable(array2table(dipoles,'VariableNames',{'x_mm','y_mm','z_mm','mx_m2','my_m2','mz_m2'}), ...
    fullfile(here,'mc_b70_equivalent_sources.csv'));
test_table=array2table([find(test),scan.xyz_mm(test,:),1e6*B(test,:),1e6*pred(test,:)], ...
    'VariableNames',{'scan_row','x_mm','y_mm','z_mm','measured_Bx_uT_A','measured_By_uT_A','measured_Bz_uT_A', ...
    'predicted_Bx_uT_A','predicted_By_uT_A','predicted_Bz_uT_A'});
writetable(test_table,fullfile(here,'mc_b70_heldout_points.csv'));
fprintf('sites=%d train=%d validation=%d test=%d selected_relative_lambda=%g\n',size(sites,1),numel(train),numel(idx_valid),nnz(test),lambda_grid(best));
disp(scores)
end

function sites=source_sites(mesh,observers)
% For each vertical ray through the STL, use its first and last intersection
% with the casing to define two internal source layers. Offsets avoid mesh edges.
vertices=mesh.Points; tri=mesh.ConnectivityList;
a=vertices(tri(:,1),:); u=vertices(tri(:,2),:)-a; v=vertices(tri(:,3),:)-a;
detxy=u(:,1).*v(:,2)-u(:,2).*v(:,1);
low=min(vertices,[],1); high=max(vertices,[],1); sites=zeros(0,3);
for x=low(1)+5.13:12:high(1)
    for y=low(2)+5.37:12:high(2)
        dx=x-a(:,1); dy=y-a(:,2);
        p=(dx.*v(:,2)-dy.*v(:,1))./detxy;
        q=(u(:,1).*dy-u(:,2).*dx)./detxy;
        inside=abs(detxy)>1e-10 & p>0 & q>0 & p+q<1;
        zz=sort(a(inside,3)+p(inside).*u(inside,3)+q(inside).*v(inside,3));
        if mod(numel(zz),2)~=0, continue; end
        for k=1:2:numel(zz)
            for f=[.3 .7]
                candidate=[x y zz(k)+f*(zz(k+1)-zz(k))]/1000;
                if min(sum((observers-candidate).^2,2))>(.006)^2
                    sites(end+1,:)=candidate; %#ok<AGROW>
                end
            end
        end
    end
end
assert(~isempty(sites))
end

function G=dipole_matrix(points,sites)
% maps [mx;my;mz] to [Bx;By;Bz], component-major columns and rows, SI units.
r={points(:,1)-sites(:,1)',points(:,2)-sites(:,2)',points(:,3)-sites(:,3)'};
r2=r{1}.^2+r{2}.^2+r{3}.^2;
inverse=1e-7*r2.^(-1.5);
G=zeros(3*size(points,1),3*size(sites,1));
for i=1:3
    for j=1:3
        block=inverse.*(3*r{i}.*r{j}./r2-double(i==j));
        rows=(i-1)*size(points,1)+(1:size(points,1));
        cols=(j-1)*size(sites,1)+(1:size(sites,1));
        G(rows,cols)=block;
    end
end
end

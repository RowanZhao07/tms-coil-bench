function run_validation(coils,root,useGPU)
% Independent closed-form circular-axis reference + structural checks.
out=fullfile(root,'results','validation');
R=0.05; I=1000; z=(0:0.005:0.1)'; p=[zeros(size(z)),zeros(size(z)),z];
ana=4*pi*1e-7*I*R^2./(2*(R^2+z.^2).^(3/2));
conv=table();
for n=[48 96 192 384 768]
    t=(0:n)'*2*pi/n; q=[R*cos(t),R*sin(t),zeros(size(t))];
    s=[q(1:end-1,:) q(2:end,:) ones(n,1)];
    b=field_segments(s,p,I); rel=max(vecnorm(b-[0*z,0*z,ana],2,2)./ana);
    conv=[conv;table(n,rel,'VariableNames',{'segments','max_axis_relative_error'})]; %#ok<AGROW>
end
assert(conv.max_axis_relative_error(end)<2e-5,'Circle analytic validation failed.');
writetable(conv,fullfile(out,'circle_analytic_convergence.csv'));
writetable(array2table([z ana b],'VariableNames',{'z_m','analytic_Bz_T','Bx_T','By_T','Bz_T'}),fullfile(out,'circle_analytic.csv'));
% Geometry discretization changes, distinct from exact/midpoint integration.
c96=build_coils(96,root); c384=build_coils(384,root);
testP=[0 0 .01;0 0 .02;0 0 .03;0 0 .05;-.03 .02 .02;.03 -.02 .05];
gt=table();
for k=1:numel(coils)
    b96=field_segments(c96(k).segments,testP,I);
    b192=field_segments(coils(k).segments,testP,I);
    b384=field_segments(c384(k).segments,testP,I);
    e96=norm(b96-b384,'fro')/norm(b384,'fro'); e192=norm(b192-b384,'fro')/norm(b384,'fro');
    refineApplicable=~strcmp(coils(k).id,'h1');
    gt=[gt;table(string(coils(k).id),e96,e192,refineApplicable,'VariableNames',{'id','relative96_vs384','relative192_vs384','geometry_refinement_applicable'})]; %#ok<AGROW>
end
writetable(gt,fullfile(out,'geometry_convergence.csv'));
% Linearity and direction are essential for superposition among winding groups.
s=coils(1).segments; b=field_segments(s,testP,I);
linearity=norm(field_segments(s,testP,2*I)-2*b,'fro')/norm(b,'fro');
rev=s; rev(:,1:3)=s(:,4:6);rev(:,4:6)=s(:,1:3);
reversal=norm(field_segments(rev,testP,I)+b,'fro')/norm(b,'fro');
assert(linearity<1e-12 && reversal<1e-10);
gpuError=NaN; gpuStatus='not_requested';
if useGPU
    bg=field_segments(s,testP,I,'exact',true);gpuError=norm(bg-b,'fro')/norm(b,'fro');
    assert(gpuError<1e-10);gpuStatus='executed';
end
% One source-document FOE case: radii 4:8:44 mm, center separation assumed 88 mm.
radii=(4:8:44)/1000; n=768; t=(0:n)'*2*pi/n; s=[];
for side=[-1 1]
    for r=radii
        q=[side*.044+r*cos(t),r*sin(t),0*t];
        s=[s;q(1:end-1,:),q(2:end,:),side*ones(n,1)]; %#ok<AGROW>
    end
end
p=[zeros(5,2),(1:5)'/100]; b=field_segments(s,p,I);
writetable(array2table([p b vecnorm(b,2,2)],'VariableNames',{'x_m','y_m','z_m','Bx_T','By_T','Bz_T','Bnorm_T'}),fullfile(out,'document_foe_case.csv'));
checks=struct('linearity_relative_error',linearity,'direction_reversal_relative_error',reversal,'gpu_status',gpuStatus,'gpu_relative_error',gpuError,'document_case_center_separation_m',.088,'document_case_center_separation_status','assumed outer loops tangent; source R0 not numeric','circle_validation_pass',true);
fid=fopen(fullfile(out,'checks.json'),'w');fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true));fclose(fid);
end

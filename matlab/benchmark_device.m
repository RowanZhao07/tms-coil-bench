function benchmark_device
% Compare CPU/GPU only when the toolbox and physical device are available.
root=fileparts(fileparts(mfilename('fullpath')));c=build_coils(96,root);s=c(1).segments;
[x,y]=meshgrid(linspace(-.1,.1,41));p=[x(:),y(:),.02*ones(numel(x),1)];
tic;b=field_segments(s,p,1000,'exact',false);cpu=toc;
r=struct('cpu_s',cpu,'gpu_s',NaN,'gpu_relative_error',NaN,'gpu_status','unavailable');
try
    g=gpuDevice;field_segments(s,p(1:5,:),1000,'exact',true);wait(g);
    tic;bg=field_segments(s,p,1000,'exact',true);wait(g);r.gpu_s=toc;
    r.gpu_relative_error=norm(bg-b,'fro')/norm(b,'fro');r.gpu_name=g.Name;
    assert(r.gpu_relative_error<1e-10);
    r.gpu_status='executed';
catch e
    r.gpu_status='failed_or_unavailable';
    r.gpu_message=e.message;
end
fid=fopen(fullfile(root,'results','device_benchmark.json'),'w');fprintf(fid,'%s',jsonencode(r,PrettyPrint=true));fclose(fid);
disp(r);
end

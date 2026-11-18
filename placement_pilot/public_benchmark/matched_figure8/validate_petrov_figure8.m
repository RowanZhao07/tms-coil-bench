function scores = validate_petrov_figure8()
% Recompute the 1% MO phantom experiment independently of the published fields.
here = fileparts(mfilename('fullpath'));
study = fileparts(fileparts(fileparts(here)));
addpath(fullfile(study,'matlab'));
cache = fullfile(study,'tmp','petrov_benchmark');
if ~isfolder(cache), mkdir(cache); end
revisions = ["5d65b6f2dd3a6658748e59c02b5fa8182b998d5a","9e8de2a708b1d0c1d527f5fb30199ce8a1354771"];
shaA = ["00e831c3f4033da39d303da6b81fd2da5163f845c5de6cd07417cacb2faf8a97","b4fa19d52c31d998c35b456e085df104f89d5e6f725ff4058c010efd147e76cd"];
shaJ = ["9ace8316c882dfde47b4796f9b827a8c54f4dec9370a10948a47611c6b45e4ee","0cb08be6d2e2b776a612c54f5fdc661bfff31290ecedf8e8e77bc7e8e410aaa7"];
shaMaps = "6f8276720f7de75aec05e68f018693f0bdf5b905507a8210b4ac2410d931ed6d";

% set_coil.srn m45 -> m42 -> m8 -> m3. A is the full external affine.
Rx = [1 0 0 0;0 cosd(-4) -sind(-4) 0;0 sind(-4) cosd(-4) 0;0 0 0 1];
lift = eye(4); lift(3,4) = 15;
for v = 1:2
    maps = read_reference(cache,revisions(v),"error_metrics/data.mat",shaMaps);
    alignment = read_reference(cache,revisions(v),"navigator_coregistration/coil_capture_align.mat",shaA(v));
    capture = read_reference(cache,revisions(v),"navigator_coregistration/coil_capture_jacobi.mat",shaJ(v));
    poses(:,:,v) = alignment.A * capture.J * lift * Rx; %#ok<AGROW>
    captures(:,:,v) = capture.J; %#ok<AGROW>
end
assert(isequal(captures(:,:,1),captures(:,:,2)));
% Counterfactual: change only the scanner-z translation of the fully composed pose.
poses(:,:,3) = poses(:,:,1);
poses(3,4,3) = poses(3,4,2);
conditions = ["v1.0","v1.1","v1.0_with_v1.1_z_only"];
% m32 is a 161 x 161 x 2 node grid. m29 maps it from indices to RAS mm.
mask = rot90(maps.data_singlecoil)~=0;
assert(nnz(mask)==10622);
[i,j,k] = ndgrid(0:160,0:160,0:1);
needed = conv2(double(maps.data_singlecoil~=0),ones(2),'full')>0;
needed = repmat(needed,1,1,2);
observers = 1e-3*[85.033-i(needed),53.931+2*k(needed),-85.401+j(needed)];

% Table 2 specifies the equivalent current in mA for each model, at 1% MO.
% Source code and SCIRun module parameters specify a 4 mm gap (centres +/-46).
models = ["BSM-811";"BSM-819";"BSM-879"];
publishedFields = ["data_singlecoil";"data_multicoil";"data_layercoil"];
currents = [25.032;3.576;0.511];
layers = [1;1;7];
audit = table();
measured_T = maps.data_MRI1;
author_T = rot90(maps.data_multicoil);
prediction_T = nan([size(mask),3]);
registered_wires_mm = cell(1,2);
for n=1:numel(models)
  for v=1:(2 + (n==2))
    segments = make_segments(poses(:,:,v),n==1,layers(n));
    if n==2 && v<=2
        registered_wires_mm{v} = 1000*segments(:,1:6);
    end
    field = field_segments(segments,observers,currents(n),'exact',false);
    nodes = nan(size(i));
    nodes(needed) = field(:,3);
    % solve_coil_field.srn interpolates the eight nodal values to each cell.
    cells = (nodes(1:end-1,1:end-1,:)+nodes(2:end,1:end-1,:) ...
        +nodes(1:end-1,2:end,:)+nodes(2:end,2:end,:))/4;
    predicted = rot90(mean(cells,3));
    if n==2, prediction_T(:,:,v) = predicted; end
    refs = {rot90(maps.(publishedFields(n))),maps.data_MRI1};
    names = ["author_SCIRun";"phantom_MRI"];
    for r=1:2
        target = refs{r}(mask); value = predicted(mask); err=value-target;
        cc = corrcoef(target,value);
        row = table(models(n), names(r),size(segments,1),currents(n),nnz(mask),...
            100*norm(err)/norm(target),1e6*sqrt(mean(err.^2)),cc(1,2),norm(value)/norm(target),...
            'VariableNames',{'model','reference','segments','equivalent_current_A',...
            'pixels','relative_L2_percent','RMSE_uT','signed_Pearson_r','rms_amplitude_ratio'});
        row.condition = conditions(v);
        audit=[audit;row]; %#ok<AGROW>
    end
  end
end
scores = audit(audit.condition=="v1.1",:);
writetable(audit,fullfile(here,'placement_version_comparison.csv'));
writetable(scores,fullfile(here,'matched_figure8_scores.csv'));
reference = ["v1.0";"v1.1"];
origin = squeeze(poses(1:3,4,1:2))';
writetable(table(reference,origin(:,1),origin(:,2),origin(:,3), ...
    'VariableNames',{'version','origin_x_mm','origin_y_mm','origin_z_mm'}),fullfile(here,'coil_origins.csv'));
save(fullfile(here,'figure8_slice_results.mat'),'mask','measured_T','author_T','prediction_T','conditions','poses','registered_wires_mm');
disp(scores)
end

function data=read_reference(cache,revision,relative,digest)
folder = fullfile(cache,revision);
if ~isfolder(folder), mkdir(folder); end
[~,name,ext] = fileparts(relative);
path = fullfile(folder,name+ext);
if ~isfile(path)
    url = "https://raw.githubusercontent.com/pip010/plosone-data-pub01/"+revision+"/phantom_experiments/"+relative;
    temporary = string(tempname(folder))+".mat";
    websave(temporary,url,weboptions('Timeout',45));
    verify_digest(temporary,digest);
    movefile(temporary,path,'f');
else
    verify_digest(path,digest);
end
data=load(path);
end

function verify_digest(path,expected)
fid = fopen(path,'rb');
assert(fid>=0);
cleanup = onCleanup(@() fclose(fid));
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(fread(fid,Inf,'*uint8'));
hex = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
assert(strcmp(hex,expected),'Petrov input SHA-256 mismatch: %s',path);
end

function S=make_segments(transform,singleLoop,layers)
% SingleloopCoilgen and MultiloopsCoilgen in SCIRun4plus v4.7.2.
% Units are mm until after applying the scanner-space transform.
S=zeros(0,7);
for side=[-1 1]
    for layer=0:layers-1
        z=layer-floor(layers/2);
        cx=46*side;
        if singleLoop
            a=linspace(0,2*pi,93)'; % Paper Table 1, 92 segments per wing.
            pts=[cx+44*cos(a),44*sin(a),repmat(z,size(a))];
        else
            lod=6;
            if layers==7, lod=4; end
            pts=zeros(0,3);
            if side<0
                for w=0:8
                    pts=[pts;arc(cx,z,26+w*2,0,w,lod);arc(cx+1,z,27+w*2,pi,w,lod)]; %#ok<AGROW>
                end
                pts=[pts;cx+44,0,z];
            else
                for w=9:-1:1
                    pts=[pts;arc(cx,z,26+w*2,pi,w,lod);arc(cx+1,z,25+w*2,0,w,lod)]; %#ok<AGROW>
                end
                pts=[pts;cx-26,0,z];
            end
        end
        xyz = [pts,ones(size(pts,1),1)]*transform';
        xyz = 1e-3*xyz(:,1:3);
        S=[S;xyz(1:end-1,:),xyz(2:end,:),repmat(-side,size(pts,1)-1,1)]; %#ok<AGROW>
    end
end
end

function xyz=arc(cx,z,radius,start,w,lod)
% GenPointsCircular omits the final point on each half-circle.
n=lod*(8+w);
a=start+pi*(0:n-1)'/n;
xyz=[cx+radius*cos(a),radius*sin(a),repmat(z,n,1)];
end

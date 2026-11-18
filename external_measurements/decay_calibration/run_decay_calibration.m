function scores = run_decay_calibration()
% Four near-depth readings calibrate CCD fields. Far observations are scoring only.
here = fileparts(mfilename('fullpath'));
parent = fileparts(here);
fr = readtable(fullfile(parent,'lyon','air_observations.csv'),'TextType','string');
de = readtable(fullfile(parent,'regensburg','hall_fig2.csv'),'TextType','string');
coils = ["Cool-B65","Cool-D-B80","C-B60"];
hashes = ["7ef96c5123e4c8b0a393177cff9e3401637025d3ceb22c1ddb2e0d340f1dacfe", ...
    "547027e48438ab67f6577bc80984306276b9bb93645368e83885d0c1e08ce637", ...
    "35913741feee6c33d626213f7c3245966d12b317f3aa87d0b9fdf4823b85eb44"];
methods = ["one_anchor","gain_only","free_offset","shrunken_offset","guarded_offset"];
scores = table(); fits = table(); points = table(); folds = table(); probe = table();
for k = 1:numel(coils)
    coil = coils(k);
    normalProbe = coil == "C-B60";
    if normalProbe
        loc = "Regensburg";
        obs = de(de.coil == coil,:);
        obsT = obs.B_peak_T;
    else
        loc = "Lyon";
        obs = fr(fr.coil == coil,:);
        obsT = obs.B_T;
    end
    depth = obs.distance_mm;
    sources = get_sources(fileparts(parent),coil,hashes(k),normalProbe);
    curve = @(d) dipole_scan(sources,d,normalProbe);
    near = ismember(depth,[10 15 20 25]);
    assert(nnz(near)==4 && all(obsT(near)>0));
    z = depth(near); y = obsT(near);
    mag = dipole_scan(sources,z,false);
    axial = dipole_scan(sources,z,true);
    probe = [probe;table(repmat(coil,numel(z),1),z,mag,axial,axial./mag, ...
        'VariableNames',{'coil','depth_mm','CCD_norm_T_per_A','CCD_abs_Bz_T_per_A', ...
        'fraction_normal'})]; %#ok<AGROW>
    fitGain = fit_profile(curve,z,y,0);
    fitFree = fit_profile(curve,z,y,1);
    fitShrunk = fit_profile(curve,z,y,2);
    e = zeros(4,2);
    for f = 1:4
        rest = true(4,1); rest(f) = false;
        gF = fit_profile(curve,z(rest),y(rest),0);
        sF = fit_profile(curve,z(rest),y(rest),2);
        e(f,:) = log([gF(1)*curve(z(f)), sF(1)*curve(z(f)+sF(2))]./y(f));
        folds = [folds;table(loc,coil,z(f),gF(1),sF(1),sF(2),e(f,1),e(f,2), ...
            'VariableNames',{'site','coil','left_out_mm','fixed_gain','offset_gain', ...
            'offset_mm','fixed_log_residual','shrunken_log_residual'})]; %#ok<AGROW>
    end
    cv = sqrt(mean(e.^2,1));
    allowShift = cv(2) < cv(1);
    guarded = fitGain;
    if allowShift, guarded = fitShrunk; end
    pars = [y(z==10)/curve(10),0;fitGain;fitFree;fitShrunk;guarded];
    common = ismember(depth,30:5:60);
    assert(nnz(common)==7);
    if normalProbe
        extra = depth>=41 & depth<=70; extraName = "41_70_dense";
    else
        extra = depth==140; extraName = "140_endpoint";
    end
    subset = near | common | extra;
    for j = 1:numel(methods)
        gain = pars(j,1); offset = pars(j,2);
        train = near;
        if j==1, train = depth==10; end
        p = nan(size(depth));
        p(subset) = gain*curve(depth(subset)+offset);
        % Distance is identified only through depth variation of the log-slope.
        dDiff = .02;
        slope = (log(curve(z+offset+dDiff))-log(curve(z+offset-dDiff)))/(2*dDiff);
        weakDirectionMM = .05 / norm(slope-mean(slope));
        fits = [fits;table(loc,coil,methods(j),nnz(train),gain,offset, ...
            sqrt(mean(log(p(train)./obsT(train)).^2)),cv(1),cv(2),allowShift, ...
            weakDirectionMM,'VariableNames',{'site','coil','method','calibration_n', ...
            'effective_gain','effective_offset_mm','train_log_rmse','gain_loo_log_rmse', ...
            'shrunken_loo_log_rmse','gate_passes','local_offset_scale_four_anchors_mm'})]; %#ok<AGROW>
        regions = {common,extra}; names = ["30_60_common",extraName];
        for i = 1:2
            t = regions{i}; residual = p(t)-obsT(t);
            scores = [scores;table(loc,coil,methods(j),names(i),nnz(t), ...
                100*norm(residual)/norm(obsT(t)),1000*sqrt(mean(residual.^2)), ...
                1000*mean(residual),100*mean(residual./obsT(t)), ...
                'VariableNames',{'site','coil','method','test_region','n','L2_percent', ...
                'RMSE_mT','mean_signed_mT','mean_signed_percent'})]; %#ok<AGROW>
        end
        n = nnz(subset);
        points = [points;table(repmat(loc,n,1),repmat(coil,n,1),repmat(methods(j),n,1), ...
            depth(subset),train(subset),common(subset),extra(subset),obsT(subset),p(subset), ...
            'VariableNames',{'site','coil','method','depth_mm','used_in_fit','common_test', ...
            'extra_test','measured_T','predicted_T'})]; %#ok<AGROW>
    end
end
writetable(scores,fullfile(here,'test_scores.csv'));
writetable(fits,fullfile(here,'calibration.csv'));
writetable(folds,fullfile(here,'near_leave_one_out.csv'));
writetable(points,fullfile(here,'predictions.csv'));
writetable(probe,fullfile(here,'probe_components.csv'));
draw_results(here,scores,coils,methods);
disp(scores(:,{'coil','method','test_region','L2_percent','RMSE_mT'}));
end

function parameters = fit_profile(curve,d,b,kind)
% For a fixed offset, eliminate log(gain) analytically.
logB = log(b); delta = 0;
if kind>0
    penalty = (kind==2)*(.05/5)^2;
    loss = @(x) norm(detrend(logB-log(curve(d+x)),0))^2 + penalty*x*x;
    % A small grid also detects the ends and multiple one-dimensional minima.
    grid = -5:.5:5; val = arrayfun(loss,grid);
    candidates = [-5,5];
    for n = find(val(2:end-1)<=val(1:end-2) & val(2:end-1)<=val(3:end))+1
        candidates(end+1) = fminbnd(loss,grid(n-1),grid(n+1)); %#ok<AGROW>
    end
    [~,best] = min(arrayfun(loss,candidates)); delta = candidates(best);
end
gain = exp(mean(logB-log(curve(d+delta))));
parameters = [gain,delta];
end

function field = dipole_scan(s,d,normalProbe)
field = zeros(numel(d),1);
for i = 1:numel(d)
    r = [0,0,d(i)/1000] - s(:,1:3);
    r2 = sum(r.^2,2); m = s(:,4:6);
    vector = 1e-7*sum((3*r.*sum(m.*r,2)./r2-m)./r2.^1.5,1);
    if normalProbe
        field(i) = abs(vector(3));
    else
        field(i) = norm(vector);
    end
end
end

function sources = get_sources(root,coil,expected,normalProbe)
if normalProbe
    file = fullfile(root,'tmp','regensburg','C-B60.ccd');
else
    file = fullfile(root,'tmp','lyon',"MagVenture_"+coil+".ccd");
end
if ~isfile(file)
    if ~isfolder(fileparts(file)), mkdir(fileparts(file)); end
    rev = "fde509f6f8de38a4b03947944cf7ee93b8bebabe";
    url = "https://raw.githubusercontent.com/simnibs/simnibs/"+rev+ ...
        "/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_"+coil+".ccd";
    websave(file,url,weboptions('Timeout',60));
end
fid = fopen(file,'rb'); assert(fid>=0);
guard = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('SHA-256'); md.update(fread(fid,Inf,'*uint8'));
assert(strcmpi(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]),expected));
clear guard
sources = readmatrix(file,'FileType','text','NumHeaderLines',3);
end

function draw_results(here,scores,coils,methods)
fig = figure('Color','w','Visible','off','Position',[50 80 1030 390]);
l = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
labels = ["10 mm only","Four-point gain","Gain + offset","Shrunken offset","Guarded offset"];
common = scores(scores.test_region=="30_60_common",:);
ax = nexttile(l); M = zeros(numel(coils),numel(methods));
for c = 1:numel(coils)
    for j = 1:numel(methods)
        M(c,j) = common.L2_percent(common.coil==coils(c) & common.method==methods(j));
    end
end
bar(ax,M,'grouped'); xticks(ax,1:3); xticklabels(ax,coils);
ylabel(ax,'Held-out relative L2 error (%)');title(ax,'Matched distances (30–60 mm)');
ax.YGrid='on'; legend(ax,labels,'FontSize',8,'Box','off','Location','northwest');
ax = nexttile(l);
for c = 1:2
    for j = 1:numel(methods)
        ix = scores.coil==coils(c) & scores.method==methods(j) & scores.test_region=="140_endpoint";
        M(c,j) = scores.mean_signed_percent(ix);
    end
end
bar(ax,M(1:2,:),'grouped'); xticks(ax,1:2); xticklabels(ax,coils(1:2)); ax.YGrid='on';
ylabel(ax,'Signed error (%)'); title(ax,'Far endpoint (140 mm, France)');
exportgraphics(fig,fullfile(here,'decay_ablation.png'),'Resolution',165); close(fig);
end

function run_finland(cache)
% Audit Matsuda et al. archive ID 15398706. Inputs use archive units as given.
here=fileparts(mfilename('fullpath'));
if nargin==0, cache=fullfile(here,'..','..','tmp','independent_lab'); end
if ~isfolder(cache), mkdir(cache); end
options=delimitedTextImportOptions('NumVariables',4,'Delimiter',',');
options.DataLines=[2 Inf];
options.VariableNames={'filename','bytes','md5','url'};
options.VariableTypes={'string','double','string','string'};
manifest=readtable(fullfile(here,'finland_source_manifest.csv'),options);
probePeaks=nan(100,30);
coilPeaks=nan(100,5,2);
for i=1:height(manifest)
    name=manifest.filename(i); path=fullfile(cache,name);
    if ~isfile(path)
        temp=string(tempname(cache));
        websave(temp,manifest.url(i),weboptions('Timeout',90));
        checkFile(temp,manifest.bytes(i),manifest.md5(i));
        movefile(temp,path);
    end
    checkFile(path,manifest.bytes(i),manifest.md5(i));
    if ~endsWith(name,'.txt'), continue; end
    lines=readlines(path,'EmptyLineRule','skip');
    ids=str2double(lines(1:2:end));
    assert(isequal(ids,(0:99)'));
    profile=zeros(100,1);
    for j=1:100
        trace=sscanf(lines(2*j),'%f');
        assert(numel(trace)==400 && all(isfinite(trace)));
        profile(j)=max(abs(trace-median(trace(1:8))));
    end
    coil=regexp(name,'^ch(\d)_coil(\d)_','tokens','once');
    if ~isempty(coil)
        coilPeaks(:,str2double(coil{2}),str2double(coil{1})+1)=profile;
    else
        target=regexp(name,'robot_target_(\d+)\.txt$','tokens','once');
        assert(~isempty(target));
        probePeaks(:,str2double(target{1}))=profile;
    end
end
assert(all(isfinite(probePeaks),'all') && all(isfinite(coilPeaks),'all'))
typical=median(probePeaks,2);
targetScores=table((1:30)',max(probePeaks,[],1)', ...
    100*vecnorm(probePeaks-typical)'/norm(typical), ...
    'VariableNames',{'target_id','max_baseline_adjusted_archive_units','L2_vs_median_profile_pct'});
profileCorrelations=corrcoef([typical,probePeaks]);
targetScores.Pearson_vs_median=profileCorrelations(1,2:end)';
writetable(targetScores,fullfile(here,'finland_target_profiles.csv'));
coilSummary=table;
for n=1:5
    for k=1:2
        pk=coilPeaks(:,n,k);
        [value,ix]=max(pk);
        coilSummary=[coilSummary; table(n,k-1,value,ix-1,'VariableNames', ...
            {'coil','channel','max_baseline_adjusted_archive_units','peak_position_id'})]; %#ok<AGROW>
    end
end
writetable(coilSummary,fullfile(here,'finland_coil_channel_peaks.csv'));
% The workbook uses rows 2:32 in its own original AVERAGE and STDEV.S formulas.
values=readmatrix(fullfile(cache,'targeting+mtms.xlsx'),'Sheet',1,'Range','B2:D32');
assert(isequal(size(values),[31 3]) && all(isfinite(values),'all'));
workbook=table(["centroid_x_mm";"centroid_y_mm";"peak_E_V_per_m"],mean(values,1)',std(values,0,1)', ...
    'VariableNames',{'quantity','mean','sample_SD'});
writetable(workbook,fullfile(here,'finland_targeting_workbook.csv'));
fig5=readtable(fullfile(here,'finland_fig5_peaks.csv'));
fig5.relative_difference_pct=100*abs(fig5.with_robot_V_per_m-fig5.without_robot_V_per_m)./fig5.without_robot_V_per_m;
% Each printed value has half of its 0.1 V/m display step as rounding margin.
a=fig5.without_robot_V_per_m; b=fig5.with_robot_V_per_m;
fig5.rounding_lower_pct=100*max(0,(b-.05)./(a+.05)-1);
fig5.rounding_upper_pct=100*max(abs((b+.05)./(a-.05)-1),abs((b-.05)./(a+.05)-1));
writetable(fig5,fullfile(here,'finland_attachment_comparison.csv'));
save(fullfile(here,'finland_profiles.mat'),'probePeaks','coilPeaks','typical','values');
f=figure('Visible','off','Color','w','Position',[0 0 1030 395]);
t=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(t); imagesc(ax,1:30,0:99,probePeaks); axis(ax,'xy');
c=colorbar(ax); ylabel(c,'Peak, archived signal units');
xlabel(ax,'Robot target ID'); ylabel(ax,'Probe position ID'); title(ax,'Channel 0, baseline-adjusted');
ax=nexttile(t); scatter(ax,values(:,1),values(:,2),42,values(:,3),'filled','MarkerEdgeColor',[.25 .25 .25]);
c=colorbar(ax); ylabel(c,'Peak E (V/m)'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'Centroid x offset (mm)');ylabel(ax,'Centroid y offset (mm)');title(ax,'Workbook, 31 entries');
exportgraphics(f,fullfile(here,'finland_repositioning.png'),'Resolution',165);close(f);
disp(workbook);disp(fig5); disp([median(targetScores.L2_vs_median_profile_pct),min(targetScores.Pearson_vs_median)])
end

function checkFile(path,bytes,expected)
info=dir(path); assert(info.bytes==bytes);
fid=fopen(path,'rb'); assert(fid>=0); guard=onCleanup(@()fclose(fid));
digest=java.security.MessageDigest.getInstance('MD5');
digest.update(fread(fid,Inf,'*uint8'));
actual=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
assert(strcmp(actual,expected),'Archive digest mismatch for %s.',path);
end

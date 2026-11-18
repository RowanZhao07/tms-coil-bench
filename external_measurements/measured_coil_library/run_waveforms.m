function scores = run_waveforms
% Compare each published pickup waveform with the publisher's own paired fit.
here=fileparts(mfilename('fullpath'));
cache=fullfile(here,'..','..','tmp','measured_coil_library');
manifest=readtable(fullfile(here,'waveform_manifest.csv'),'TextType','string');
scores=table();
for k=1:height(manifest)
    path=fullfile(cache,manifest.filename(k));
    fid=fopen(path,'rt'); assert(fid>0);
    names=split(strtrim(fgetl(fid)),char(9)); fclose(fid);
    assert(numel(names)==3 && names(3)==names(2)+"_fit");
    d=readmatrix(path,'FileType','text','NumHeaderLines',1,'Delimiter',char(9));
    t=d(:,1); y=d(:,2); fitted=d(:,3);
    assert(all(isfinite(d),'all') && all(diff(t)>0));
    cc=corrcoef(y,fitted);
    row=table(erase(manifest.filename(k),".tsv"),names(2),numel(t), ...
        1e6*median(diff(t)),1e6*t(1),1e6*t(end),100*norm(fitted-y)/norm(y), ...
        norm(fitted)/norm(y),cc(1,2), ...
        'VariableNames',{'coil','stimulator','samples','dt_us','tmin_us','tmax_us', ...
        'relative_L2_percent','rms_ratio','pearson_r'});
    scores=[scores;row]; %#ok<AGROW>
end
writetable(scores,fullfile(here,'waveform_scores.csv'));
fprintf('N=%d timepoints=%d median relative L2=%.4f%%, range %.4f to %.4f%%\n', ...
    height(scores),sum(scores.samples),median(scores.relative_L2_percent), ...
    min(scores.relative_L2_percent),max(scores.relative_L2_percent));
end

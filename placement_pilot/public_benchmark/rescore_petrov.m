function results = rescore_petrov(file)
% Re-score the authors' predictions against their measured Bz phantom slice.
% Source and their rotation / ROI conventions are recorded in README.md.
here = fileparts(mfilename('fullpath'));
if nargin == 0
    file = fullfile(here,'..','..','tmp','petrov_v1_metrics.mat');
    if ~isfile(file)
        if ~isfolder(fileparts(file))
            mkdir(fileparts(file));
        end
        websave(file,['https://raw.githubusercontent.com/pip010/' ...
            'plosone-data-pub01/v1.0/phantom_experiments/error_metrics/data.mat']);
    end
end
maps = load(file);
measured = maps.data_MRI1;
mask = rot90(maps.data_singlecoil) ~= 0;
assert(isequal(size(measured),[160,160]) && nnz(mask)==10622);
measured = measured(mask);
assert(all(isfinite(measured)));
names = ["data_singlecoil"; "data_multicoil"; "data_layercoil"];
models = ["BSM-811"; "BSM-819"; "BSM-879"];
L2_percent = zeros(3,1);
RMSE_uT = zeros(3,1);
for k = 1:numel(names)
    simulated = rot90(maps.(names(k)));
    residual = simulated(mask)-measured;
    assert(all(isfinite(residual)));
    L2_percent(k) = 100*norm(residual)/norm(measured);
    RMSE_uT(k) = sqrt(mean(residual.^2))*1e6;
end
results = table(models,repmat(nnz(mask),3,1),L2_percent,RMSE_uT, ...
    'VariableNames',{'published_model','pixels','L2_vs_MRI_percent','RMSE_uT'});
writetable(results,fullfile(here,'published_predictions_vs_mri.csv'));
disp(results)
end

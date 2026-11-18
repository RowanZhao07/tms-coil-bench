function [results, summary] = run_italy()
% Public peak B measurements for the Magstim HP 90 mm 9784-00 circular coil.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','..','matlab'))
results = readtable(fullfile(here,'italy_measured.csv'));
observers_m = [results.radial_horizontal_cm, results.height_cm - 113, results.axial_cm] / 100;
current_A = 6800; % Peak current used in D'Agostino et al. for the 70% MSO comparison.
% Place 14 centerlines at equal intervals in the stated winding annulus.
radii_m = 0.035 + ((1:14)-0.5)*(0.061-0.035)/14;
spec_radii_m = 0.0335 + ((1:14)-0.5)*(0.063-0.0335)/14;
field = @(r,n,p) 1e6*vecnorm(field_segments(turns(r,n),p,current_A,'exact',false),2,2);
results.predicted_peak_uT = field(radii_m,1536,observers_m);
results.original_annulus_uT = field(spec_radii_m,1536,observers_m);
results.mesh_change_pct = 100*(field(radii_m,768,observers_m)./results.predicted_peak_uT-1);
results.delta_minus5mm_uT = field(radii_m,1536,observers_m-[0,0,.005]);
results.delta_plus5mm_uT = field(radii_m,1536,observers_m+[0,0,.005]);
axial = observers_m(:,1)==0 & observers_m(:,2)==0;
exact_axis = 1e6*sum((4*pi*1e-7*current_A/2)*radii_m.^2 ./ ...
    (radii_m.^2 + observers_m(axial,3).^2).^(3/2),2);
assert(max(abs(results.predicted_peak_uT(axial)./exact_axis-1)) < 1e-5);
assert(max(abs(results.mesh_change_pct)) < .002);
results.relative_error_pct = 100*(results.predicted_peak_uT./results.measured_peak_uT-1);
results.in_reported_5pct_band = abs(results.relative_error_pct)<=5;
results.published_error_pct = 100*(results.published_sim4life_uT./results.measured_peak_uT-1);
writetable(results,fullfile(here,'italy_results.csv'));
groups = {true(height(results),1),isfinite(results.published_sim4life_uT),~isfinite(results.published_sim4life_uT)};
summary = table(["all_seven";"five_reprinted";"two_additional_from_instrument_report"], ...
    'VariableNames',{'points'});
for k=1:numel(groups)
    ix=groups{k};
    delta=results.predicted_peak_uT(ix)-results.measured_peak_uT(ix);
    err=results.relative_error_pct(ix);
    summary.count(k)=nnz(ix);
    summary.relative_L2_percent(k)=100*norm(delta)/norm(results.measured_peak_uT(ix));
    summary.mean_abs_point_percent(k)=mean(abs(err));
    summary.max_abs_point_percent(k)=max(abs(err));
    summary.in_5pct_band(k)=nnz(abs(err)<=5);
end
writetable(summary,fullfile(here,'italy_summary.csv'));
f=figure('Visible','off','Color','w','Position',[0 0 1050 410]);
t=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(t); hold(ax,'on'); k=results.point;
semilogy(ax,k,results.measured_peak_uT,'ko','MarkerFaceColor','k');
semilogy(ax,k,results.predicted_peak_uT,'bs','LineWidth',1.3);
semilogy(ax,k,results.published_sim4life_uT,'r+','LineWidth',1.5,'MarkerSize',9);
set(ax,'YScale','log','XTick',k); xlim(ax,[.5,7.5]); grid(ax,'on');
xlabel(ax,'Instrument-report point'); ylabel(ax,'Peak B (µT)'); title(ax,'Absolute field');
legend(ax,{'Measured','MATLAB','Published Sim4Life'},'Location','northeast');
ax=nexttile(t); hold(ax,'on');
patch(ax,[.5 7.5 7.5 .5],[-5 -5 5 5],[.89 .93 .89],'EdgeColor','none');
plot(ax,k,results.relative_error_pct,'bs','LineWidth',1.3);
plot(ax,k,results.published_error_pct,'r+','LineWidth',1.5,'MarkerSize',9);
yline(ax,0,'k:'); set(ax,'XTick',k); xlim(ax,[.5,7.5]);grid(ax,'on');
xlabel(ax,'Instrument-report point'); ylabel(ax,'(Prediction − measurement) / measurement (%)');
title(ax,'Signed pointwise error');legend(ax,{'Reported ±5%','MATLAB','Published Sim4Life'},'Location','best');
exportgraphics(f,fullfile(here,'italy_comparison.png'),'Resolution',165);close(f);
disp(results(:,{'point','measured_peak_uT','predicted_peak_uT','relative_error_pct','published_error_pct'}));
disp(summary)
end

function s=turns(radii,segmentsPerTurn)
a=(0:segmentsPerTurn)'*(2*pi/segmentsPerTurn);
x=cos(a)*radii; y=sin(a)*radii;
z=zeros((numel(a)-1)*numel(radii),1);
s=[reshape(x(1:end-1,:),[],1),reshape(y(1:end-1,:),[],1),z, ...
    reshape(x(2:end,:),[],1),reshape(y(2:end,:),[],1),z,ones(size(z))];
end

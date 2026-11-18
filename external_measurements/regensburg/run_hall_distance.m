function [summary,profiles] = run_hall_distance()
% Two published commercial-coil distance profiles. MATLAB, no toolbox required.
here=fileparts(mfilename('fullpath'));
obs=readtable(fullfile(here,'hall_fig2.csv'),'TextType','string');
src=read_dipoles(here);
summary=table();
profiles=table();
readout=table();
coils=["C-B60","D-B80"];
for coil=coils
    take=obs.coil==coil;
    distance=obs.distance_mm(take);
    B=obs.B_peak_T(take);
    % Split along distance instead of interleaving neighbouring points.
    training=distance>=5 & distance<=40;
    held=distance>=41 & distance<=70;
    dh=distance(training); bh=B(training);
    tau=fminbnd(@(t) sum((fit_exp(dh,bh,t,dh)-bh).^2),2,150);
    predicted=fit_exp(dh,bh,tau,distance);
    summary=[summary;score(coil,"distance_exponential",B(held),predicted(held),tau,41,70)]; %#ok<AGROW>
    % Sensitivity to a second reading rule on the same pixels, not sensor noise.
    alternate=(708.35-obs.strip_median_y_px(take))*.7/(708.35-10.90);
    t2=fminbnd(@(t) sum((fit_exp(dh,alternate(training),t,dh)-alternate(training)).^2),2,150);
    p2=fit_exp(dh,alternate(training),t2,distance);
    delta=sort(1e3*abs(B-alternate));
    r=table(coil,median(delta),delta(ceil(.95*numel(delta))),max(delta),...
        1e3*sqrt(mean((p2(held)-alternate(held)).^2)),...
        'VariableNames',{'coil','median_marker_vs_strip_mT','p95_marker_vs_strip_mT','max_marker_vs_strip_mT','alternate_exp_test_RMSE_mT'});
    readout=[readout;r]; %#ok<AGROW>
    dipole_z=nan(size(B));dipole_total=nan(size(B));
    if coil=="C-B60"
        field=dipole_field(src,distance);
        % The article reports a single Hall component, with the sensitive plane
        % parallel to the coil surface. Coordinate origin is the model centre.
        dipole_z=abs(field(:,3));
        dipole_total=vecnorm(field,2,2);
        % The stimulator's peak coil current is not reported. Use ONE distance.
        % This tests the curve shape. No orientation/gain fitted on test points.
        anchor=distance==10;
        dipole_z=dipole_z*B(anchor)/dipole_z(anchor);
        summary=[summary;score(coil,"CCD_normal_anchored_10mm",B(held),dipole_z(held),NaN,41,70)];
    end
    p=table(repmat(coil,size(B)),distance,B,predicted,dipole_z,dipole_total, ...
        'VariableNames',{'coil','distance_mm','digitized_peak_T','distance_exp_T','CCD_normal_scaled_T','CCD_norm_T_per_A'});
    profiles=[profiles;p]; %#ok<AGROW>
end
writetable(summary,fullfile(here,'held_out_metrics.csv'));
writetable(profiles,fullfile(here,'distance_profiles.csv'));
writetable(readout,fullfile(here,'readout_sensitivity.csv'));
draw_figure(here,profiles);
disp(summary);
end

function p=fit_exp(d,b,t,z)
e=exp(-d/t);scale=dot(e,b)/dot(e,e);p=scale*exp(-z/t);
end

function s=score(coil,method,y,p,tau,start_mm,end_mm)
s=table(coil,method,numel(y),start_mm,end_mm,tau, ...
    1e3*sqrt(mean((p-y).^2)),100*norm(p-y)/norm(y), ...
    'VariableNames',{'coil','method','test_points','first_test_mm','last_test_mm','decay_length_mm','RMSE_mT','relative_L2_percent'});
end

function src=read_dipoles(here)
cache=fullfile(here,'..','..','tmp','regensburg');
if ~isfolder(cache),mkdir(cache);end
file=fullfile(cache,'C-B60.ccd');
rev='fde509f6f8de38a4b03947944cf7ee93b8bebabe';
url=['https://raw.githubusercontent.com/simnibs/simnibs/',rev, ...
    '/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_C-B60.ccd'];
if ~isfile(file),websave(file,url,weboptions('Timeout',60));end
fid=fopen(file,'rb'); assert(fid>=0);
guard=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(fread(fid,Inf,'*uint8'));
sha=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
assert(strcmp(sha,'35913741feee6c33d626213f7c3245966d12b317f3aa87d0b9fdf4823b85eb44'));
clear guard
src=readmatrix(file,'FileType','text','NumHeaderLines',3);
assert(isequal(size(src),[1960,6]));
end

function b=dipole_field(src,distance)
% SI field of an ideal magnetic dipole, normalized to 1 A coil current.
mom=src(:,4:6);b=zeros(numel(distance),3);
for i=1:numel(distance)
    r=[0,0,.001*distance(i)]-src(:,1:3);
    r2=sum(r.^2,2);
    b(i,:)=1e-7*sum((3*r.*sum(mom.*r,2)./r2-mom)./r2.^1.5,1);
end
end

function draw_figure(here,prof)
f=figure('Visible','off','Color','w','Position',[0 0 1100 405]);
set(f,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei');
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
coils=["C-B60","D-B80"];
for n=1:2
    p=prof(prof.coil==coils(n),:);
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on');
    patch(ax,[40.5 70.5 70.5 40.5],[0 0 0.63 0.63],[0.92 .95 .96], ...
        'EdgeColor','none','HandleVisibility','off');
    h=plot(ax,p.distance_mm,p.digitized_peak_T,'.','MarkerSize',8,'Color',[.15 .15 .15]);
    j=plot(ax,p.distance_mm,p.distance_exp_T,'-','Color',[.12 .38 .59],'LineWidth',1.8);
    if n==1
        beyond=p.distance_mm>=10;
        k=plot(ax,p.distance_mm(beyond),p.CCD_normal_scaled_T(beyond),'--','Color',[.67 .26 .12],'LineWidth',1.4);
        legend(ax,[h j k],{'霍尔峰值（图上读取）','指数模型','等效源法向分量（10 mm 对齐）'}, ...
            'Location','northeast','Box','off');
    else
        legend(ax,[h j],{'霍尔峰值（图上读取）','指数模型'},'Location','northeast','Box','off');
    end
    xlabel(ax,'距线圈中心表面的位移（mm）'); ylabel(ax,'峰值读数（T）');
    title(ax,coils(n),'FontWeight','normal'); xlim(ax,[0 100]);ylim(ax,[0 .63]);
end
exportgraphics(f,fullfile(here,'hall_profiles.png'),'Resolution',165);close(f);
end

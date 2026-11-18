function summary=run_lyon_distance()
% Quesada air data only. Use one amplitude anchor and no position adjustment.
here=fileparts(mfilename('fullpath'));
obs=readtable(fullfile(here,'air_observations.csv'),'TextType','string');
coils=["Cool-B65","Cool-D-B80"];
sha=["7ef96c5123e4c8b0a393177cff9e3401637025d3ceb22c1ddb2e0d340f1dacfe", ...
    "547027e48438ab67f6577bc80984306276b9bb93645368e83885d0c1e08ce637"];
summary=table();
sensitivity=table();
f=figure('Visible','off','Color','w','Position',[0 0 1060 410]);
set(f,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei');
t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for i=1:2
    coil=coils(i); which=obs.coil==coil; z=obs.distance_mm(which); y=obs.B_T(which);
    src=read_source(here,coil,sha(i));
    b=local_field(src,z); v=vecnorm(b,2,2);
    anchor=(z==10); p=v*y(anchor)/v(anchor);
    obs.model_T_per_A(which)=v;
    obs.Bx_T_per_A(which)=b(:,1);obs.By_T_per_A(which)=b(:,2);obs.Bz_T_per_A(which)=b(:,3);
    obs.model_scaled_T(which)=p;
    test=z>=15 & z<=60; far=z==140;
    error=p(test)-y(test);
    % Calculate the model's far point independently. No interpolation across 60–140 mm.
    row=table(coil,sum(test),1e3*sqrt(mean(error.^2)),100*norm(error)/norm(y(test)), ...
        y(far),p(far),100*(p(far)/y(far)-1), ...
        'VariableNames',{'coil','test_points_15_60','test_RMSE_mT','test_relative_L2_percent', ...
        'observed_140_T','predicted_140_T','signed_140_error_percent'});
    summary=[summary;row]; %#ok<AGROW>
    % The active-centre offset from the contact tip is not recorded. Compare
    % fixed offsets without selecting a best shift using the experimental data.
    for offset=[-5,0,5]
        shifted=vecnorm(local_field(src,z(test)+offset),2,2);
        ref=norm(local_field(src,10+offset));
        candidate=shifted*y(anchor)/ref;
        s=table(coil,offset,100*norm(candidate-y(test))/norm(y(test)), ...
            'VariableNames',{'coil','assumed_axial_offset_mm','test_relative_L2_percent'});
        sensitivity=[sensitivity;s]; %#ok<AGROW>
    end
    ax=nexttile(t); hold(ax,'on');grid(ax,'on');
    fine_z=(0:140)'; fine=vecnorm(local_field(src,fine_z),2,2)*y(anchor)/v(anchor);
    h=plot(ax,fine_z,fine,'Color',[.17 .40 .61],'LineWidth',1.7);
    s=plot(ax,z,y,'o','Color',[.10 .10 .10],'MarkerSize',4);
    a=plot(ax,z(anchor),y(anchor),'s','MarkerFaceColor',[.77 .46 .10],'MarkerEdgeColor',[.77 .46 .10],'MarkerSize',7);
    set(ax,'YScale','log'); ylim(ax,[.002 1.5]);xlim(ax,[-4 144]);
    xlabel(ax,'探头从接触点起的位移（mm）');ylabel(ax,'磁场强度（T，对数坐标）');
    legend(ax,[s h a],{'空气台架平均读数','同型号等效源预测','10 mm 标定点'},'Box','off','Location','northeast');
    title(ax,coil,'FontWeight','normal');
end
exportgraphics(f,fullfile(here,'air_profiles.png'),'Resolution',170);close(f);
writetable(obs,fullfile(here,'air_profiles.csv'));writetable(summary,fullfile(here,'air_metrics.csv'));
writetable(sensitivity,fullfile(here,'axial_origin_sensitivity.csv'));
disp(summary);
end

function dipoles=read_source(here,coil,expected)
dir=fullfile(here,'..','..','tmp','lyon');
if ~isfolder(dir),mkdir(dir);end
name="MagVenture_"+coil+".ccd"; file=fullfile(dir,name);
url="https://raw.githubusercontent.com/simnibs/simnibs/fde509f6f8de38a4b03947944cf7ee93b8bebabe/"+...
    "simnibs/resources/coil_models/Drakaki_BrainStim_2022/"+name;
if ~isfile(file),websave(file,url,weboptions('Timeout',60));end
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));
m=java.security.MessageDigest.getInstance('SHA-256');m.update(fread(fid,Inf,'*uint8'));
assert(strcmpi(reshape(dec2hex(typecast(m.digest(),'uint8'),2)',1,[]),expected));
clear guard
dipoles=readmatrix(file,'FileType','text','NumHeaderLines',3);
end

function b=local_field(src,z)
b=zeros(numel(z),3); m=src(:,4:6);
for k=1:numel(z)
    r=[0,0,.001*z(k)]-src(:,1:3);r2=sum(r.^2,2);
    b(k,:)=1e-7*sum((3*r.*sum(m.*r,2)./r2-m)./r2.^1.5,1);
end
end

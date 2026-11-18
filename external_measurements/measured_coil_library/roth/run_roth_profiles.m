function result = run_roth_profiles
% Recalculate laboratory-scale E-field profiles from published plotted markers.
here=fileparts(mfilename('fullpath'));
result=readtable(fullfile(here,'figure_marker_pixels.csv'),'TextType','string');
for fig=[8 9]
    sel=result.figure==fig;
    % Image coordinates in the original embedded JPEG, counting from zero.
    % Endpoints are identified from labeled axis ticks, not coil curves.
    if fig==8
        xmin=250; xmax=697; depthMax_cm=10; y0=961; ytop=493; Emax=160;
    else
        xmin=264; xmax=968; depthMax_cm=8; y0=919; ytop=497; Emax=180;
    end
    result.distance_cm(sel)=(result.x_px(sel)-xmin)*depthMax_cm/(xmax-xmin);
    result.field_V_per_m(sel)=(y0-result.y_px(sel))*Emax/(y0-ytop);
end
writetable(result,fullfile(here,'digitized_e_fields.csv'));
for name=["H1_left","figure8_left"]
    data=result(result.coil==name,:);
    E=interp1(data.distance_cm,data.field_V_per_m,4.5);
    fprintf('%s E on the shared left lateral line at 4.5 cm: %.2f V/m\n',name,E);
end
fig=figure('Visible','off','Color','white','Position',[40 40 920 380]);
layout=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
for panel=1:2
    ax=nexttile(layout);
    if panel==1
        names=["H1_left","figure8_left"]; labels=["H1","Figure-8"];
        markers=["s","x"]; heading="A   Shared left measurement line";
    else
        names="figure8_normal"; labels="Figure-8"; markers="d";
        heading="B   Normal to figure-8 center";
    end
    hold(ax,'on');
    for k=1:numel(names)
        data=result(result.coil==names(k),:);
        plot(ax,data.distance_cm,data.field_V_per_m,'Marker',markers(k), ...
            'MarkerSize',6,'LineWidth',1.4,'DisplayName',labels(k));
    end
    if panel==1
        xline(ax,4.5,':','4.5 cm','HandleVisibility','off');
    end
    ax.Box='off'; ax.FontSize=10;
    xlabel(ax,'Distance along measurement line (cm)');
    ylabel(ax,'Induced electric field (V/m)');
    title(ax,heading,'FontWeight','normal'); 
    ylim(ax,[0 170]); legend(ax,'Location','northeast','Box','off');
end
exportgraphics(fig,fullfile(here,'e_field_profiles.png'),'Resolution',180);
close(fig);
end

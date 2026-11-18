function [measurements, scores] = rerun_abbasi()
% Abbasi et al., Sensors 2024, Figs 7 and 8. Only isolated measurement
% markers are inputs. Pixel coordinates in the CSV have a top-left zero origin.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','..','matlab'));
measurements = readtable(fullfile(here,'published_markers.csv'),'TextType','string');
% Centres of first and last numbered y-axis labels in the publisher PNG.
zeroPx = [1418,1014.5]; topPx = [56,114.5]; top_mT = [300,400];
outer_mm = [62,76]; turns = [9,6]; % two tangent circular wings, one effective radius
for k = 1:2
    rows = measurements.figure == 6+k;
    y = measurements.pixel_y_zero_based(rows);
    measurements.Hall_mT(rows) = (zeroPx(k)-y)*top_mT(k)/(zeroPx(k)-topPx(k));
    imageData = read_graph(here,k);
    currentsInFigure = measurements.peak_current_A(rows);
    indices = find(rows);
    for j = 1:numel(indices)
        idx = indices(j);
        x = measurements.pixel_x_zero_based(idx);
        y = measurements.pixel_y_zero_based(idx);
        ink = 166*(currentsInFigure(j)==min(currentsInFigure));
        assert(imageData(y+1,x+1)==ink,'Selected marker does not have the expected ink.');
        % Sample the author's continuous line away from the measured marker.
        offsets = [-32,32];
        if k==1 && x<600, offsets=[32,56]; end
        if k==1 && x>2200, offsets=[-56,-32]; end
        yLine = zeros(size(offsets));
        for q=1:2
            validRows = (round(topPx(k))+20):(round(zeroPx(k))-12);
            inkRows = validRows(imageData(validRows+1,x+offsets(q)+1)==ink);
            assert(~isempty(inkRows) && max(inkRows)-min(inkRows)<30);
            yLine(q)=median(inkRows);
        end
        yAtMarker=interp1(offsets,yLine,0,'linear','extrap');
        measurements.author_COMSOL_mT(idx)=(zeroPx(k)-yAtMarker)*top_mT(k)/(zeroPx(k)-topPx(k));
    end
    radius_m = outer_mm(k)/4000;
    phi = linspace(0,2*pi,385)';
    seg = zeros(0,7);
    for wing = [-1,1]
        coords = radius_m*[wing + cos(phi),sin(phi),zeros(size(phi))];
        seg = [seg;coords(1:end-1,:),coords(2:end,:), ...
            repmat(wing*turns(k),numel(phi)-1,1)]; %#ok<AGROW>
    end
    z = measurements.depth_mm(rows)/1000;
    field_T_per_A = field_segments(seg,[radius_m*ones(size(z)),zeros(size(z)),z],1);
    measurements.thin_loop_mT(rows) = 1e3*abs(field_T_per_A(:,3)).*measurements.peak_current_A(rows);
end
writetable(measurements,fullfile(here,'hall_depths.csv'));
scores = table();
transfer = table();
for coil = ["A","B"]
    series = measurements(measurements.coil==coil,:);
    currents = unique(series.peak_current_A);
    for current = currents'
        select = series.peak_current_A==current;
        truth = series.Hall_mT(select);
        pred = series.thin_loop_mT(select);
        author = series.author_COMSOL_mT(select);
        scores = [scores;table(coil,current,sum(select),100*norm(pred-truth)/norm(truth),100*norm(author-truth)/norm(truth), ...
            'VariableNames',{'coil','peak_current_A','depths','thin_loop_L2_percent','author_curve_L2_percent'})]; %#ok<AGROW>
    end
    low = series(series.peak_current_A==currents(1),:);
    high = series(series.peak_current_A==currents(2),:);
    [depth,lo,hi] = intersect(low.depth_mm,high.depth_mm);
    predicted = low.Hall_mT(lo)*currents(2)/currents(1);
    repeated = high.Hall_mT(hi);
    transfer = [transfer;table(repmat(coil,numel(depth),1),depth,low.Hall_mT(lo),predicted,repeated, ...
        'VariableNames',{'coil','depth_mm','low_drive_mT','predicted_high_drive_mT','measured_high_drive_mT'})]; %#ok<AGROW>
end
writetable(scores,fullfile(here,'thin_loop_scores.csv'));
writetable(transfer,fullfile(here,'current_transfer.csv'));
disp(scores)
for coil=["A","B"]
    t=transfer(transfer.coil==coil,:);
    fprintf('%s: held-out high-current L2 %.3f%% on %d depths\n',coil, ...
        100*norm(t.predicted_high_drive_mT-t.measured_high_drive_mT)/norm(t.measured_high_drive_mT),height(t));
end
end

function g=read_graph(here,k)
cache=fullfile(here,'..','..','tmp','abbasi_published_graphs');
if ~isfolder(cache), mkdir(cache); end
name=sprintf('sensors-24-01584-g%03d.png',6+k);
p=fullfile(cache,name);
if ~isfile(p)
    websave(p,['https://mdpi-res.com/sensors/sensors-24-01584/article_deploy/html/images/' name],weboptions('Timeout',45));
end
hashes=["3a3f5a73bc3b57cde6beda0608e43a55856802692936e99ac3195e6ad15db304", ...
    "61b8d0f4c7511a27d451af787f574f4823195b4f94ac3b7dab25967632c3ec37"];
fid=fopen(p,'rb'); cleanup=onCleanup(@() fclose(fid));
hash=java.security.MessageDigest.getInstance('SHA-256');
hash.update(fread(fid,Inf,'*uint8'));
assert(strcmpi(reshape(dec2hex(typecast(hash.digest(),'uint8'),2).',1,[]),hashes(k)));
[raw,map]=imread(p);
if isempty(map), g=raw(:,:,1); else, g=uint8(round(255*ind2rgb(raw,map))); g=g(:,:,1); end
end

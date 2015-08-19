clear all;
close all;

basin = 'pend';

basedir='/raid3/muxiao/bpa_inverse_routing/data/inverse_ro/input/';
outpdir='/raid3/muxiao/bpa_inverse_routing/data/inverse_ro/output/';

% basin grid georeferences
f_header = [basedir basin '.inputs/' basin '.header'];
fid = fopen(f_header, 'r');
while ~feof(fid)
   fline = fgets(fid);
   flins = textscan(fline, '%s');
   fcont = flins{1};
   comnd = [fcont{1} '=' fcont{2} ';'];
   eval(comnd);
end
fclose(fid);


tmpa_flag = false;
usgs_flag = true;

if (tmpa_flag)
    initname = 'tmpa';
else
    initname = 'null';
    tmpa_mean = 1.4745;
end

if (usgs_flag)
    strmname = 'usgs';
else
    strmname = 'nldas';
end

ssteps = 70;

load([outpdir basin '/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/all_data.mat']);

% load rainfall data
%
fin1 = fopen('runoff_precip/data_nldas.bin', 'r');

prec1_save = zeros(ncells, ssteps, nwins);

for s=1:ksteps
    prec1 = fread(fin1, [ncols nrows], 'float32');
    runoff1 = fread(fin1, [ncols nrows], 'float32');
end

for w=1:nwins
    for s=1:ssteps
        
        prec1 = fread(fin1, [ncols nrows], 'float32');
        runoff1 = fread(fin1, [ncols nrows], 'float32');
        
        prec1 = flipud(prec1').*basin_mask;
        prec1(prec1<0) = 0;
        prec1_compact = prec1(:);
        prec1_compact(isnan(prec1_compact)) = [];
        
        prec1_save(:, s, w) = prec1_compact;
        
    end
end

fclose(fin1);


gauge_lat = gauge_list(:, 3);
gauge_lon = gauge_list(:, 4) - 360;

gauge_i = nrows-floor((gauge_lat-yllcorner)/cellsize);
gauge_j = floor((gauge_lon-xllcorner)/cellsize)+1;

%% plots

for w=1:nwins
    
    runoff1_combine = squeeze(runoff1_save(:, w));
    runoff2_combine = squeeze(runoff2_save(:, w));
    runoff2_combine_post = squeeze(runoff2_save_post(:, w));
    
    for s_sel = 1:ssteps;

        runoff1_sel = runoff1_combine((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff1_sel_flat = runoff1_sel(cmap);
        runoff1_sel_2d = reshape(runoff1_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;

        runoff2_sel = runoff2_combine((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff2_sel_flat = runoff2_sel(cmap);
        runoff2_sel_2d = reshape(runoff2_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;

        runoff2_post_sel = runoff2_combine_post((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff2_post_sel_flat = runoff2_post_sel(cmap);
        runoff2_post_sel_2d = reshape(runoff2_post_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;
        
        diff_prio = runoff1_sel_2d(:) - runoff2_sel_2d(:);
        diff_post = runoff1_sel_2d(:) - runoff2_post_sel_2d(:);
        
        runoff_rms_prio(ssteps*(w-1)+s_sel) = sqrt(nanmean(diff_prio.*diff_prio));
        runoff_rms_post(ssteps*(w-1)+s_sel) = sqrt(nanmean(diff_post.*diff_post));

        runoff_bias_prio(ssteps*(w-1)+s_sel) = nanmean(-diff_prio);
        runoff_bias_post(ssteps*(w-1)+s_sel) = nanmean(-diff_post);
        
        if (w==1&&s_sel==59)
        close;
        figure('Position',[1 1 1200 800]);
        
        shx = 0.15;
        shy = 0.1;
        cmax = max([max(max(runoff2_sel_2d)) max(max(runoff1_sel_2d)) max(max(runoff2_post_sel_2d))]) * 0.3;
        % cmax = max([max(max(runoff2_sel_2d)) max(max(runoff1_sel_2d)) max(max(runoff2_post_sel_2d))]) * 0.6;
        cmax = min([cmax 10]);
        
        increment = runoff2_post_sel_2d-runoff2_sel_2d;

        %colormap(buildcmap('wygbr'));
        colormap(othercolor('PuBuGn9'));
        % colormap('jet');
        %colormap([[1 1 1]; othercolor('PuBuGn9')]);
        subplot(2, 2, 1);
        % runoff2_sel_2d(isnan(runoff2_sel_2d)) = -9999;
        h = imagescnan(runoff2_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff2_sel_2d));
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 cmax]);
        xlim([1 120]);
        if (tmpa_flag)
            text(10, 4, 'a) Initial Guess (TMPA-derived)', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Times');
        else
            text(10, 4, ['a) Initial Guess (Null, ' num2str(tmpa_mean, '%.3f') ' mm/day)'], 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Times');
        end
        freezeColors;
        %cbfreeze(colorbar('FontSize', 13));
    
        subplot(2, 2, 2);
        %runoff1_sel_2d(isnan(runoff1_sel_2d)) = -9999;
        h = imagescnan(runoff1_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff1_sel_2d));
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 cmax]);
        xlim([1 120]);
        if (~usgs_flag)
            text(10, 4, 'b) Synthetic Truth (NLDAS-derived)', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Times');
        else
            text(10, 4, 'b) NLDAS-derived Reference', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Times');
        end
        i_row = 1; i_col = 2;
        curpos=get(gca, 'Position');
        curpos(1)=curpos(1)-shx*(i_col-1);
        curpos(2)=curpos(2)+shy*(i_row-1);
        set(gca, 'Position', curpos);
        freezeColors;
        h = cbfreeze(colorbar('FontSize', 12, 'Location', 'East'));
        xlabel(h, 'mm/day');

        subplot(2, 2, 3);
        hold on;
        %runoff2_post_sel_2d(isnan(runoff2_post_sel_2d)) = -9999;
        h = imagescnan(runoff2_post_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff2_post_sel_2d));
        plot(gauge_j, gauge_i, '.', 'MarkerSize', 15, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        set(gca,'YDir','reverse');
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 cmax]);
        xlim([1 120]);
        text(10, 4, 'c) Inverted', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Times');
        i_row = 2; i_col = 1;
        curpos=get(gca, 'Position');
        curpos(1)=curpos(1)-shx*(i_col-1);
        curpos(2)=curpos(2)+shy*(i_row-1);
        set(gca, 'Position', curpos);
        freezeColors;
        %cbfreeze(colorbar('FontSize', 13));

        %colormap([[1 1 1]; othercolor('BrBu_12')]);
        colormap(othercolor('BrBu_12'));
        subplot(2, 2, 4);
        hold on;
        %increment(isnan(increment)) = -9999;
        h = imagescnan(increment);
        %set(h, 'AlphaData', ~isnan(runoff2_post_sel_2d));
        plot(gauge_j, gauge_i, '.', 'MarkerSize', 15, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        set(gca,'YDir','reverse');
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([-cmax cmax]);
        xlim([1 120]);
        text(10, 4, 'd) Inversion Increment', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Times');
        i_row = 2; i_col = 2;
        curpos=get(gca, 'Position');
        curpos(1)=curpos(1)-shx*(i_col-1);
        curpos(2)=curpos(2)+shy*(i_row-1);
        set(gca, 'Position', curpos);
        freezeColors;
        h = cbfreeze(colorbar('FontSize', 12, 'Location', 'East'));
        xlabel(h, 'mm/day');

        set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf, ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/inversion_day_' int2str(w) '_' int2str(s_sel) '_withgauges'], 'png');
        print(gcf, '-depsc2', '-painters', ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/inversion_day_' int2str(w) '_' int2str(s_sel) '_withgauges.eps']);
        end
    end

end

%% plot time series of basin mean bias and rmse

close;
figure('Position',[1 1 1200 900]);

fs=15;

subplot(3, 1, 1);
hold on;
plot(ksteps+1:ksteps+ssteps*nwins, runoff_bias_prio(1:ssteps*nwins), '-', 'Color', [0 0 1], 'LineWidth', 1.5);
plot(ksteps+1:ksteps+ssteps*nwins, runoff_bias_post(1:ssteps*nwins), '-', 'Color', [1 0 0], 'LineWidth', 1.5);
plot([ksteps+1 nsteps], [0 0], '-k');

set(gca, 'FontSize', fs);

avebias_prio = mean(abs(runoff_bias_prio(1:ssteps*nwins)));
avebias_post = mean(abs(runoff_bias_post(1:ssteps*nwins)));

xlim([ksteps+1 ksteps+ssteps*nwins]);
yr = max(abs([runoff_bias_prio(1:ssteps*nwins); runoff_bias_post(1:ssteps*nwins)]));
ylim([-yr yr]);

%xlabel('Day in 2009', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bias (mm/day)', 'FontSize', fs, 'FontWeight', 'bold');
if (tmpa_flag)
    leg = legend('Initial Guess (TMPA-derived) $$r$$', 'Inverted $$r$$', 'Location', 'SouthEast');
else
    leg = legend('Initial Guess (Null) $$r$$', 'Inverted $$r$$', 'Location', 'SouthEast');
end
set(leg, 'FontSize', fs, 'FontWeight', 'normal', 'FontName', 'Times', 'Interpreter', 'latex');
text(ksteps+5, yr*0.9, 'a) Basin Mean Bias', 'FontSize', fs, 'FontWeight', 'normal', 'FontName', 'Times');
text(ksteps+11, yr*0.65, ['Time Average (Absolute Bias) = ' num2str(avebias_prio, '%.3f') ' mm/day (Initial Guess) and ' num2str(avebias_post, '%.3f') ' mm/day (Inverted)'], ...
    'FontSize', fs, 'FontWeight', 'normal', 'FontName', 'Times');

subplot(3, 1, 2);
hold on;
plot(ksteps+1:ksteps+ssteps*nwins, runoff_rms_prio(1:ssteps*nwins), '-', 'Color', [0 0 1], 'LineWidth', 1.5);
plot(ksteps+1:ksteps+ssteps*nwins, runoff_rms_post(1:ssteps*nwins), '-', 'Color', [1 0 0], 'LineWidth', 1.5);

set(gca, 'FontSize', fs);

averms_prio = mean(runoff_rms_prio(1:ssteps*nwins));
averms_post = mean(runoff_rms_post(1:ssteps*nwins));

xlim([ksteps+1 ksteps+ssteps*nwins]);
xlabel('Day in 2009', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('RMSE (mm/day)', 'FontSize', fs, 'FontWeight', 'bold');
yyll=ylim;
text(ksteps+5, yyll(2)*0.9, 'b) Basin Root Mean Squared Errors', 'FontSize', fs, 'FontWeight', 'normal', 'FontName', 'Times');
text(ksteps+11, yyll(2)*0.77, ['Time Average = ' num2str(averms_prio, '%.3f') ' mm/day (Initial Guess) and ' num2str(averms_post, '%.3f') ' mm/day (Inverted)'], ...
    'FontSize', fs, 'FontWeight', 'normal', 'FontName', 'Times');

i_row = 2;
curpos=get(gca, 'Position');
% curpos(1)=curpos(1)-0.08*(i_col-1);
curpos(2)=curpos(2)+0.05*(i_row-1);
set(gca, 'Position', curpos);

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/runoff_rms_bias'], 'png');
print(gcf, '-depsc2', '-painters', ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/runoff_rms_bias.eps']);

%% plot some hydrographs
close;
figure('Position',[1 1 1200 900]);
letters='abcdefghijklmnopqrstuvwxyz';

tstep    = 86400;    % timestep size    [s]

fs = 15;

%g_sels = [80 40 20 60];
if (~usgs_flag)
    g_sels = [28 12 40 60];
else
    g_sels = [28 15 42 61];
end
g_sels = [1 2 3 4];

i_row = 0;
%for g_sel=20:20:ngauges
for g_sel=g_sels
    
    i_row = i_row+1;
    subplot(length(g_sels), 1, i_row);
    hold on;
    if (~usgs_flag)
        plot(ksteps+1:nsteps, streamflow_gauge1_all(g_sel, :)/tstep, '-', 'Color', [0 1 0], 'LineWidth', 4);
        plot(ksteps+1:nsteps, streamflow_gauge2_all(g_sel, :)/tstep, '-', 'Color', [0 0 1], 'LineWidth', 1.5);
        plot(ksteps+1:nsteps, streamflow_gauge2_post_all(g_sel, :)/tstep, '-', 'Color', [1 0 0], 'LineWidth', 1.5);
    else
        plot(ksteps+1:nsteps, streamflow_usgs(g_sel, ksteps+1:nsteps)/tstep, '-', 'Color', [0 1 0], 'LineWidth', 4);
        plot(ksteps+1:nsteps, streamflow_gauge2_all(g_sel, :)/tstep, '-', 'Color', [0 0 1], 'LineWidth', 1.5);
        plot(ksteps+1:nsteps, streamflow_gauge2_post_all(g_sel, :)/tstep, '-', 'Color', [1 0 0], 'LineWidth', 1.5);
    end
    set(gca, 'FontSize', fs);
    
    ts=sprintf('USGS 0%d, %.0f km^2', gauge_list(g_sel, 2), gauge_area(g_sel)/1e6);
    yyll=ylim;
    text(20, yyll(2)*0.9, [letters(i_row) ') ' ts], 'FontSize', fs, 'FontWeight', 'normal', 'FontName', 'Times');
    xlim([ksteps+1 ksteps+ssteps*nwins]);
    if (i_row==2)
        ylabel('Streamflow (m^3/s)', 'FontSize', fs, 'FontWeight', 'bold');
    end
    if (i_row==length(g_sels)-3)
        if (tmpa_flag)
            leg = legend('Synthetic Truth $$Q$$ (NLDAS-derived)', 'Initial Guess $$Q$$ (TMPA-derived)', '$$Q$$ Reconstructed from Inverted $$r$$', 'Location', 'North', 'FillColor', 'none');
        else
            leg = legend('Synthetic Truth $$Q$$ (NLDAS-derived)', 'Initial Guess $$Q$$ (Null)', '$$Q$$ Reconstructed from Inverted $$r$$', 'Location', 'NorthEast');
        end
        set(leg, 'FontSize', fs-2, 'FontWeight', 'normal', 'FontName', 'Times', 'Interpreter', 'latex');
        lpos = get(leg, 'Position');
        lpos(1) = lpos(1)*1.05;
        set(leg, 'Position', lpos);
    end
    if (i_row==length(g_sels))
        xlabel('Day in 2009', 'FontSize', fs, 'FontWeight', 'bold');
    end
    
    curpos=get(gca, 'Position');
    % curpos(1)=curpos(1)-0.08*(i_col-1);
    curpos(2)=curpos(2)+0.025*(i_row-1);
    set(gca, 'Position', curpos);

end

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/streamflow_truth_vs_guess_post'], 'png');
print(gcf, '-depsc2', '-painters', ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/streamflow_truth_vs_guess_post.eps']);

if (usgs_flag)
    
    close;
    figure('Position',[1 1 1200 900]);

    i_row = 0;
    for g_sel=g_sels
    
        i_row = i_row+1;
        subplot(length(g_sels), 1, i_row);
        hold on;
        plot(ksteps+1:nsteps, streamflow_gauge1_all(g_sel, :)/tstep, 'k.', 'MarkerSize', 10);
        plot(ksteps+1:nsteps, streamflow_usgs(g_sel, ksteps+1:nsteps)/tstep, '-', 'Color', [0 1 0], 'LineWidth', 1.5);
        set(gca, 'FontSize', fs);
        
        ts=sprintf('USGS 0%d, %.0f km^2', gauge_list(g_sel, 2), gauge_area(g_sel)/1e6);
        yyll=ylim;
        text(20, yyll(2)*0.9, [letters(i_row) ') ' ts], 'FontSize', fs, 'FontWeight', 'normal', 'FontName', 'Times');
        xlim([ksteps+1 ksteps+ssteps*nwins]);
        if (i_row==2)
            ylabel('Streamflow (m^3/s)', 'FontSize', fs, 'FontWeight', 'bold');
        end
        if (i_row==length(g_sels)-3)
            leg = legend('USGS Measurements', 'Synthetic Truth $$Q$$ (NLDAS-derived)', 'Location', 'North');
            set(leg, 'FontSize', fs-1, 'FontWeight', 'normal', 'FontName', 'Times', 'Interpreter', 'latex');
            lpos = get(leg, 'Position');
            lpos(1) = lpos(1)*1.1;
            set(leg, 'Position', lpos);
        end
        if (i_row==length(g_sels))
            xlabel('Day in 2009', 'FontSize', fs, 'FontWeight', 'bold');
        end
    
        curpos=get(gca, 'Position');
        % curpos(1)=curpos(1)-0.08*(i_col-1);
        curpos(2)=curpos(2)+0.025*(i_row-1);
        set(gca, 'Position', curpos);

    end

    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/streamflow_usgs_vs_truth'], 'png');
    print(gcf, '-depsc2', '-painters', ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/streamflow_usgs_vs_truth.eps']);

end

%%

close;
figure('Position',[1 1 1600 500]);

wsels = [1 2 3 5];
nwsels = length(wsels);

i_col = 0;
for w=wsels
    
    runoff1_combine = squeeze(runoff1_save(:, w));
    runoff2_combine = squeeze(runoff2_save(:, w));
    runoff2_combine_post = squeeze(runoff2_save_post(:, w));
    
    s_sel = 21;

        runoff1_sel = runoff1_combine((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff1_sel_flat = runoff1_sel(cmap);
        runoff1_sel_2d = reshape(runoff1_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;

        runoff2_sel = runoff2_combine((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff2_sel_flat = runoff2_sel(cmap);
        runoff2_sel_2d = reshape(runoff2_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;

        runoff2_post_sel = runoff2_combine_post((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff2_post_sel_flat = runoff2_post_sel(cmap);
        runoff2_post_sel_2d = reshape(runoff2_post_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;
        
        shx = 0.075;
        shy = 0.1;
        cmax = max([max(max(runoff2_sel_2d)) max(max(runoff1_sel_2d)) max(max(runoff2_post_sel_2d))]) * 0.3;
        cmax = min([cmax 10]);
        cmax = 3;
        
        colormap(othercolor('PuBuGn9'));
        
        i_col = i_col + 1;
        
        i_row = 1; 
        pn = i_col+(i_row-1)*nwsels;
        
        subplot(2, nwsels, pn);
        %runoff1_sel_2d(isnan(runoff1_sel_2d)) = -9999;
        h = imagescnan(runoff1_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff1_sel_2d));
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 cmax]);
        xlim([1 120]);
        if (w==1)
            ylabel('Synthetic Truth', 'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Times', 'Color', 'black');
        end
        title(['Day ' num2str(ksteps+(w-1)*ssteps+s_sel, '%.0f')], 'FontSize', 14);
        curpos=get(gca, 'Position');
        curpos(1)=curpos(1)-shx*(i_col-1);
        curpos(2)=curpos(2)+shy*(i_row-1);
        set(gca, 'Position', curpos);
        freezeColors;
        if (i_col==nwsels)
            h = cbfreeze(colorbar('FontSize', 13, 'Location', 'East'));
            xlabel(h, 'mm/day');
        end

        i_row = 2; 
        pn = i_col+(i_row-1)*nwsels;
        
        subplot(2, nwsels, pn);
        %runoff2_post_sel_2d(isnan(runoff2_post_sel_2d)) = -9999;
        h = imagescnan(runoff2_post_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff2_post_sel_2d));
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 cmax]);
        xlim([1 120]);
        if (i_col==1)
            ylabel('Inverted', 'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Times', 'Color', 'black');
        end
        curpos=get(gca, 'Position');
        curpos(1)=curpos(1)-shx*(i_col-1);
        curpos(2)=curpos(2)+shy*(i_row-1);
        set(gca, 'Position', curpos);
        freezeColors;
        %cbfreeze(colorbar('FontSize', 13));

end

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/inversion_allwins_' int2str(s_sel)], 'png');
print(gcf, '-depsc2', '-painters', ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/inversion_allwins_' int2str(s_sel) '.eps']);


return;

%%

close;
figure('Position',[1 1 1600 900]);

wsels = [1 2 3 5];
accd = [10 7 7 7 7];
nwsels = length(wsels);

i_col = 0;
for w=wsels
    
    runoff1_combine = squeeze(runoff1_save(:, w));
    runoff2_combine = squeeze(runoff2_save(:, w));
    runoff2_combine_post = squeeze(runoff2_save_post(:, w));
    
    s_sel = 21;

        runoff1_sel = runoff1_combine((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff1_sel_flat = runoff1_sel(cmap);
        runoff1_sel_2d = reshape(runoff1_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;

        runoff2_sel = runoff2_combine((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff2_sel_flat = runoff2_sel(cmap);
        runoff2_sel_2d = reshape(runoff2_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;

        runoff2_post_sel = runoff2_combine_post((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff2_post_sel_flat = runoff2_post_sel(cmap);
        runoff2_post_sel_2d = reshape(runoff2_post_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;
        
        prec1_sel = squeeze(mean(prec1_save(:, s_sel-accd(w):s_sel, w), 2));
        prec1_sel_flat = prec1_sel(cmap);
        prec1_sel_2d = reshape(prec1_sel_flat, nrows, ncols) .* basin_mask;
        
        shx = 0.075;
        shy = 0.07;
        cmax = max([max(max(runoff2_sel_2d)) max(max(runoff1_sel_2d)) max(max(runoff2_post_sel_2d))]) * 0.3;
        cmax = min([cmax 10]);
        cmax = 3;
        
        colormap(othercolor('PuBuGn9'));
        
        i_col = i_col + 1;
        
        i_row = 1; 
        pn = i_col+(i_row-1)*nwsels;
        
        subplot(3, nwsels, pn);
        %runoff1_sel_2d(isnan(runoff1_sel_2d)) = -9999;
        h = imagescnan(runoff1_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff1_sel_2d));
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 cmax]);
        xlim([1 120]);
        if (w==1)
            ylabel('Synthetic Truth', 'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Times', 'Color', 'black');
        end
        title(['Day ' num2str(ksteps+(w-1)*ssteps+s_sel, '%.0f')], 'FontSize', 14);
        curpos=get(gca, 'Position');
        curpos(1)=curpos(1)-shx*(i_col-1);
        curpos(2)=curpos(2)+shy*(i_row-1);
        set(gca, 'Position', curpos);
        freezeColors;
        if (i_col==nwsels)
            h = cbfreeze(colorbar('FontSize', 13, 'Location', 'East'));
            xlabel(h, 'mm/day');
        end

        i_row = 2; 
        pn = i_col+(i_row-1)*nwsels;
        
        subplot(3, nwsels, pn);
        %runoff2_post_sel_2d(isnan(runoff2_post_sel_2d)) = -9999;
        h = imagescnan(runoff2_post_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff2_post_sel_2d));
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 cmax]);
        xlim([1 120]);
        if (i_col==1)
            ylabel('Inverted', 'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Times', 'Color', 'black');
        end
        curpos=get(gca, 'Position');
        curpos(1)=curpos(1)-shx*(i_col-1);
        curpos(2)=curpos(2)+shy*(i_row-1);
        set(gca, 'Position', curpos);
        freezeColors;
        %cbfreeze(colorbar('FontSize', 13));

        colormap(othercolor('RdPu9'));
        
        i_row = 3;
        pn = i_col+(i_row-1)*nwsels;
        
        subplot(3, nwsels, pn);
        %runoff2_post_sel_2d(isnan(runoff2_post_sel_2d)) = -9999;
        h = imagescnan(prec1_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff2_post_sel_2d));
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 9]);
        xlim([1 120]);
        if (i_col==1)
            ylabel('Precipitation', 'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Times', 'Color', 'black');
        end
        curpos=get(gca, 'Position');
        curpos(1)=curpos(1)-shx*(i_col-1);
        curpos(2)=curpos(2)+shy*(i_row-1);
        set(gca, 'Position', curpos);
        freezeColors;
        %cbfreeze(colorbar('FontSize', 13));
        
        if (i_col==nwsels)
            h = cbfreeze(colorbar('FontSize', 13, 'Location', 'East'));
            xlabel(h, 'mm/day');
        end

end

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/inversion_allwins_prec_' int2str(s_sel)], 'png');
print(gcf, '-depsc2', '-painters', ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/inversion_allwins_prec_' int2str(s_sel) '.eps']);


return;

%% examine some of the H's

% g_sel = 20;
% s_sel = 3;
% 
% H_sel = HH(g_sel, :, s_sel);
% 
% H_sel_flat = H_sel(cmap);
% H_sel_2d = reshape(H_sel_flat, nrows, ncols) .* basin_mask;
% 
% imagesc(H_sel_2d);

%% plot some hydrographs
close;

tstep    = 86400;    % timestep size    [s]

%g_sels = [80 40 20 60];
if (~usgs_flag)
    g_sels = [28 12 40 60];
else
    g_sels = [28 15 42 61];
end

i_row = 0;
%for g_sel=20:20:ngauges
for g_sel=g_sels
    
    i_row = i_row+1;
    close;
    figure('Position',[1 1 1200 300]);
    hold on;
    plot(ksteps+1:nsteps, streamflow_gauge1_all(g_sel, :)/tstep, '-', 'Color', [1 0 0], 'LineWidth', 2);
    yyll=ylim;
    xlim([ksteps+1 ksteps+ssteps*nwins]);
    
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/streamflow_truth_' int2str(i_row)], 'png');
    print(gcf, '-depsc2', '-painters', ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/streamflow_truth_' int2str(i_row) '.eps']);

end

%%

close all;

wsels = [1 2 3 5];
nwsels = length(wsels);

i_col = 0;
for w=wsels
    
    figure('Position',[1 1 800 600]);

    runoff1_combine = squeeze(runoff1_save(:, w));
    runoff2_combine = squeeze(runoff2_save(:, w));
    runoff2_combine_post = squeeze(runoff2_save_post(:, w));
    
    s_sel = 21;

        runoff1_sel = runoff1_combine((ssteps-s_sel)*ncells+1:(ssteps-s_sel+1)*ncells);
        runoff1_sel_flat = runoff1_sel(cmap);
        runoff1_sel_2d = reshape(runoff1_sel_flat, nrows, ncols) .* basin_mask ./grid_area *1000;
        
        shx = 0.075;
        shy = 0.1;
        cmax = max(max(runoff1_sel_2d)) * 0.3;
        cmax = min([cmax 10]);
        cmax = 3;
        
        colormap(othercolor('PuBuGn9'));
        
        i_col = i_col + 1;
        
        i_row = 1; 
        pn = i_col+(i_row-1)*nwsels;
        
        %runoff1_sel_2d(isnan(runoff1_sel_2d)) = -9999;
        h = imagescnan(runoff1_sel_2d);
        %set(h, 'AlphaData', ~isnan(runoff1_sel_2d));
        set(gca,'xcolor',get(gcf,'color'));
        set(gca,'ycolor',get(gcf,'color'));
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        caxis([0 cmax]);

    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/inversion_win_' int2str(i_col)], 'png');
    print(gcf, '-depsc2', '-painters', ['imgs/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/inversion_win_' int2str(i_col) '.eps']);

end



%% post save data into ascii

data_all=runoff2_save_post(:);
data_2_w=ones(ncells,length(data_all)/ncells)*-1;
fill_area=grid_area(:);
fill_area(isnan(fill_area))=[];
save_steps=size(runoff2_save_post,1)/ncells*nwins
for i=1:save_steps
    st=1+(i-1)*ncells;
    ed=ncells+(i-1)*ncells;
    data_2_w(:,i)=data_all(st:ed)./fill_area*1000;
end

% fliplr the output
for i=1:nwins
    st=1+(i-1)*ssteps;
    ed=i*ssteps;
    data_2_w(:,st:ed)=fliplr(data_2_w(:,st:ed));
end


[bm,bn]=size(basin_mask);
[m,n]=find(~isnan(basin_mask));
lon=0.5*cellsize+xllcorner+(n-1)*cellsize;
lat=0.5*cellsize+yllcorner+(bn-1)*cellsize-(m-1)*cellsize;
lon_lat_data=[lon,lat,data_2_w];

dlmwrite(['./results/',basin,'/data_all_day_',num2str(ssteps)],lon_lat_data,...
    'delimiter', '\t', 'newline', 'unix','precision', 9);


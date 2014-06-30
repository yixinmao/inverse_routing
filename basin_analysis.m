% do all the basin analysis

function [ basin_mask grid_lat grid_lon flow_length grid_area flow_accum ...
           ds_i ds_j dmax ds_length gauge_lat gauge_lon gauge_i gauge_j ...
           gauge_area gauge_accum basin_mask_gauge gauge_list flow_length_gauge ] ...
           = basin_analysis(flow_dir, georef, gauge_list, flow_vel, tstep, flag_plot)
      
r_earth = 6371.0072 * 1000;  % earth radius [m]

% basin grid georeferences
xllcorner    = georef(1);
yllcorner    = georef(2);
cellsize     = georef(3);
nodata_value = georef(4);

[nrows ncols] = size(flow_dir);
flow_dir(flow_dir==nodata_value) = NaN;

%% routing basin analysis

% basin mask
basin_mask = flow_dir*0+1;
ncells = nansum(nansum(basin_mask));

% lat and lon grids
grid_lat = repmat(flipud((yllcorner+((1:nrows)-0.5)*cellsize)'), 1, ncols).*basin_mask;
grid_lon = repmat(xllcorner+((1:ncols)-0.5)*cellsize, nrows, 1).*basin_mask;

% grid area
grid_area = cellsize/180*pi*r_earth*cellsize/180*pi*r_earth*cos(grid_lat/180*pi);

% lat and lon of 1st order downstream, for stream plotting purpose
ds1_lat = grid_lat;
ds1_lat(flow_dir==1|flow_dir==2|flow_dir==8) = ds1_lat(flow_dir==1|flow_dir==2|flow_dir==8) + cellsize;
ds1_lat(flow_dir==4|flow_dir==5|flow_dir==6) = ds1_lat(flow_dir==4|flow_dir==5|flow_dir==6) - cellsize;

ds1_lon = grid_lon;
ds1_lon(flow_dir==2|flow_dir==3|flow_dir==4) = ds1_lon(flow_dir==2|flow_dir==3|flow_dir==4) + cellsize;
ds1_lon(flow_dir==6|flow_dir==7|flow_dir==8) = ds1_lon(flow_dir==6|flow_dir==7|flow_dir==8) - cellsize;

lats = [reshape(grid_lat(basin_mask==1), ncells, 1) reshape(ds1_lat(basin_mask==1), ncells, 1)];
lons = [reshape(grid_lon(basin_mask==1), ncells, 1) reshape(ds1_lon(basin_mask==1), ncells, 1)];

% downstream flow distance, great circle distance using haversine formula
hslat = (1-cos((grid_lat-ds1_lat)/180*pi))/2;
hslon = (1-cos((grid_lon-ds1_lon)/180*pi))/2;
flow_dist = 2*r_earth*asin(sqrt(hslat + cos(grid_lat/180*pi).*cos(ds1_lat/180*pi).*hslon));

% downstream grid search at all orders
search_max = 200;
ds_i = zeros(nrows, ncols, search_max);
ds_j = zeros(nrows, ncols, search_max);
ds_length = zeros(nrows, ncols, search_max);

% initialize
ds_i(:, :, 1) = repmat((1:nrows)', 1, ncols);
ds_j(:, :, 1) = repmat(1:ncols, nrows, 1);
ds_flow_dir = flow_dir;
flow_dir_flat = flow_dir(:);
basin_mask_flat = basin_mask(:);

% flow length to outlet
flow_length = basin_mask*0;
flow_length_flat = flow_length(:);
flow_dist_flat = flow_dist(:);

% contributing area
grid_area_flat = grid_area(:);
flow_accum = grid_area;
flow_accum_flat = flow_accum(:);

for d=2:search_max
    
    ds_i_tmp = squeeze(ds_i(:, :, d-1));
    ds_j_tmp = squeeze(ds_j(:, :, d-1));
    
    ds_i_flat = ds_i_tmp(:);
    ds_j_flat = ds_j_tmp(:);
    
    % accumulate flow length
    flow_length_flat = flow_length_flat + flow_dist_flat(ds_i_flat+(ds_j_flat-1)*nrows);
    ds_length(:, :, d) = reshape(flow_length_flat, nrows, ncols);
    
    % move downstream
    ds_i_tmp(ds_flow_dir==1|ds_flow_dir==2|ds_flow_dir==8) = ds_i_tmp(ds_flow_dir==1|ds_flow_dir==2|ds_flow_dir==8) - 1;
    ds_i_tmp(ds_flow_dir==4|ds_flow_dir==5|ds_flow_dir==6) = ds_i_tmp(ds_flow_dir==4|ds_flow_dir==5|ds_flow_dir==6) + 1;
    
    ds_j_tmp(ds_flow_dir==2|ds_flow_dir==3|ds_flow_dir==4) = ds_j_tmp(ds_flow_dir==2|ds_flow_dir==3|ds_flow_dir==4) + 1;
    ds_j_tmp(ds_flow_dir==6|ds_flow_dir==7|ds_flow_dir==8) = ds_j_tmp(ds_flow_dir==6|ds_flow_dir==7|ds_flow_dir==8) - 1;
    
    ds_i(:, :, d) = ds_i_tmp;
    ds_j(:, :, d) = ds_j_tmp;
    
    ds_i_tmp_flat = ds_i_tmp(:);
    ds_j_tmp_flat = ds_j_tmp(:);
    
    % key step -- find the flow direction of the downstream cell
    ds_flow_dir_flat = flow_dir_flat(ds_i_tmp_flat+(ds_j_tmp_flat-1)*nrows);
    ds_flow_dir = reshape(ds_flow_dir_flat, nrows, ncols);
        
    % accumulate contributing area
    % too bad the following doesn't work:
    % flow_accum_flat(ds_i_flat+(ds_j_flat-1)*nrows) = flow_accum_flat(ds_i_flat+(ds_j_flat-1)*nrows) + grid_area_flat;
    for i=1:numel(ds_i_flat)
        if ( ~isnan(ds_i_flat(i)) && (ds_i_flat(i)~=ds_i_tmp_flat(i) || ds_j_flat(i)~=ds_j_tmp_flat(i)) )
            flow_accum_flat(ds_i_flat(i)+(ds_j_flat(i)-1)*nrows) = flow_accum_flat(ds_i_flat(i)+(ds_j_flat(i)-1)*nrows) + grid_area_flat(i);
        end
    end

    % all reach the outlet?
    if (nanstd(ds_i_tmp_flat.*basin_mask_flat)==0 && nanstd(ds_j_tmp_flat.*basin_mask_flat)==0)
        dmax = d;
        outlet_i = nanmean(ds_i_tmp_flat.*basin_mask_flat);
        outlet_j = nanmean(ds_j_tmp_flat.*basin_mask_flat);
        break;
    end
    
end

flow_length = reshape(flow_length_flat, nrows, ncols);
flow_accum = reshape(flow_accum_flat, nrows, ncols);

% delete extra downstream map
ds_i(:, :, dmax+1:search_max) = [];
ds_j(:, :, dmax+1:search_max) = [];
ds_length(:, :, dmax+1:search_max) = [];

%% extract sub-basins for gauge stations

gauge_lat = gauge_list(:, 3);
gauge_lon = gauge_list(:, 4) - 360;

gauge_i = nrows-floor((gauge_lat-yllcorner)/cellsize);
gauge_j = floor((gauge_lon-xllcorner)/cellsize)+1;

gauge_accum = flow_accum(gauge_i+(gauge_j-1)*nrows);
gauge_area = gauge_list(:, 5)*1609*1609;
%plot(gauge_area, gauge_accum, '.');

% remove stations with too inaccurate flow area
bad_gauges = gauge_area./gauge_accum>2 | gauge_area./gauge_accum<0.5;

% gauge_list(bad_gauges, :) = [];
% gauge_lat(bad_gauges) = [];
% gauge_lon(bad_gauges) = [];
% 
% gauge_i(bad_gauges) = [];
% gauge_j(bad_gauges) = [];
% 
% gauge_accum(bad_gauges) = [];
% gauge_area(bad_gauges) = [];

ngauges = numel(gauge_lat);

% sub-basin flow length

basin_mask_gauge = zeros(nrows, ncols, ngauges);
flow_length_gauge = zeros(nrows, ncols, ngauges);

for g=1:ngauges
    basin_mask_gauge(:, :, g) = basin_mask*0;
end

for i=1:nrows
    for j=1:ncols
        
        if (isnan(basin_mask(i, j)))
            continue;
        end
        
        ds_i_cell = squeeze(ds_i(i, j, :));
        ds_j_cell = squeeze(ds_j(i, j, :));
        
        for g=1:ngauges
            if (~isempty(find(ds_i_cell==gauge_i(g)&ds_j_cell==gauge_j(g),1)))
                basin_mask_gauge(i, j, g) = 1;
            end
        end
        
    end
end

for g=1:ngauges
    flow_length_gauge(:, :, g) = (flow_length - flow_length(gauge_i(g), gauge_j(g))) .* squeeze(basin_mask_gauge(:, :, g));
end

%% make some plots

if (flag_plot~=1)
    return
end

close all;

figure('Position',[1 1 1200 900]);

% flow length background
colormap('cool');
h = imagesc(flow_length/flow_vel/tstep);
set(h, 'AlphaData', ~isnan(basin_mask));
set(gca,'xcolor',get(gcf,'color'));
set(gca,'ycolor',get(gcf,'color'));
set(gca,'ytick',[]);
set(gca,'xtick',[]);
colorbar('Location', 'east', 'FontSize', 15);
text(106, 45,'Travel Time (days)', 'Rotation', 90, 'FontSize', 20, 'FontWeight', 'bold');

% plot the stream network

hold on;
stream_width = flow_accum_flat(~isnan(flow_accum_flat));
stream_width = sqrt(sqrt(stream_width/max(stream_width))) * 7;

for c = 1:ncells
    % line(lons(c,:), lats(c,:), 'LineWidth', 1, 'Color', 'black');
    line((lons(c,:)-xllcorner)/cellsize+0.5, nrows-(lats(c,:)-yllcorner)/cellsize+0.5, 'LineWidth', stream_width(c), 'Color', 'white');
end

plot3(outlet_j, outlet_i, 0.1, 'pentagram', 'MarkerSize', 15, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, 'imgs/flow_network', 'png');

% plot a sample flow path
%colorbar('off');
path_i = squeeze(ds_i(10, 98, :));
path_j = squeeze(ds_j(10, 98, :));

% line(xllcorner+(path_j-0.5)*cellsize, yllcorner+(nrows-path_i+0.5)*cellsize, 'Color', 'black', 'LineWidth', 1);
line(path_j, path_i, 'Color', 'black', 'LineWidth', 2.5);

plot3(98, 10, 0.1, '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black');

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, 'imgs/sample_path', 'png');

hold off;

%% make a gauge plot
figure('Position',[1 1 1200 900]);

% background
colormap([0.7 0.7 0.7]);
h = imagesc(flow_length*0);
set(h, 'AlphaData', ~isnan(basin_mask));
set(gca,'xcolor',get(gcf,'color'));
set(gca,'ycolor',get(gcf,'color'));
set(gca,'ytick',[]);
set(gca,'xtick',[]);

hold on;
for c = 1:ncells
    % line(lons(c,:), lats(c,:), 'LineWidth', 1, 'Color', 'black');
    line((lons(c,:)-xllcorner)/cellsize+0.5, nrows-(lats(c,:)-yllcorner)/cellsize+0.5, 'LineWidth', stream_width(c), 'Color', 'white');
end

plot3(gauge_j, gauge_i, gauge_j*0+0.1, '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, 'imgs/gauge_map', 'png');

%% make a sub-basin map

%for g_sel = 1:ngauges
for g_sel = 1:1
    
    basin_mask_sel = squeeze(basin_mask_gauge(:, :, g_sel));
    flow_length_sel = squeeze(flow_length_gauge(:, :, g_sel));

    close;
    figure('Position',[1 1 1200 900]);

    %colormap([0.7 0.7 0.7; 0 0 0]);
    colormap([[0.7 0.7 0.7]; colormap('cool')]);

    h = imagesc(flow_length_sel);
    set(h, 'AlphaData', ~isnan(basin_mask_sel));
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);

    hold on;
    
    for c = 1:ncells
        % line(lons(c,:), lats(c,:), 'LineWidth', 1, 'Color', 'black');
        line((lons(c,:)-xllcorner)/cellsize+0.5, nrows-(lats(c,:)-yllcorner)/cellsize+0.5, 'LineWidth', stream_width(c), 'Color', 'white');
    end

    plot3(gauge_j(g_sel), gauge_i(g_sel), 0.1, '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');

    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, ['imgs/gauge_basin_0' int2str(gauge_list(g_sel,2))], 'png');
    
end

%% make travel time map for one sub basin

g_sel = 37;

basin_mask_sel = squeeze(basin_mask_gauge(:, :, g_sel));
flow_length_sel = squeeze(flow_length_gauge(:, :, g_sel));
flow_days_sel = flow_length_sel/flow_vel/tstep;
flow_days_sel(flow_days_sel==0) = -1;

max_days = ceil(max(max(flow_days_sel)));

for d=0:max_days
    
    close;
    figure('Position',[1 1 1200 900]);

    if (d==0)
        colormap([[0.7 0.7 0.7]; colormap('cool')]);
        h = imagesc(flow_days_sel);
        colorbar('Location', 'east', 'YLim', [0 max(max(flow_days_sel))], 'FontSize', 15);
        text(106, 45,'Travel Time (days)', 'Rotation', 90, 'FontSize', 20, 'FontWeight', 'bold');
    else
        colordat = colormap('cool');
        colormap([[0.7 0.7 0.7]; repmat([0.5 0.5 0.5], length(colordat)+5, 1); colordat]);
        flow_days_tmp = flow_days_sel;
        flow_days_tmp(flow_days_tmp<d-1&flow_days_tmp>=0) = 0;
        flow_days_tmp(flow_days_tmp>d) = 0;
        h = imagesc(flow_days_tmp);
    end
        
    set(h, 'AlphaData', ~isnan(basin_mask_sel));
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    
    hold on;
    
    for c = 1:ncells
        % line(lons(c,:), lats(c,:), 'LineWidth', 1, 'Color', 'black');
        if (basin_mask_sel(nrows-(lats(c,2)-yllcorner)/cellsize+0.5, (lons(c,2)-xllcorner)/cellsize+0.5) == 1)
            line((lons(c,:)-xllcorner)/cellsize+0.5, nrows-(lats(c,:)-yllcorner)/cellsize+0.5, 'LineWidth', stream_width(c), 'Color', 'white');
        end
    end

    plot3(gauge_j(g_sel), gauge_i(g_sel), 0.1, '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
    
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, ['imgs/gauge_basin_0' int2str(gauge_list(g_sel,2)) '_' int2str(d)], 'png');

end

%% make gauge basin selection maps

g_sels1 = [37 36 35 26 25 24 20 16 12 2];
g_sels2 = [46 57 62 76 85 96 ];
g_sels3 = [112 106 101];
g_sels = [g_sels1 g_sels2 g_sels3];

for g=1:length(g_sels)
    
    mask_tmp = squeeze(basin_mask_gauge(:, :, g_sels(g)));
    mask_tmp(mask_tmp==0) = NaN;
    
    close;
    figure('Position',[1 1 1200 900]);
    
    colordat = colormap('hsv');
    
    colormap(colordat(randi(length(colordat)),:));
    
    h = imagesc(mask_tmp);
    set(h, 'AlphaData', ~isnan(mask_tmp));
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    
    hold on;
    for c = 1:ncells
        % line(lons(c,:), lats(c,:), 'LineWidth', 1, 'Color', 'black');
        if (mask_tmp(nrows-(lats(c,2)-yllcorner)/cellsize+0.5, (lons(c,2)-xllcorner)/cellsize+0.5) == 1)
            line((lons(c,:)-xllcorner)/cellsize+0.5, nrows-(lats(c,:)-yllcorner)/cellsize+0.5, 'LineWidth', stream_width(c), 'Color', 'white');
        end
    end
    plot3(gauge_j(g_sels(g)), gauge_i(g_sels(g)), gauge_j(g_sels(g))*0+0.1, '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');

    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, ['imgs/gauge_basins_white/gauge_sel_' int2str(g)], 'png');

end

%% plot some sample Green's functions

x = 1:0.5:100;
C = 1.5;
D = 0.5;

close;
figure('Position',[1 1 1200 900]);

for t=10:10:50
    
    subplot(5, 1, t/10);
    
    plot(x, x/2/t/sqrt(t*D).*exp(-(C*t-x).*(C*t-x)/4/D/t), '-', 'LineWidth', 3);
    
    %set(gca,'xcolor',get(gcf,'color'));
    %set(gca,'ycolor',get(gcf,'color'));
    set(gca, 'Box', 'off');
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    
    ylim([0 0.4]);
    
end

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, 'imgs/greens_function_sample', 'png');

return
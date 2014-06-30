flow_vel = 1.4;      % flow velocity    [m/s]
flow_dif = 0.1;      % flow diffusivity [m^2/s]
tstep    = 86400;    % timestep size    [s]

basin = 'oh';

% basin grid georeferences
%ncols        = 104;
%nrows        = 72;
xllcorner    = -90;
yllcorner    = 34;
cellsize     = 0.125;
nodata_value = -1;

georef = [xllcorner yllcorner cellsize nodata_value];

% read flow directions
flow_dir = dlmread([basin '.dir']);
flow_dir(flow_dir<0) = NaN;
[nrows ncols] = size(flow_dir);
ncells = nansum(nansum(flow_dir*0+1));

% read gauge info
gauge_list = dlmread([basin '_gauge_list_75.txt']);

r_earth = 6371.0072 * 1000;  % earth radius [m]

% basin grid georeferences
xllcorner    = georef(1);
yllcorner    = georef(2);
cellsize     = georef(3);
nodata_value = georef(4);

[nrows ncols] = size(flow_dir);
flow_dir(flow_dir==nodata_value) = NaN;

% routing basin analysis

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

% extract sub-basins for gauge stations

gauge_lat = gauge_list(:, 3);
gauge_lon = gauge_list(:, 4) - 360;

gauge_i = nrows-floor((gauge_lat-yllcorner)/cellsize);
gauge_j = floor((gauge_lon-xllcorner)/cellsize)+1;

gauge_accum = flow_accum(gauge_i+(gauge_j-1)*nrows);
gauge_area = gauge_list(:, 5)*1609*1609;
%plot(gauge_area, gauge_accum, '.');

% remove stations with too inaccurate flow area
bad_gauges = gauge_area./gauge_accum>2 | gauge_area./gauge_accum<0.5;

gauge_list(bad_gauges, :) = [];
gauge_lat(bad_gauges) = [];
gauge_lon(bad_gauges) = [];

gauge_i(bad_gauges) = [];
gauge_j(bad_gauges) = [];

gauge_accum(bad_gauges) = [];
gauge_area(bad_gauges) = [];

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

close all;

figure('Position',[1 1 1200 900]);

% flow length background
colormap('cool');
colormap(buildcmap('mbcgyr'));
flow_length(isnan(basin_mask)) = NaN;
imagescnan(flow_length/flow_vel/tstep);
%set(h, 'AlphaData', ~isnan(basin_mask));
set(gca,'xcolor',get(gcf,'color'));
set(gca,'ycolor',get(gcf,'color'));
set(gca,'ytick',[]);
set(gca,'xtick',[]);
freezeColors;
cbfreeze(colorbar('Location', 'east', 'FontSize', 20));
text(106, 45,'Travel Time (days)', 'Rotation', 90, 'FontSize', 22, 'FontWeight', 'normal');

% plot the stream network

hold on;
stream_width = flow_accum_flat(~isnan(flow_accum_flat));
stream_width = sqrt(sqrt(stream_width/max(stream_width))) * 7;

for c = 1:ncells
    % line(lons(c,:), lats(c,:), 'LineWidth', 1, 'Color', 'black');
    line((lons(c,:)-xllcorner)/cellsize+0.5, nrows-(lats(c,:)-yllcorner)/cellsize+0.5, 'LineWidth', stream_width(c), 'Color', 'white');
end

plot3(outlet_j, outlet_i, 0.1, 'pentagram', 'MarkerSize', 20, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black');

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, 'imgs/flow_network_travel_time', 'png');
print(gcf, '-depsc2', '-painters', 'imgs/flow_network_travel_time.eps');


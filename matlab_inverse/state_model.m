% calculate all the linear operatores for the state routing model
% F maps runoff to streamflow over all points
% H maps runoff to streamflow over gauges only

function HH = state_model(basin_mask, basin_mask_gauge, flow_length_gauge, ...
    flow_vel, tstep)

[nrows ncols] = size(basin_mask);
ncells = nansum(nansum(basin_mask));
ksteps = ceil(max(max(max(flow_length_gauge)/flow_vel/tstep)));
ngauges = size(flow_length_gauge, 3);

HH = zeros(ngauges, ncells, ksteps);

% create the mapping for compacting cells
% tmp_mask_flat = basin_mask(:);
% tmp_mask_flat(isnan(tmp_mask_flat)) = 0;
% cmap = cumsum(tmp_mask_flat);
% gauge_loc = cmap(gauge_i+(gauge_j-1)*nrows);

for g=1:ngauges
    
    flow_length_gauge_compact = reshape(squeeze(flow_length_gauge(:, :, g)), nrows*ncols, 1);
    flow_length_gauge_compact(isnan(flow_length_gauge_compact)) = [];
    
    travel_steps = round(flow_length_gauge_compact/flow_vel/tstep) + 1;
    
    basin_mask_gauge_compact = reshape(squeeze(basin_mask_gauge(:, :, g)), nrows*ncols, 1);
    basin_mask_gauge_compact(isnan(basin_mask_gauge_compact)) = [];
    
    for s=1:ksteps
        HH(g, travel_steps==s, s) = 1;
        HH(g, basin_mask_gauge_compact~=1, s) = 0;
    end
    
end

return
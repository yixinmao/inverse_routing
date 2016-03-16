function route_main(basin, nst, sst, skpst)
flow_vel = 1.5;      % flow velocity    [m/s]
flow_dif = 0.1;      % flow diffusivity [m^2/s]
tstep    = 86400;    % timestep size    [s]

% basin = 'pend';

% basin grid georeferences
basedir='/raid2/ymao/VIC_RBM_east_RIPS/inverse_routing/input/';
outpdir='/raid2/ymao/VIC_RBM_east_RIPS/inverse_routing/output/';

%% check output dir
if exist([],'dir')
   error();
end


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


georef = [xllcorner yllcorner cellsize NODATA_value];

% read flow directions
flow_dir = dlmread([basedir basin '.inputs/' basin '.dir'  ]);
flow_dir(flow_dir<0) = NaN;
[nrows ncols] = size(flow_dir);
ncells = nansum(nansum(flow_dir*0+1));

% note that: ksteps = k+1 and ssteps = s+1 in the equations

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

% nsteps = 12005;
% ssteps = 80;

ssteps = str2num(sst)
spteps = str2num(skpst)
nsteps = str2num(nst)-spteps

checkd = [outpdir,basin,'/',initname,'_init_',strmname,...
          '_strm/smooth',int2str(ssteps)]
%% check output dir
if ~exist(checkd,'dir')
   error('output directory does not exist');
end


% read gauge info
gauge_list = dlmread([basedir basin '.inputs/' basin '.stn.list']);

streamflow_usgs = dlmread([basedir basin '.inputs/' basin '.stn.obs']);  % read in gauge flow data
streamflow_usgs(1:spteps,:) = [];  % delete the first #spteps days to skip
streamflow_usgs = streamflow_usgs'/35.3146662127; %% [n_gauges * nsteps]; convert cfs to m3/s

%% do basin and gauge analysis
[ basin_mask grid_lat grid_lon flow_length grid_area flow_accum ...
	ds_i ds_j dmax ds_length gauge_lat gauge_lon gauge_i gauge_j ...
	gauge_area gauge_accum basin_mask_gauge gauge_list flow_length_gauge ] ...
	= basin_analysis(flow_dir, georef, gauge_list, flow_vel, tstep, 0);

ngauges = length(gauge_i);

if (usgs_flag)
    for g=1:ngauges
        % streamflow_usgs(g, :) = streamflow_usgs(g, :) .* gauge_area(g) ./ (gauge_list(g, 5)*1609*1609) * tstep;
        % gauge_area(g) ./ (gauge_list(g, 5)*1609*1609) == 1
        streamflow_usgs(g, :) = streamflow_usgs(g, :) .* tstep;
        % This step converts streamflow_usgs to [m3/day]
    end
end

%% build the state routing model

HH = state_model(basin_mask, basin_mask_gauge, flow_length_gauge, ...
    flow_vel, tstep);

ksteps = ceil(max(max(max(flow_length_gauge)/flow_vel/tstep)));
ksteps

H = reshape(HH, ngauges, ncells*ksteps);

% create the mapping for compacting cells
tmp_mask_flat = basin_mask(:);
tmp_mask_flat(isnan(tmp_mask_flat)) = 0;
cmap = cumsum(tmp_mask_flat);
cmap(cmap==0) = 1;


%% start to route

% the first ksteps are model spin-up and the last ksteps are incompletely updated
% effective length of each smoothing window is ssteps-ksteps
nwins = floor((nsteps-ksteps*2)/(ssteps-ksteps));  % number of windows

runoff1_compact = zeros(ncells*ksteps, 1);
runoff2_compact = zeros(ncells*ksteps, 1);

runoff1_combine = zeros(ncells*ssteps, 1);
runoff2_combine = zeros(ncells*ssteps, 1);

runoff2_save_post = zeros(ncells*(ssteps-ksteps), nwins);

Hprime = sparse(ngauges*ssteps, ncells*ssteps);

streamflow_gauge1 = zeros(ngauges, ssteps);
streamflow_gauge2 = zeros(ngauges, ssteps);
streamflow_gauge2_post = zeros(ngauges, ssteps);

streamflow_gauge1_all = zeros(ngauges, nsteps-ksteps);
streamflow_gauge2_all = zeros(ngauges, nsteps-ksteps);
streamflow_gauge2_post_all = zeros(ngauges, nsteps-ksteps);


% assemble H'
for s=1:ssteps
    if (s+ksteps<=ssteps)
        % too bad the following works too slow
        Hprime(ngauges*(s-1)+1:ngauges*s, ncells*(s-1)+1:ncells*(s-1+ksteps)) = H;
    else
        % too bad the following works too slow
        Hprime(ngauges*(s-1)+1:ngauges*s, ncells*(s-1)+1:ncells*ssteps) = H(:, 1:ncells*(ssteps-s+1));
    end
end

%
[basedir, basin, '.inputs/', basin, '.runoff']  % hacked by Yixin
fin2 = fopen([basedir, basin, '.inputs/', basin, '.runoff'], 'r');

% read everything all in once to save seek/rewind time

% skip the first spteps days
for skipi=1:spteps
    runoff2 = fscanf(fin2, '%f', [nrows ncols]);
end

runoff2_input = zeros(ncells, nsteps);
for s=1:nsteps

    % read runoff
    runoff2 = fscanf(fin2, '%f', [nrows ncols]);
    runoff2 = runoff2.* basin_mask/1000.*grid_area;  % [m3/day]
    fillval = mean(runoff2(runoff2>0));
    runoff2(runoff2<=0) = fillval;

    runoff2_tmp_compact = runoff2(:);
    runoff2_tmp_compact(isnan(runoff2_tmp_compact)) = [];  % [ncells*1]
    runoff2_input(:, s) = runoff2_tmp_compact;
    
end

fclose(fin2);

% model spin-up

for s=1:ksteps
    
    % shift for one step
    runoff2_compact = circshift(runoff2_compact, [ncells 0]);
    runoff2_compact(1:ncells) = runoff2_input(:, s);
    
end

% routing and assimilation
% nwins=4;
for w=1:nwins
    
    w
    % record routing initial condition
    if (w==1)
        runoff2_init_compact = runoff2_compact;
    else   % if w>1, rewind for ksteps
        runoff2_compact = runoff2_init_compact;
    end
    
    for s=1:ssteps
    
        s_global = ksteps+(ssteps-ksteps)*(w-1)+s;
        
        % shift for one step
        runoff2_compact = circshift(runoff2_compact, [ncells 0]);
        runoff2_compact(1:ncells) = runoff2_input(:, s_global);
    
        streamflow_gauge2(:, s) = H * runoff2_compact;  % this is one "big row" in ( H'xt' + L'x(t-k)' )
                                                        % (each "big row" has ngauges number of rows)
               % When putting into streamflow_gauge2, the "big row" becomes
               % "big column"

        % filling up H' and x'
        runoff2_combine(ncells*(ssteps-s)+1:ncells*(ssteps-s+1)) = runoff2_input(:, s_global);
                % This is the current day;

%         % save the runoff data vector ksteps prior to the end of the smoothing window
%         if (s==ssteps-ksteps)
%             runoff2_compact_last = runoff2_compact;
%         end
        
    end

    % Kalman

    % errors are proportional to runoff magnitude        
    P = sparse(1:ncells*ssteps, 1:ncells*ssteps, runoff2_combine.*runoff2_combine+fillval*fillval*100, ncells*ssteps, ncells*ssteps, ncells*ssteps);

    % Kalman gain K is: K = P * Hprime' * inv(Hprime * P * Hprime');    
    A = sparse(P * Hprime');
    B = sparse(Hprime * A);
    K = inv(B);
    K = A * K;
    
    if (~usgs_flag)
%%%        innov = reshape(fliplr(streamflow_gauge1-streamflow_gauge2), ngauges*ssteps, 1);
    else
        innov = reshape(fliplr(streamflow_usgs(:, ksteps+(ssteps-ksteps)*(w-1)+1:ksteps+(ssteps-ksteps)*w+ksteps)-streamflow_gauge2), ngauges*ssteps, 1);
    end

    % Kalman Update
    runoff2_combine_post = runoff2_combine + K*innov;
    
%    runoff2_combine_post(runoff2_combine_post<0) = 0;  % Hacked by Yixin
    
    % recalculate streamflow using updated runoff fields and
    % reset routing initial condition for next smoothing window
    runoff2_compact = runoff2_init_compact;
    for s=1:ssteps
        % shift for one step
        runoff2_compact = circshift(runoff2_compact, [ncells 0]);
        runoff2_compact(1:ncells) = runoff2_combine_post((ssteps-s)*ncells+1:(ssteps-s+1)*ncells);
        streamflow_gauge2_post(:, s) = H * runoff2_compact;
        
        % save the runoff data vector ksteps prior to the end of the smoothing window
        if (s==ssteps-ksteps)
            runoff2_init_compact = runoff2_compact;
        end
    end
    
    if (w==1)
        streamflow_gauge2_all(:, ssteps*(w-1)+1:ssteps*w) = streamflow_gauge2; % this is ( H'xt' + L'x(t-k)' ), 
                                                                            % but reshaped to [ngauges*ssteps]
        streamflow_gauge2_post_all(:, ssteps*(w-1)+1:ssteps*w) = streamflow_gauge2_post;
    else
        streamflow_gauge2_all(:, ksteps+(ssteps-ksteps)*(w-1)+1:ksteps+(ssteps-ksteps)*w) = streamflow_gauge2(:, ksteps+1:ssteps);
        streamflow_gauge2_post_all(:, ksteps+(ssteps-ksteps)*(w-1)+1:ksteps+(ssteps-ksteps)*w) = streamflow_gauge2_post(:, ksteps+1:ssteps);
    end
        
    % test: resetting it to "true" initial condition for next window
    % note: in real assimilatin experiment, this is not possible.
    %if (tmpa_flag)
    %    runoff2_compact = runoff1_compact;
    %end
    
    % Final runoff results: ignore the last ksteps in the window
    runoff2_save_post(:, w) = runoff2_combine_post(ncells*ksteps+1:ncells*ssteps);
    
end

save([outpdir basin '/' initname '_init_' strmname '_strm/smooth' int2str(ssteps) '/all_data.mat'], ...
    'streamflow_gauge1_all', 'streamflow_gauge2_all', 'streamflow_gauge2_post_all', 'streamflow_usgs', ...
    'ncols', 'nrows', 'ncells', 'flow_vel', 'tstep', 'nsteps', 'ksteps', 'ssteps', 'nwins', 'cmap', ...
    'basin_mask', 'grid_area', 'gauge_list', 'gauge_area', 'initname', 'strmname');

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
    st=1+(i-1)*(ssteps-ksteps);
    ed=i*(ssteps-ksteps);
    data_2_w(:,st:ed)=fliplr(data_2_w(:,st:ed));
end

[bm,bn]=size(basin_mask);
[m,n]=find(~isnan(basin_mask));
lon=0.5*cellsize+xllcorner+(n-1)*cellsize;
lat=0.5*cellsize+yllcorner+(bm-1)*cellsize-(m-1)*cellsize;
lon_lat_data=[lon,lat,data_2_w];

dlmwrite([outpdir,basin,'/data_all_day_',num2str(ssteps),...
          '_skip',skpst],lon_lat_data,'delimiter', '\t', ...
          'newline', 'unix','precision', 9);

end


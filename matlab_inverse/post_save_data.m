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

%% fliplr the output
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
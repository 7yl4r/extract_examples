% RSCLASS_EXTRACT.M
% Written by Dan Otis, March 2019; revised August 2023
% Generic extraction script to extract satellite data
% Extracts data from a directory/list of .nc files
% User inputs:
% Input/output directories, parameter to examine, locations to extract
% Location(s):
%    a. Point(s)
%    b. Polygon(s)
%    c. Shapefile(s)

clear

% Set filepaths, lat/lon limits, and x/y sizes
path_main='/Users/imars_mbp/Box Sync/ML_working';

% Define product(s) to be extracted
prod='adg_443';

% INPUT AND OUTPUT PATHS
file_path='/Users/imars_mbp/Box Sync/RS_Class_2023/Share/NetCDF_examples/4km_adg443_MO_means/';
path_out='/Users/imars_mbp/Box Sync/RS_Class_2023/Share/NetCDF_examples/4km_adg443_MO_means/out/';

% Get info from filenames
eval(['dirinfo_tmp=struct2cell(dir(''' file_path '/*.nc''));'])
filenames_tmp=dirinfo_tmp(1,:)'; % filenames are contained in top row
filenames_str=char(filenames_tmp);

% CONVERT TIME
year_st = str2num(filenames_str(:,12:15));
month_st = str2num(filenames_str(:,16:17));
day_st = str2num(filenames_str(:,18:19));

mltime = datenum(year_st,month_st,day_st,0,0,0);
numfiles = size(filenames_str,1);

% Get lat/lon from test file
file=filenames_str(1,:);
cd(file_path)
[lat]=open_nc(file,'lat');
lat=lat'; % Transpose latitude
[lon]=open_nc(file,'lon');
cd(path_main)

%%%%% LOCATIONS FOR DATA EXTRACTION %%%%%%
% Options for extraction:
% 1. Discreet point(s) - entered manually or loaded from file
% 2. Polygon (from shapefile)
% 3. Polygon (points entered manually)
% 4. Polygon or points loaded from a file (.csv, .mat, etc.)

% Extract from single point sites here (use 3x3 or 5x5 pixel boxes)
% Four points defined here (N. Atl., S. Atl., GoM, Southern Ocean)
loc_lat=[42.88,-31.04,25.65,-65.51];
loc_lon=[-28.62,-13.46,-89.85,28.56];

% Convert geo coordinates to image (array) coordinates
[loc_pix,loc_line]=latlon2pixline(loc_lat,loc_lon,lat(:,1),lon(1,:));

% Select offset to extract a box
bx_offset=50; % Gives 51x51 pixel box. Use 1 for 3x3 box; 2 for 5x5 box, etc.

% Create a set of "extraction masks" where we want to extract data
for s=1:length(loc_lat)
% Create array of zeros the same size as image
loc_pts_tmp=false(size(lat,1),size(lon,2));
% Define a "box" around the center pixel (size based on offset)
loc_pts_tmp(loc_line(s)-bx_offset:loc_line(s)+bx_offset,loc_pix(s)-bx_offset:loc_pix(s)+bx_offset)=1;
% Create 3-D array of extraction "masks" (ones where we want data, zeroes where we don't)
loc_pts(:,:,s)=loc_pts_tmp;
end

% Extract from polygon (user must define points)
lat_poly=[50,50,30,30,50];
lon_poly=[-170,-150,-150,-170,-170];

[lon_matrix,lat_matrix]=meshgrid(lon,lat);
poly_pts = inpolygon(lon_matrix,lat_matrix,lon_poly,lat_poly); % poly_pts is created as a mask like the others

% Read data from a shapefile (Moneterey Bay National Marine Sanvctuary boundaries)
% M = m_shaperead('/Users/IMaRS_iMAC/Box Sync/RS_Class_2019/Share/mbnms_py2/mbnms_py');
M = m_shaperead('/Users/imars_mbp/Box Sync/RS_Class_2023/Share/mbnms_py2/mbnms_py');

% Define points (points are contained in "ncst" portion of structure variable(M)
shp_tmp=M.ncst{1,1};
shp_lat=shp_tmp(:,2);
shp_lon=shp_tmp(:,1);
shp_pts = inpolygon(lon_matrix,lat_matrix,shp_lon,shp_lat); % inpolygon finds points within the given polygon

% Quick plot with test file to check location(s) of extracted data
% Define test image
test_file=filenames_str(1,:);
cd(file_path)
img_test=double(open_nc(test_file,'adg_443'));
% Apply scale factor and offset
ncid = netcdf.open(test_file,'NC_NOWRITE'); % Define "ncid" as the file identifier
p_id = netcdf.inqVarID(ncid,'adg_443');
scl = netcdf.getAtt(ncid,p_id,'scale_factor','double');
offset = netcdf.getAtt(ncid,p_id,'add_offset','double');
img_scl = img_test.*(scl) + offset;
img_test(img_test <= 0) = NaN; % Remove neg. vals.
cd(path_main)

% Loop on point locations
for s=1:length(loc_lat)
img_scl(loc_pts(:,:,s))=1; % Set image vals to 1 over point extraction masks
end
img_scl(poly_pts)=1; % Polygon extraction mask
img_scl(shp_pts)=1; % Shapfile extraction mask

figure('Position', [800, 100, 1600, 800]);
colormap(jet)
clims = [0.001 0.05];
imagesc(img_scl(:,:,1),clims);


%%%% MAIN LOOP %%%%
%%%% LOOP ON FILES TO EXTRACT %%%%
for i=1:numfiles
file=filenames_str(i,:);
cd(file_path)
eval(['[img_tmp]=double(open_nc(file,''' prod '''));'])
img_tmp = img_tmp.*(scl) + offset;
img_tmp(img_tmp <= 0) = NaN; % Remove negative vals.

% Point location extraction (4 locations)
for s=1:length(loc_lat)
pt_mean=nanmedian(img_tmp(loc_pts(:,:,s)));
pt_std=nanstd(img_tmp(loc_pts(:,:,s)));
pt_ts_mn(i,s)=pt_mean;
pt_ts_std(i,s)=pt_std;
end

% Polygon extraction (N. Pacific)
poly_ts_mn(i)=nanmedian(img_tmp(poly_pts));
poly_ts_std(i)=nanstd(img_tmp(poly_pts));

% Shapefile extraction (Monterey Bay NMS)
shp_ts_mn(i)=nanmedian(img_tmp(shp_pts));
shp_ts_std(i)=nanstd(img_tmp(shp_pts));
end

% Clean up
clear img_tmp

% Plotting
% Convert time to datetime
dt_time=datetime(mltime,'ConvertFrom','datenum');

figure('Position', [800, 100, 1600, 800],'Name','2015 monthly Chl time-series');
subplot(3,1,1)
plot(dt_time,pt_ts_mn,'.-','Linewidth',2)
legend('N. Atl.','S. Atl.','GoM','Southern Ocean')
subplot(3,1,2)
plot(dt_time,poly_ts_mn,'.-','Linewidth',2)
title('N. Pacific polygon')
subplot(3,1,3)
plot(dt_time,shp_ts_mn,'.-','Linewidth',2)
title('Monterey Bay NMS shapefile')

% Output
% Many options here:
% .csv file(s) - Use csvwrite
% .xls file(s) - Use xlswrite
% Text file(s) - Can explicitly define output using fprintf
% .mat file(s) - Use "save" command

% cd(path_out)
% save RS_extract_out.mat pt_ts_mn pt_ts_std poly_ts_mn poly_ts_std shp_ts_mn shp_ts_std

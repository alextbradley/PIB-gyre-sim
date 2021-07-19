%Make a data output file for the gyre train simulation (ISOMIP run with shelf beyond 180km removed) and extended (1000km) domain.
% Data file will contain:
% 'X':  	(1 x nx) array
%		X co-ordinates of grid points (nx: number of grid pts in x)
%
% 'Y':  	(1 x ny) array 
%		Y co-ordinates of grid points (ny: number of grid pts in y)
%
% 'Z'		(1 x nz) array
% 		Z co-ordinates of grid points (nz: number of grid pts in z)
%
% 'bathy': 	(nx x ny) array
%		Sea bathymetry. Zeros indicate locations of solid walls 
%
% 'topo':	(nx x ny) array
%		Ice shelf draft. Zeros indicate locations where no ice shelf
%
% 'melt': 	(nx x ny) array 
% 		Simulated melt rate at grid points (units: m/yr)
%
% 'U': 		(nx x ny) array 
%		Simulated depth averaged velocities in x direction at grid points
%
% 'V':  	(nx x ny) array
% 		Simulated depth averaged velocities in y direction at grid points
%
% 'bsf':	(nx x ny) array
% 		Computed barotropic streamfunction from simulation output
%
% All simulated output averaged over the final 3 years of a 36 year simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Preliminaires
%
%input info
savename = 'gyre_train_data_shiC1e-3.mat';
run_nos = "036";
extent = 18;


% Data locations 
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/AISOMIP_';
topodir = '/data/hpcdata/users/aleey/mitgcm/matlab/interp_AISOMIP/extended_domain_files/extended_snapped';
bathy_path = '/data/hpcdata/users/aleey/mitgcm/matlab/interp_AISOMIP/extended_domain_files/bathy_extendedDom.shice';

%run info and flags
state2D_fname = strcat(rootdir, run_nos, '/run/MITout_2D.nc');
state3D_fname =  strcat(rootdir, run_nos, '/run/MITout_3D.nc');
nt0 = 66; %first time step out 
nt  = 4;  %number of timesteps out (67, 68, 69, 70, 71, 72)

%grid 
nx=480; % number of grid cells along longitudinal direction
ny=42; % number of grid cells along latitudinal direction
nz=36; % number of vertical grid cells
dx=2000;
dy=2000;
dz=20;
X = 0:dx:(nx-1)*dx;
Y = 0:dx:(ny-1)*dy;
Z = 0:dz:(nz-1)*dz;

%parameters
secs_per_year = 360*24*60*60;
density_ice = 918.0;

% 
% Get input data
%

%bathymetry 
fid = fopen(bathy_path);
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx, ny]);

%topography
topo_fname = ['shelfice_topo_extendedDom_extent' num2str(extent) 'km.bin'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx, ny]);


%
% Get simulated data 
%

%velocities
U3d =  ncread(state3D_fname, 'UVEL', [1, 1, 1, nt0], [Inf, Inf,Inf, nt]);
U3d = mean(U3d,4); %time average
U   = mean(U3d,3); %depth average
V3d =  ncread(state3D_fname, 'VVEL', [1, 1, 1, nt0], [Inf, Inf,Inf, nt]);
V3d = mean(V3d,4); %time average
V   = mean(V3d,3); %depth average

%melt rate
melt = ncread(state2D_fname, 'SHIfwFlx', [1,1,nt0], [Inf,Inf,nt]);
melt = mean(melt, 3);%time average
melt = -melt * secs_per_year / density_ice;

%bsf
bsf=zeros(nx,ny);
vv = squeeze(sum(V3d, 3)) * dz; 
bsf(nx,:)=vv(nx,:)*dx;
for p=nx-1:-1:1
 bsf(p,:)=bsf(p+1,:) + vv(p,:)*dx;
end
bsf = bsf / 1e6; %convert to sv

%
% save 
%
save(savename, 'bathy', 'bsf', 'melt', 'topo', 'U', 'V', 'X', 'Y', 'Z');





%

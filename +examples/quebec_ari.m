% Script to play for Ariane
%found it!!
clear all
clc

%% Path management
RF = 'C:\Users\arian\Documents\internship'; %RootFolder
addpath(genpath(strcat(RF,'/git/adcptools'))); %path to ADCPTools
addpath(genpath("C:\Users\arian\Documents\internship\git\adcptools\post_processing"))
% addpath(genpath(strcat(RF,'Tools\adcptools'))); %possible other folders

%% Quick documentation walkthrough - comment out

%open_adcptools_documentation()

%% Constituents

constituents = {'M2', 'M4'};

%% Loading in the data
% addpath('./Donnees_validation'); %path to data
%dat = rdi.readDeployment('rijn', './data');
%dat = rdi.readDeployment('Lauzon_0_0', './data/Lauzon_0');

% addpath('./data'); %path to data
%dat = rdi.readDeployment('rijn', './data');
dat = rdi.readDeployment('Quebec_0_0','C:\Users\arian\Documents\internship\Donnees_validation\2009\ADCP 2009\Quebec_0');

%% test dat frequency results
iseven = 1:50:30430;
 dat.FileNumber = dat.FileNumber(:,iseven);
 dat.firmver = dat.firmver(:,iseven);
 dat.firmrev = dat.firmrev(:,iseven);
 dat.sysconf = dat.sysconf(:,iseven);
 dat.sysconfstr = dat.sysconfstr(:,iseven);
 dat.SymData = dat.SymData(:,iseven);
 dat.LagLength = dat.LagLength(:,iseven);
 dat.usedbeams = dat.usedbeams(:,iseven);
 dat.nbins = dat.nbins(:,iseven);
 dat.pingperens = dat.pingperens(:,iseven);
 dat.binsize = dat.binsize(:,iseven);
 dat.blnk = dat.blnk(:,iseven);
 dat.minthrsh = dat.minthrsh(:,iseven);
 dat.ncodrep = dat.ncodrep(:,iseven);
 dat.minpercgood = dat.minpercgood(:,iseven);
 dat.maxerrvel = dat.maxerrvel(:,iseven);
 dat.Tbetweenpng = dat.Tbetweenpng(:,iseven);
 dat.corinfo = dat.corinfo(:,iseven);
 dat.corstr = dat.corstr(:,iseven);
 dat.headalign = dat.headalign(:,iseven);
 dat.headbias = dat.headbias(:,iseven);
 dat.sensource = dat.sensource(:,iseven);
 dat.senavail = dat.senavail(:,iseven);
 dat.distmidbin1 = dat.distmidbin1(:,iseven);
 dat.lngthtranspulse = dat.lngthtranspulse(:,iseven);
 dat.watrefbins = dat.watrefbins(:,iseven);
 dat.mintarget = dat.mintarget(:,iseven);
 dat.lowlattrig = dat.lowlattrig(:,iseven);
 dat.distpulse = dat.distpulse(:,iseven);
 dat.cpuserial = dat.cpuserial(:,iseven);
 dat.bandwidth = dat.bandwidth(:,iseven);
 dat.syspower = dat.syspower(:,iseven);
 dat.basefreqid = dat.basefreqid(:,iseven);
 dat.serial = dat.serial(:,iseven);
 dat.HADCPbeamangle = dat.HADCPbeamangle(:,iseven);
 dat.VEL = dat.VEL(:,iseven,:);
 dat.ECHO =  dat.ECHO(:,iseven,:);
 dat.CORR = dat.CORR(:,iseven,:);
 dat.PERC = dat.PERC(:,iseven,:);
 dat.btpingperens = dat.btpingperens(:,iseven);
 dat.reacqdelay = dat.reacqdelay(:,iseven);
 dat.mincormag = dat.mincormag(:,iseven);
 dat.minevampl = dat.minevampl(:,iseven);
 dat.btminpergood = dat.btminpergood(:,iseven);
 dat.btmode = dat.btmode(:,iseven);
 dat.btmaxerrv = dat.btmaxerrv(:,iseven);
 dat.btrange = dat.btrange(iseven,:);
 dat.btvel = dat.btvel(iseven,:);
 dat.btcor = dat.btcor(iseven,:);
 dat.btevampl = dat.btevampl(iseven,:);
 dat.btpercgood = dat.btpercgood(iseven,:);
 dat.minlyrsize = dat.minlyrsize(:,iseven);
 dat.nearbnd = dat.nearbnd(:,iseven);
 dat.farbnd = dat.farbnd(:,iseven);
 dat.reflyrvel = dat.reflyrvel(iseven,:);
 dat.reflyrcor = dat.reflyrcor(iseven,:);
 dat.reflyrint = dat.reflyrint(iseven,:);
 dat.reflyrpergood = dat.reflyrpergood(iseven,:);
 dat.maxdepth = dat.maxdepth(:,iseven);
 dat.rssiamp = dat.rssiamp(iseven,:);
 dat.gain = dat.gain(:,iseven);
 dat.sp_btrange_fract = dat.sp_btrange_fract(iseven,:);
 dat.ensnum = dat.ensnum(:,iseven);
 dat.BITcheck = dat.BITcheck(iseven,:);
 dat.speedsound = dat.speedsound(:,iseven);
 dat.depthtransd = dat.depthtransd(:,iseven);
 dat.heading = dat.heading(:,iseven);
 dat.pitch = dat.pitch(:,iseven);
 dat.roll = dat.roll(:,iseven);
 dat.salinity = dat.salinity(:,iseven);
 dat.temperature = dat.temperature(:,iseven);
 dat.prepingT = dat.prepingT(:,iseven);
 dat.headstd = dat.headstd(:,iseven);
 dat.pitchstd = dat.pitchstd(:,iseven);
 dat.rollstd = dat.rollstd(:,iseven);
 dat.ADC = dat.ADC(iseven,:);
 dat.errorstat1 = dat.errorstat1(iseven,:);
 dat.errorstat2 = dat.errorstat2(iseven,:);
 dat.errorstat3 = dat.errorstat3(iseven,:);
 dat.errorstat4 = dat.errorstat4(iseven,:);
 dat.pressure = dat.pressure(:,iseven);
 dat.pressurevar = dat.pressurevar(:,iseven);
 dat.timeV = dat.timeV(iseven,:);
 dat.timeV1C = dat.timeV1C(iseven,:);
 dat.NMEAGGA.deltaT = dat.NMEAGGA.deltaT(iseven,:);
 dat.NMEAGGA.msgHeader = dat.NMEAGGA.msgHeader(iseven,:);
 dat.NMEAGGA.UTC = dat.NMEAGGA.UTC(iseven,:);
 dat.NMEAGGA.Lat = dat.NMEAGGA.Lat(iseven,:);
 dat.NMEAGGA.SN = dat.NMEAGGA.SN(iseven,:);
 dat.NMEAGGA.Long = dat.NMEAGGA.Long(iseven,:);
 dat.NMEAGGA.EW = dat.NMEAGGA.EW(iseven,:);
 dat.NMEAGGA.Qual = dat.NMEAGGA.Qual(iseven,:);
 dat.NMEAGGA.NSat = dat.NMEAGGA.NSat(iseven,:);
 dat.NMEAGGA.HDOP = dat.NMEAGGA.HDOP(iseven,:);
 dat.NMEAGGA.Alt = dat.NMEAGGA.Alt(iseven,:);
 dat.NMEAGGA.AltUnit = dat.NMEAGGA.AltUnit(iseven,:);
 dat.NMEAGGA.Geoid = dat.NMEAGGA.Geoid(iseven,:);
 dat.NMEAGGA.GeoidUnit= dat.NMEAGGA.GeoidUnit(iseven,:);
 dat.NMEAGGA.AgeDGPS = dat.NMEAGGA.AgeDGPS(iseven,:);
 dat.NMEAGGA.RefStatID = dat.NMEAGGA.RefStatID(iseven,:);
 dat.nFiles.GGA.utc = dat.nFiles.GGA.utc(iseven,:);
 dat.nFiles.GGA.latitude = dat.nFiles.GGA.latitude(iseven,:);
 dat.nFiles.GGA.longitude = dat.nFiles.GGA.longitude(iseven,:);
 dat.nFiles.GGA.quality = dat.nFiles.GGA.quality(iseven,:);
 dat.nFiles.GGA.numsat = dat.nFiles.GGA.numsat(iseven,:);
 dat.nFiles.GGA.hdop = dat.nFiles.GGA.hdop(iseven,:);
 dat.nFiles.GGA.alt = dat.nFiles.GGA.alt(iseven,:);
 dat.nFiles.GGA.geoid = dat.nFiles.GGA.geoid(iseven,:);
 dat.nFiles.GGA.age_dgps = dat.nFiles.GGA.age_dgps(iseven,:);
 dat.nFiles.GGA.ref_station_id = dat.nFiles.GGA.ref_station_id(iseven,:);
%% Load water level data
load('C:\Users\arian\Documents\internship\Donnees_validation\2009\marégraphes_h_2009_HNE_NMM_3min\marégraphes_h_2009_HNE_NMM_3min\3250Lauzon2009_HNE_NMM_3min.mat')


%% waterlevel
filt = ~isnan(h);
water_level = VaryingWaterLevel(datetime(t(filt), 'ConvertFrom', 'datenum'), h(filt));
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();


%% Modify the following code to analyze the data

V = rdi.VMADCP(dat);
% V.horizontal_position_provider = HorizontalPositionFromBottomTracking; % possibly modify

 V.water_level_object = water_level;  % return

B = BathymetryScatteredPoints(V);

Bfilt = find(B.known(2,:)>0);

B.interpolator.span = .001;
figure;
B.plot

V.filters = Filter;
%V.shipvel_provider = ShipVelocityFromBT; % possibly modify


%figure;
%hold on
%V.plot_all

[ef, xs] = cross_section_selector(V);

%% Mesh for plotting

mesh_makers = SigmaZetaMeshFromVMADCP(ef, xs, B, 'NoExpand', V);

% input preferred mresh size

hver = 5; % depth mesh cell in m
hhor = 25; %width mesh cell in m

%% select max & min in that order

figure
plot(V.horizontal_position(1,:))
[~, x] = ginput;
maxx = x(1,1);
minx = x(2,1);


plot(V.horizontal_position(2,:));
[~, y] = ginput;
maxy = y(1,1);
miny = y(2,1);

%% calculations

lengthriv = sqrt((maxy-miny)^2+(maxx-minx)^2); 

n = round(lengthriv/hhor);

acthor = lengthriv/n;

depthriv = mean((max(V.bt_vertical_range,[], 'omitnan')));

z = round(depthriv/hver);
actver = depthriv/z;

fprintf('Used horizontal mesh size is: %.2f\n', acthor);
fprintf('Used vertical mesh size is: %.2f\n', actver);

mesh = mesh_makers.get_mesh(resn = n, resz = z);

%% Model
opts = SolverOptions(extrapolate_vert = 0, lat_weight_factor = 10); % possibly modify
%opts.force_zero = [1 1 1 1 1];

% Empirical model: VelocityModel;
flow_model = TaylorTidalVelocityModel; % possibly modify to enter desired empirical model formulation
flow_model.constituents = constituents;

% or TaylorVelocityModel
flow_model.s_order =[1 1 1]; %u
flow_model.n_order = [1 1 1]; %v
flow_model.sigma_order = [1 1 1]; %sja


%Solver options and regularization
flow_regs = regularization.Velocity.get_all_regs(mesh, B, xs, flow_model, opts, 'NoExpand', V);


% Bulk regularization parameter % possibly modify
lc = 0.1;
flow_regs(1).weight =  lc;
flow_regs(2).weight =  lc;
flow_regs(3).weight =  lc;
flow_regs(4).weight =  lc;
flow_regs(5).weight =  lc;



% Solve for the flow
flow_solv = LocationBasedVelocitySolver(mesh, B, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
flow_solv(1).rotation = xs.angle;
flow = flow_solv.get_solution(); % possibly modify
%figure
%spy(flow.M)
% Plot the state vector
flow.plot_solution()

%%
%     err = cross_validate_0D(flow);
% [RMSE] = plot_mrse_mesh(mesh, 'V', 2, pars_V, t_plot, constituents, xs, V, channel);
%     cross_validate_1D(flow, 0, 1000, 1000)
for i = 8
 flow.solver.opts.cv_mode = 'track';
 [rho, lambda,eta] = cross_validate_1D_track(flow, 0, 1, 1000,i);
end
  
%%
rho = cell2mat(rho);
eta = cell2mat(eta);
[r,e,l] = l_corner(rho2(:,1),eta2(:,1),lambda2(:,1),2,18);
%%
[CV, rpc, rps] = flow.cross_validate_2D_track([0,0], [1, 1e12], [20,20], 2,2);

figure;
%make use of the ordering of the constraints: first continuity

contourf(helpers.symlog(rpc), helpers.symlog(rps), reshape([CV{:,1}]./CV{1,1}, 20,20), 100)
colorbar
colormap(flipud(helpers.cmaps('velmap')))
% clim([0,2]) % 0 - very good (too good to be true) 1
xlabel('cont lambda (symlog10)')
ylabel('smoothness lambda (symlog10)')
title('2D cross-validation: generalization error')
%% Post-Processing - focus on decomposition of the solution
addpath('C:\Users\Ariane.VandePas\Documents\GitHub\adcptools')
addpath('C:\Users\Ariane.VandePas\Documents\GitHub\adcptools\+post_processing\cartesian')
addpath(genpath(("C:\Users\Ariane.VandePas\Documents\GitHub\adcptools\+post_processing\plot")))
tim = flow.solver.adcp.time;

Tlim(1)= min(tim);

M2T = flow.solver.model.periods(1,1)/(3600*24);
Tlim(2) = Tlim(1) + M2T;
Tlimn = datenum(Tlim);

limu = {[Tlimn(1), Tlimn(2)];... %days (!)
    [min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5];...
    [0,1]};

% evaluation resolution
tres = 60;
evres = [tres, flow.solver.mesh.nverticals, flow.solver.mesh.max_ncells_vertical]; % t, y , sigma
X = get_coords(limu, evres);


reg_idx = 1; % only relevant if multiple regularization parameter settings are entered upon model fitting.
u = get_var(flow, X, reg_idx); % Vector variable on regular sigma grid. Three cells are the three velocity components.

% Conversion to 'semi-Cartesian' coordinates: time, lateral, sigma - t, y,
% sig

[H, Wl, Zb] = get_H(X, flow, 0);

if size(H, 3) == 1
    H = repmat(H, [1,1,numel(X.sig)]);
end
X.Z = Zb + X.Sig.*H;


D = post_processing.Decomposition(X = X, H = H, wl = Wl(:,1), zb = Zb(1,:)');

% Plot some variables
name = 'flow550_01reg.gif';
sav = 1;
animate_solution(u{1}, X, name, sav)

[u_decomp, u_avg] = D.decompose_function(u{2}); % U-Flow

D.plot_components(u_decomp, 'velmap')
D.plot_components(u_avg, 'velmap')


%% Post-Processing - focus on regularization parameters - NEEDS REVISION
% Preliminary K-fold cross-validation. How well does the solution fit
% unseen data?

% Regularization does:
% - make the solution more 'regular', often more smooth
% - help in interpolating and extrapolating
% - improve robustness to noise (according to some metrics)
% - introduce bias wrt the data alone, to reduce variance (see bias-variance decomposition)

% Regularization does not:
% - fix outliers - these still can play a large role
% - obtain a more 'true' solution
% - magically fix all problems
% - the regularization constraints also contain numerical and
%       discretization errors

reg_pars_mat = repmat([0, logspace(-5,3,2)]', 1, 5);

flow.opts.training_perc = .66;
flow.opts.cv_iter = 1; %1-fold cross validation (See Brunton & Kutz)
CV = flow.cross_validate_single(reg_pars_mat); % todo: wrapper for 2D sensitivity like figs from ADCPpaper


reg_pars_plot = reg_pars_mat(:,1);
reg_pars_plot(1) = reg_pars_plot(2)/10;
%workaround for plotting 0 -> set at logscale factor 10 back from second
%smallest value.
figure;
semilogx(reg_pars_plot', [CV{:,1}])
xlabel('reg pars')
ylabel('generalization error')
title('lambda vs scaled generalization error')

%% trying to plot raw data

    vel_pos={V.depth_cell_position};
    vel_pos=cellfun(@(x) mean(x(:,:,:,3),3,'omitnan'), vel_pos,...
        'UniformOutput', false);

    vel_xy={V.water_velocity(CoordinateSystem.Earth)};
    tim=V.time;

    vel_sn = {nan(size(vel_xy{1}))};
    vel_sn{1}(:,:,3) = vel_xy{1}(:,:,3);
    vel_sn{1}(:,:,4) = vel_xy{1}(:,:,4);

    R = [xs.direction_orthogonal(1), xs.direction_orthogonal(2); xs.direction(1),  xs.direction(2)];

    for i = 1:size(vel_xy{1}(:,:,1),1)
        for j = 1:size(vel_xy{1}(:,:,1),2)
            velocity_U = vel_xy{1}(i,j,1);
            velocity_V = vel_xy{1}(i,j,2);
            vel_sn{1}(i,j,1) = R(1,1) * velocity_U + R(1,2) * velocity_V;
            vel_sn{1}(i,j,2) = R(2,1) * velocity_U + R(2,2) * velocity_V;

            % disp(['velocity_U: ', num2str(velocity_U), ' m/s']);
            % disp(['velocity_V: ', num2str(velocity_V), ' m/s']);
            % disp(['U vel_U_sn: ', num2str(vel_sn{1}(i,j,1)), ' m/s']);
            % disp(['V vel_V_sn: ', num2str(vel_sn{1}(i,j,2)), ' m/s']);
        end
    end



    %% stoopid code

function X = get_coords(lim, res)
t = linspace(lim{1,1}(1), lim{1,1}(2), res(1)); % in days
n = linspace(lim{2,1}(1), lim{2,1}(2), res(2));
sig = linspace(lim{3,1}(1), lim{3,1}(2), res(3));

[T, N, Sig] = ndgrid(t, n, sig);

X.t = t; X.y = n; X.sig = sig;
X.T = T; X.Y = N; X.Sig = Sig;
end


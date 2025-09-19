% Script to play for Ariane
%found it!!
clear all
clc

%% Path management
RF = 'C:\Users\arian\Documents\internship'; %RootFolder
addpath(genpath(strcat(RF,'/git/adcptools'))); %path to ADCPTools
% addpath(genpath(strcat(RF,'Tools\adcptools'))); %possible other folders

%% Quick documentation walkthrough - comment out

%open_adcptools_documentation()

%% Constituents

constituents = {'M2', 'M4'};

%% Loading in the data
addpath('./Donnees_validation'); %path to data
%dat = rdi.readDeployment('rijn', './data');
%dat = rdi.readDeployment('Lauzon_0_0', './data/Lauzon_0');

% addpath('./data'); %path to data
%dat = rdi.readDeployment('rijn', './data');
dat = rdi.readDeployment('Quebec_0_0', './2009/ADCP 2009/Quebec_0');
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
V.horizontal_position_provider = HorizontalPositionFromBottomTracking; % possibly modify

V.water_level_object = water_level;

B = BathymetryScatteredPoints(V);

%Bfilt = find(B.known(2,:)>0);

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

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, B, 'NoExpand', V);

mesh = mesh_maker.get_mesh(resn = 50, resz = 15);

%% Model
opts = SolverOptions(extrapolate_vert = 0, lat_weight_factor = 10); % possibly modify
%opts.force_zero = [1 1 1 1 1];

% Empirical model: VelocityModel;
flow_model = TaylorTidalVelocityModel; % possibly modify to enter desired empirical model formulation
flow_model.constituents = constituents;

%or TaylorVelocityModel
flow_model.n_order = [1 1 1];
flow_model.s_order = [1 1 1];
flow_model.sigma_order = [1 1 1];


%Solver options and regularization
flow_regs = regularization.Velocity.get_all_regs(mesh, B, xs, flow_model, opts, 'NoExpand', V);


% Bulk regularization parameter % possibly modify
lc = 0.0;
flow_regs(1).weight =  lc;
flow_regs(2).weight =  lc;
flow_regs(3).weight =  lc;
flow_regs(4).weight =  lc;
flow_regs(5).weight =  lc;


% Solve for the flow
flow_solv = LocationBasedVelocitySolver(mesh, B, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
flow_solv.rotation = xs.angle;
flow = flow_solv.get_solution(); % possibly modify
%figure
%spy(flow.M)
% Plot the state vector
flow.plot_solution()

%% Post-Processing - focus on decomposition of the solution
addpath(genpath(strcat(RF,'git\adcptools\post_processing')))
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
name = 'flow';
sav = 0;
animate_solution(u{1}, X, name, sav)

[u_decomp, u_avg] = D.decompose_function(u{1}); % U-Flow

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
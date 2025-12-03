%% complete script
clear all
clc

%% Path management
RF = 'C:\Users\arian\Documents\internship'; %RootFolder
addpath(genpath(strcat(RF,'/git/adcptools'))); %path to ADCPTools
addpath(genpath("C:\Users\arian\Documents\internship\git\adcptools\post_processing")) %% add postprocessing analysis
dir = 'C:\Users\Ariane.VandePas\Documents\results\Portneuf';
cd dir

%% load data
dat = open("C:\Users\Ariane.VandePas\Documents\Donnees_validation\2009\ADCP 2009\Portneuf_0\Portneuf.mat");
wl = open("C:\Users\Ariane.VandePas\Documents\Donnees_validation\2009\marégraphes_h_2009_HNE_NMM_3min\marégraphes_h_2009_HNE_NMM_3min\3300Portneuf2009_HNE_NMM_3min.mat");

%% Information for saving figures

transect = 'Portneuf';

%% adjustable input values

% starting regularization values
reg_we = [1,1,1,1,1];

% regularization method
regmet = 'track';

% mesh sizes
hor = 50; %[5,5, 15, 25, 25, 50, 100, 100, 150];
ver =  10; %[1,2.5, 2.5, 2.5, 5, 5, 5, 10, 10];

% temporary
track = 6;

% FINDING LAMBDA

% Curvature method
scalingfactorcur = 10000;
maxregcur = 0.0001;
regstepcur = 100;
pcur = [0.8,1,0.999]; %smoothing paramethers for smoothing spline fit (loglogrho, second rho, eta)

%Decay method
scalingfactordec = 10000000;
maxregdec = 0.000000015;
regstepdec = 100;
pdec = 0.9999;

% minimum error
max1 = 100000;
max2 = 100000;
step1 = 25;
step2 = 25;

%% Constituents

constituents = {'M2', 'M4'};

% waterlevel
filt = ~isnan(wl.h);
water_level = VaryingWaterLevel(datetime(wl.t(filt), 'ConvertFrom', 'datenum'), wl.h(filt));
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();

% Modify the following code to analyze the data

V = rdi.VMADCP(dat.dat);
% V.horizontal_position_provider = HorizontalPositionFromBottomTracking; % possibly modify

 V.water_level_object = water_level;  % return

B = BathymetryScatteredPoints(V);

% delete points too close to the surface
Bfilt = find(B.known(2,:)>0);

B.interpolator.span = .001;
figure;
B.plot

V.filters = Filter;

[ef, xs] = cross_section_selector(V);

% Mesh for plotting (only as start)

mesh_makers = SigmaZetaMeshFromVMADCP(ef, xs, B, 'NoExpand', V);

% input preferred mresh size

hver = 5; % depth mesh cell in m
hhor = 25; %width mesh cell in m

%% select max & min in that order

figure
plot(V.horizontal_position(1,:))
title('identify maximum and minimum with mouseclick and push enter')
[~, x] = ginput;%%check this part pls
maxx = x(1,1);
minx = x(2,1);


plot(V.horizontal_position(2,:));
title('identify maximum and minimum with mouseclick and push enter')
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

%% multiple meshsizes

    r = zeros(length(track),length(hor));
    e = zeros(length(track),length(hor));
    l = zeros(length(track),length(hor));
    l2 = zeros(length(track),length(hor));

for i = 1:length(hor)
    %install mesh size
    hhor = hor(i);
    hver = ver(i);

    n = round(lengthriv/hhor);
    z = round(depthriv/hver);
    acthor = lengthriv/n;
    actver = depthriv/z;
    
    fprintf('Used horizontal mesh size is: %.2f\n', acthor);
    fprintf('Used vertical mesh size is: %.2f\n', actver);
    
    mesh = mesh_makers.get_mesh(resn = n, resz = z);

    % solve initial flow conditions

    flow = get_tidal_model(V, constituents, mesh, B, xs, ef, reg_we);

    % plot solution

    figure(222)
    flow.plot_solution();
    
    stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'initrun');
    saveas(figure(222), stringfig1)

    % find lambda with the decay function and maximum curvature
    %track = unique(V.fileid)

    for j = 1:length(track)
         j = track %remove after tests

         % curvature method
         flow.solver.opts.cv_mode = regmet;
         [rho, lambda,eta] = cross_validate_1D_track(flow, 0, maxregcur, regstepcur,j);
         rho = cell2mat(rho);
         eta = cell2mat(eta);
         lambda = lambda * scalingfactorcur;
         [r(j,i),e(j,i),l(j,i)] = l_corner(rho(:,1),eta(:,1),lambda(:,1),pcur(1),pcur(2),pcur(3),300);
         l(j,i) = l(j,i)/scalingfactorcur;

         %decay method
          [rho, lambda,eta] = cross_validate_1D_track(flow, 0, maxregdec,regstepdec,j);
          rho = cell2mat(rho);
          eta = cell2mat(eta);
          lambda = lambda*scalingfactordec;
          [l2(j,i)] = expon(lambda(:,1),rho(:,1),pdec);
          l2(j,i) = l2(j,i)/scalingfactordec;

          % optimal regularization
          [CV, rpc, rps] = flow.cross_validate_2D_track([0,0], [max1,max2], [step1,step2], j);
          
          contourf(helpers.symlog(rpc), helpers.symlog(rps), reshape([CV{:,1}]./CV{1,1}, [step1,step2]), 100)
          colorbar
          colormap(flipud(helpers.cmaps('velmap')))
          % clim([0,2]) % 0 - very good (too good to be true) 1
          xlabel('cont lambda (symlog10)')
          ylabel('smoothness lambda (symlog10)')
          title('2D cross-validation: generalization error')
    end
        

    
end












function flow = get_tidal_model(V, constituents, mesh, bathy, xs, ef, reg_weights)

    % Create tidal model (incl. regularisation)
    flow_model = TaylorTidalVelocityModel();
    flow_model.constituents = constituents;

    % Initiate model
    flow_model.s_order = [1,1,1];       % Calculated only first order derivatives in all directions
    flow_model.n_order = [1,1,1];
    flow_model.sigma_order = [1,1,1];
    % velocity_model.get_names();
    opts = SolverOptions();

    % Initialize regularization
    flow_regs = regularization.Velocity.get_all_regs(mesh, bathy, xs, flow_model, opts, 'NoExpand', V);

    flow_regs(1).weight =  reg_weights(1);
    flow_regs(2).weight =  reg_weights(2);
    flow_regs(3).weight =  reg_weights(3);
    flow_regs(4).weight =  reg_weights(4);
    flow_regs(5).weight =  reg_weights(5);

    % % Solve
    flow_solv = LocationBasedVelocitySolver(mesh, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs);
    flow_solv.rotation = xs.angle;

    flow = flow_solv.get_solution();
end

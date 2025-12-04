%% complete script
clear all
clc

%% Path management
RF = 'C:\Users\arian\Documents\internship'; %RootFolder
addpath(genpath(strcat(RF,'/git/adcptools'))); %path to ADCPTools
addpath(genpath("C:\Users\arian\Documents\internship\git\adcptools\post_processing")) %% add postprocessing analysis
dir = 'C:\Users\Ariane.VandePas\Documents\results\Portneuf';
cd(dir)

%% Path management ECCC laptop
% RF = 'C:\Users\arian\Documents\internship'; %RootFolder
addpath(genpath('C:\Users\Ariane.VandePas\Documents\GitHub\adcptools')); %path to ADCPTools
addpath(genpath("C:\Users\Ariane.VandePas\Documents\GitHub\adcptools\+post_processing")) %% add postprocessing analysis
addpath('C:\Users\Ariane.VandePas\Documents\GitHub\adcptools\+post_processing\cartesian')
addpath('C:\Users\Ariane.VandePas\Documents\GitHub\adcptools\+post_processing\plot')
 
dir = 'C:\Users\Ariane.VandePas\Documents\results\Portneuf';
cd(dir)

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
hor = 150; %[5,5, 15, 25, 25, 50, 100, 100, 150];
ver =  10; %[1,2.5, 2.5, 2.5, 5, 5, 5, 10, 10];

% temporary
track = 1;

% FINDING LAMBDA

% Curvature method
scalingfactorcur = 10000;
maxregcur = 0.0001;
regstepcur = 200;
pcur = [0.8,1,0.999]; %smoothing paramethers for smoothing spline fit (loglogrho, second rho, eta)

%Decay method
scalingfactordec = 10000000;
maxregdec = 0.0000001;
regstepdec = 200;
pdec = 0.99999;

% minimum error
max1 = 100000;
max2 = 100000;
step1 = 25;
step2 = 25;

% add continuity regulrization for curvature and decay methods
con = false;

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
cd(dir)
    r = zeros(length(track),length(hor));
    e = zeros(length(track),length(hor));
    l = zeros(length(track),length(hor));
    l2 = zeros(length(track),length(hor));
    l3 = cell(2,length(hor));

    L = zeros(1,length(hor));
    L2 = zeros(1,length(hor));
    L3 = cell(2,length(hor));

for i = 1:length(hor)
    
    %install mesh size
    hhor = hor(i);
    hver = ver(i);

    dirstr = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver));
    mkdir(dirstr);
    cd(dirstr) 

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

    figure(299)
    flow.plot_solution();
    
    stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'initrun');
    saveas(figure(299), stringfig1)
    saveas(figure(229), append(stringfig1,'.jpg'))


    % find lambda with the decay function and maximum curvature
    %track = unique(V.fileid)

    for j = 1:length(track)
         j = track %remove after tests

         % curvature method
         flow.solver.opts.cv_mode = regmet;
         [rho, lambda,eta] = cross_validate_1D_track(flow, 0, maxregcur, regstepcur,j,con);
         rho = cell2mat(rho);
         eta = cell2mat(eta);
         lambda = lambda * scalingfactorcur;
         [r(j,i),e(j,i),l(j,i)] = l_corner(rho(:,1),eta(:,1),lambda(:,3),pcur(1),pcur(2),pcur(3),300);
         l(j,i) = l(j,i)/scalingfactorcur;

         %save figures
         stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'loglogcurvature_track_',num2str(j));
         saveas(figure(302), stringfig1)
         saveas(figure(302), append(stringfig1,'.jpg'))

         stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'fitlogcurvature_track_',num2str(j));
         saveas(figure(304), stringfig1)
         saveas(figure(304), append(stringfig1,'.jpg'))

         stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'fitetacurvature_track_',num2str(j));
         saveas(figure(305), stringfig1)
         saveas(figure(305), append(stringfig1,'.jpg'))

         stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'maxcurvature_track_',num2str(j));
         saveas(figure(306), stringfig1)
         saveas(figure(306), append(stringfig1,'.jpg'))

         stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'etarho_track_',num2str(j));
         saveas(figure(307), stringfig1)
         saveas(figure(307), append(stringfig1,'.jpg'))

          %decay method
          [rho, lambda,eta] = cross_validate_1D_track(flow, 0, maxregdec,regstepdec,j,con);
          rho = cell2mat(rho);
          eta = cell2mat(eta);
          lambda = lambda*scalingfactordec;
          [l2(j,i)] = expon(lambda(:,3),rho(:,1),pdec);
          l2(j,i) = l2(j,i)/scalingfactordec;

          %save figure
          stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'decayrate_track_',num2str(j));
          saveas(figure(308), stringfig1)
          saveas(figure(308), append(stringfig1,'.jpg'))

          % optimal regularization
          [CV, rpc, rps] = flow.cross_validate_2D_track([0,0], [max1,max2], [step1,step2], j);
          
          figure(309)
          contourf(helpers.symlog(rpc), helpers.symlog(rps), reshape([CV{:,1}]./CV{1,1}, [step1,step2]), 100)
          colorbar
          colormap(flipud(helpers.cmaps('velmap')))
          % clim([0,2]) % 0 - very good (too good to be true) 1
          xlabel('cont lambda (symlog10)')
          ylabel('smoothness lambda (symlog10)')
          title('2D cross-validation: generalization error')

         stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), '2D_',num2str(j));
         saveas(figure(309), stringfig1)
         saveas(figure(309), append(stringfig1,'.jpg'))
    
         [val,idx] = min(cell2mat(CV(:,1))./cell2mat(CV(1,1)));
         smoparidx = rem(idx,step2);
         contparidx = (idx-smoparidx)/step1;

         l3{1,i}(j) = rpc(smoparidx,contparidx);
         l3{2,i}(j) = rps(smoparidx,contparidx);
    end

% calculate mean lambda values

L(1,i) = mean(l(:,i));
L2(1,i) = mean(l2(:,i));
L3{1,i} = mean(l3{1,i}(:));
L3{2,i} = mean(l3{1,i}(:));
    
% plot results

for k = 1:3 %once for each method

    if k == 1
        method = 'curvature';
        if con
            reg = [L(1,i), L(1,i), L(1,i), L(1,i), L(1,i)];
        else
            reg = [0,0, L(1,i), L(1,i), 0];
        end

    elseif k ==2
        method = 'decay';
        if con
            reg = [L2(1,i), L2(1,i), L2(1,i), L2(1,i), L2(1,i)];
        else
            reg = [0,0, L2(1,i), L2(1,i), 0];
        end

    elseif k == 3
        method = 'minerr';
        reg = [L3(1,i), L3(1,i), L3(2,i), L3(2,i), L3(1,i)];
    end

% get model
flow = get_tidal_model(V, constituents, mesh, B, xs, ef, reg);

%preparation for analyzing individual cells
[pars_U, pars_V, pars_W] = sort_flow_output(flow); % define to get model for single direction

t0 = (datenum(V.time(1)))*86400; % determine time for plotting single measurements 
t_end = (datenum(V.time(end)))*86400;  
t_plot = (t0:10:t_end);

% prepare plotting results
tim = flow.solver.adcp.time; % get time for plotting
Tlim(1)= min(tim);
M2T = flow.solver.model.periods(1,1)/(3600*24);
Tlim(2) = Tlim(1) + M2T;
Tlimn = datenum(Tlim);
limu = {[Tlimn(1), Tlimn(2)];... %days
        [min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5];...
        [0,1]};
    
% evaluation resolution
tres = 60;
evres = [tres, flow.solver.mesh.nverticals, flow.solver.mesh.max_ncells_vertical]; % t, y , sigma
X = get_coords(limu, evres);
    
reg_idx = 1; % only relevant if multiple regularization parameter settings are entered upon model fitting.
u = get_var(flow, X, reg_idx); % Vector variable on regular sigma grid. Three cells are the three velocity components.
    
[H, Wl, Zb] = get_H(X, flow, 0);
    
    if size(H, 3) == 1
        H = repmat(H, [1,1,numel(X.sig)]);
    end
X.Z = Zb + X.Sig.*H;
    
D = post_processing.Decomposition(X = X, H = H, wl = Wl(:,1), zb = Zb(1,:)');
    

% analyze 3 velocity directions

    for d = 1:3
    if d == 1
        dir = 'U';
        pars = pars_U;
    elseif d == 2
        dir = 'V';
        pars = pars_V;
    else
        dir = 'W';
        pars = pars_W;
    end
    
    
    % RMSE on transect
    [RMSE_cell] = plot_mrse_mesh(mesh, dir, pars,t_plot,constituents,xs,V);
    
    stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'method__',method,'RMSE_dir_',dir);
    saveas(figure(308 + 2*d), stringfig1)
    saveas(figure(308 + 2*d), append(stringfig1,'.jpg'))
    
    % model fits for evenly spaced cells
    plot_ts_random_cells(mesh, dir, pars, t_plot, constituents, xs, V);
    
    stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'method__',method,'cells_dir_',dir);
    saveas(figure(309 + 2*d), stringfig1)
    saveas(figure(309 + 2*d), append(stringfig1,'.jpg'))
    
    % get results for each direction
    % Plot some variables
    name = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'method__',method,'dir_',dir,'results.gif');
    sav = 1;
    animate_solution(u{d}, X, name, sav)
    
    [decomp, avg] = D.decompose_function(u{d}); % U-Flow
    
    figure(315 + d)
    D.plot_components(decomp, 'velmap')
    stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'method__',method,'_dir_',dir, '_resultsdecom');
    saveas(figure(315+d), stringfig1)
    saveas(figure(315+d), append(stringfig1,'.jpg'))

    figure(318 + d)
    D.plot_components(avg, 'velmap')
    stringfig1 = append(transect, '_hor_', num2str(hhor), '_ver_', num2str(hver), 'method__',method,'_dir_',dir, '_resultsdecomavg');
    saveas(figure(318+d), append(stringfig1,'.jpg'))
    end

end

end



%% functions


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

function [pars_U, pars_V, pars_W] = sort_flow_output(flow)
    pars = flow.pars;
    pars_names = flow.solver.model.all_names;
    pars_U_names = {};
    pars_U_columns = [];
    pars_V_names = {};
    pars_V_columns = [];
    pars_W_names = {};
    pars_W_columns = [];

    for i = 1:length(pars_names)
        if contains(pars_names{i}, 'u', 'IgnoreCase', true)
            pars_U_names{end+1} = pars_names{i};  % Append item to uItems
            pars_U_columns(end+1) = i;        % Store the index
        elseif contains(pars_names{i}, 'v', 'IgnoreCase', true)
            pars_V_names{end+1} = pars_names{i};  % Append item to vItems
            pars_V_columns(end+1) = i;        % Store the index
        elseif contains(pars_names{i}, 'w', 'IgnoreCase', true)
            pars_W_names{end+1} = pars_names{i};  % Append item to wItems
            pars_W_columns(end+1) = i;        % Store the index
        end
    end
    pars_U = pars(:,pars_U_columns);
    pars_V = pars(:,pars_V_columns);
    pars_W = pars(:,pars_W_columns);
end

function plot_ts_random_cells(mesh, lett ,pars_U,t_plot,constituents,xs,V)
    % Tidal components' periods (in hours)
    M2_period = 12.4206012;          % Semi-diurnal component
    M4_period = 6.210300601;         % M4 (fourth diurnal component)
    M6_period = 4.140200401;
    M1_period = 24.84120241;
    M3_period = 8.280400802;


    oM2 = 1/(M2_period*3600)*2*pi;       % rad/s
    oM4 = 1/(M4_period*3600)*2*pi;       % rad/s
    oM6 = 1/(M6_period*3600)*2*pi;       % rad/s
    oM1 = 1/(M1_period*3600)*2*pi;       % rad/s
    oM3 = 1/(M3_period*3600)*2*pi;       % rad/s


    random_cells = round(linspace(20,mesh.ncells-20,8));
   
    for i = 1:length(random_cells)
        cell = random_cells(i);
        if mesh.domains(cell,1) == 4 || mesh.domains(cell,1) == 3 || mesh.domains(cell,1) == 2
            cell = cell + 1;
        end
        random_cells(i) = cell;
    end
    
    if lett == 'U'
        figure(311);
    elseif lett == 'V'
        figure(313)
    else
        figure(315)
    end

    for k = 1:8
        CellID = random_cells(k);
        subplot(4, 4, k);

        calc_vel_U = pars_U(CellID,1);
        for constituent = 1:length(constituents)
            cur_constiturent = constituents{constituent};
            pars_col = constituent*2;
            period = eval(['o', cur_constiturent]);
            calc_vel_U = calc_vel_U + pars_U(CellID,pars_col)*cos(period.*t_plot) + pars_U(CellID,pars_col+1)*sin(period.*t_plot);
        end

        plot(datetime(t_plot/86400, 'ConvertFrom', 'datenum'), calc_vel_U);

        [time_in_column, vel_in_cell] = extract_measured_velocity_mesh_cell(CellID,V,xs,mesh, lett);
        hold on

        plot(time_in_column,vel_in_cell,'k.','MarkerSize',4)

        title(strcat("CellID = ",string(CellID)))
    end
    subplot(2, 1, 2);
    patch(mesh.n_patch,mesh.z_patch,'white')
    hold on
    patch(mesh.n_patch(:,random_cells),mesh.z_patch(:,random_cells),'r')

    plot(mesh.nb_all,mesh.zb_all,'k','Linewidth',2);
    plot(mesh.nw,mesh.nw*0+mesh.water_level,'b','Linewidth',2);
end

function [RMSE_cell] = plot_mrse_mesh(mesh, lett, pars_U,t_plot,constituents,xs,V)

    for j = 1:mesh.ncells
     % Tidal components' periods (in hours)
    M2_period = 12.4206012;          % Semi-diurnal component
    M4_period = 6.210300601;         % M4 (fourth diurnal component)
    M6_period = 4.140200401;
    M1_period = 24.84120241;
    M3_period = 8.280400802;


    oM2 = 1/(M2_period*3600)*2*pi;       % rad/s
    oM4 = 1/(M4_period*3600)*2*pi;       % rad/s
    oM6 = 1/(M6_period*3600)*2*pi;       % rad/s
    oM1 = 1/(M1_period*3600)*2*pi;       % rad/s
    oM3 = 1/(M3_period*3600)*2*pi;       % rad/s


        CellID = j;

        calc_vel_U = pars_U(CellID,1);
        for constituent = 1:length(constituents)
            cur_constiturent = constituents{constituent};
            pars_col = constituent*2;
            period = eval(['o', cur_constiturent]);
            calc_vel_U = calc_vel_U + pars_U(CellID,pars_col)*cos(period.*t_plot) + pars_U(CellID,pars_col+1)*sin(period.*t_plot);
        end

         [time_in_column, vel_in_cell] = extract_measured_velocity_mesh_cell(CellID,V,xs,mesh, lett);
         
         RMSE = []; 
         for h = 1:length(time_in_column)
         [val,idx] = min(abs(datetime(t_plot/86400, 'ConvertFrom', 'datenum') - time_in_column(h)));

         u_vel = (calc_vel_U(1,idx));
         RMSE(:,h) = ((vel_in_cell(:,h)-u_vel).^2);
         end
         RMSE_cell(j) = sqrt(mean(RMSE, 'all', 'omitnan'));

         fprintf('RMSE percentage: %2.2f percent \n', (j/mesh.ncells)*100)
    end
 
    if lett == 'U'
        figure(310);
    elseif lett == 'V'
        figure(312)
    else
        figure(314)
    end

    patch(mesh.n_patch,mesh.z_patch,RMSE_cell)

    cMap = interp1([0;1],[0 1 0; 1 0 0],linspace(0,1,256));
    colormap(cMap)
    colorbar
    clim([0 0.4])
    hold on
    plot(mesh.nb_all,mesh.zb_all,'k','Linewidth',2);
    plot(mesh.nw,mesh.nw*0+mesh.water_level,'b','Linewidth',2);
    hold off

end

function [time_in_column, vel_in_cell] = extract_measured_velocity_mesh_cell(cellID,V,xs,mesh, direction)
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
        end
    end

    bounds_cellID_x = [min(mesh.x_patch(:,cellID))-1 max(mesh.x_patch(:,cellID))+1];
    bounds_cellID_y = [min(mesh.y_patch(:,cellID)) max(mesh.y_patch(:,cellID))];
    bounds_cellID_z = [min(mesh.z_patch(:,cellID))-0.25 max(mesh.z_patch(:,cellID))+0.25];


    valid_column = V.horizontal_position(2,:) >= bounds_cellID_y(1) & V.horizontal_position(2,:) <= bounds_cellID_y(2);

    depth_in_column = vel_pos{1}(:,(valid_column));

    time_in_column = tim(:,valid_column);
    if isequal(direction, 'U')
        direct = 1;
    elseif isequal(direction, 'V')
        direct = 2;
    elseif isequal(direction, 'W')
        direct = 3;
    else
        warning('Direction not correct')
    end
    vel_in_column = vel_sn{1}(:,(valid_column),direct);
    vel_in_cell = NaN(size(vel_in_column));
    for col = 1:size(vel_in_column, 2)
        valid_rows = (depth_in_column(:, col) >= bounds_cellID_z(1)) & (depth_in_column(:, col) <= bounds_cellID_z(2));
        % Extract the velocities for the valid rows in this column
        vel_in_cell(valid_rows, col) = vel_in_column(valid_rows, col);
    end
end
% Script to play for Ariane
close all
clc
clear all
%%
dir = 'C:\Users\arian\Documents\internship\datafiles matlab\portneuf';
addpath(genpath(dir))
load('portneuf2009.mat')
%% Path management
RF = 'C:\Users\arian\Documents\internship'; %RootFolder
addpath(genpath('C:\Users\arian\Documents\internship\git\adcptools')); %path to ADCPTools of Bart Vermeulen
% addpath(genpath(strcat(RF,'Tools\adcptools'))); %possible other folders
%% Quick documentation walkthrough - comment out

%open_adcptools_documentation()

%% Constituents

constituents = {'M2', 'M4'};

%% Loading in the data 
% %test for quebec script
addpath('./Donnees_validation'); %path to data
%dat = rdi.readDeployment('rijn', './data');
%dat = rdi.readDeployment('Lauzon_0_0', './data/Lauzon_0');

% addpath('./data'); %path to data
%dat = rdi.readDeployment('rijn', './data');
dat = rdi.readDeployment('Portneuf_0_0', './2009/ADCP 2009/Portneuf_0');
%% Load water level data
load('C:\Users\arian\Documents\internship\Donnees_validation\2009\marégraphes_h_2009_HNE_NMM_3min\marégraphes_h_2009_HNE_NMM_3min\3300Portneuf2009_HNE_NMM_3min.mat')


%% waterlevel
filt = ~isnan(h);
water_level = VaryingWaterLevel(datetime(t(filt), 'ConvertFrom', 'datenum'), h(filt));
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();

%% Modify the following code to analyze the data

V = rdi.VMADCP(dat);
%  V.horizontal_position_provider = HorizontalPositionFromBottomTracking; % possibly modify

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

mesh_makers = SigmaZetaMeshFromVMADCP(ef, xs, B, 'NoExpand', V);

% input preferred mresh size

hver = 2.5; % depth mesh cell in m
hhor = 15; %width mesh cell in m

%% select max & min in that order

newplot
plot(V.horizontal_position(1,:))
[~, x] = ginput;
maxx = x(1,1);
minx = x(2,1);

newplot
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

%% End 
constituents = {'M2', 'M4'};
tide = '';
channel = 'portneuf';
transect = '2009';



 reg_weights = [1,1,1,1,1];  % [0,0,0,0,0], [0.25,0.25,0.25,0.25,0.25], [1,1,1,1,1], [1,1,10,10,1][10,10,100,100,10]


%% saving results
model_name = get_model_name(channel,tide,transect,constituents,reg_weights);
% model_name = strcat('CoarseGrid_HighReg_', model_name);



%% plot mesh
fig = figure;
fig.Units = 'centimeters';
fig.Position = [0.5, 0.5, 10, 7];
mesh.plot()
ylabel('Height [m +PWD]')
xlabel('Location along cross-section [m]')
legend('River Bed','Median Water Level', Location='southeast')
%% Create tidal model (incl. regularisation)

fit_new = true;
if fit_new
    flow = get_tidal_model(V, constituents, mesh, B, xs, ef, reg_weights, 1, model_name);
else
    flow = load_tidal_model(model_name);
end



%% Sort output

[pars_U, pars_V, pars_W] = sort_flow_output(flow);

%% read single tracks
fit_new = true;
if fit_new
    flow_tracks = get_model_individual_track(V, mesh, B, xs, reg_weights, 0, model_name);
else
    flow_tracks = load_model_individual_track(model_name);
end
%% Time settings

t0 = (datenum(V.time(1)))*86400;
t_end = (datenum(V.time(end)))*86400;  
t_plot = (t0:300:t_end);

%% plot results - changes in selected cell

plot_ts_random_cells(mesh, 'W', 3, pars_W, t_plot, constituents, xs, V, channel, flow_tracks, 0, model_name, 31,100,183,277,405,517,650,736);

%% plot results - mesh video

 flow_pattern_video(flow, mesh, constituents, t_plot,  1, model_name)

%% a, b to Amplitude and phase

create_amplitude_bar_chart(constituents, pars_U, pars_V, pars_W,  1, model_name)

%% CV MSE
%Heavy function takes a long time to run!
if include_CV
    CV = CV_model_performance(flow, reg_weights, 1, model_name);
end
%% Fixed Point

% time shift based on t = s/v not included
[point_loc_x, point_loc_y, water_velocity_fixed] = compare_transect_to_fixed(channel, tide, transect, constituents, fit_new, V, mesh, xs, pars_U, pars_V, pars_W, t_plot, flow_tracks, 1, model_name);
%% Create areal orientation figure

plot_areal_image(V, xs, channel, point_loc_x, point_loc_y, 1, model_name)


%% Fixed Point select

[point_loc_x, point_loc_y, water_velocity_fixed] = compare_transect_to_fixed_simp(channel, tide, transect, fixed_point, constituents, fit_new, V, mesh, xs, pars_U, pars_V, pars_W, t_plot, flow_tracks, 1, model_name);




%% Functions
function model_name = get_model_name(channel,tide,transect,constituents,reg_weights)
    model_constituents = strjoin(constituents, '_');
    if all(reg_weights == 0)
        model_reg = 'no_reg';  % If all numbers are 0
    elseif isequal(reg_weights, [1,1,10,10,1])
        model_reg = 'high_smoothing';
    elseif isequal(reg_weights, [1,1,1,1,1])
        model_reg = 'reg_1';
    elseif isequal(reg_weights, [0.25,0.25,0.25,0.25,0.25])
        model_reg = 'reg_025';
    else
        model_reg = 'opt_reg';  % If not all numbers are 0
    end
    name_base = strcat(channel,'_',tide,'_',transect);
    model_name = strcat(model_constituents, '_',model_reg, '_', name_base);
end


function reg_weights = load_opt_reg_weights(channel,transect,tide,constituents)
    if isequal(channel,'Tetulia') & isequal(transect,'T2')
        fixed_point = {'V2'};
    elseif isequal(channel,'Tetulia') & isequal(transect,'T3')
        fixed_point = {'V1'};
    elseif isequal(channel,'Meghna')
        fixed_point = {'V1', 'V3'} ;
    end
    for f = 1:length(fixed_point)
        fixed_point_f = fixed_point{f};
        model_constituents = strjoin(constituents, '_');
        name_base = strcat(channel,'_',tide,'_',transect);
        file_name = strcat(model_constituents, '_', name_base, '_opt_weights');
        file_name = strcat(file_name,'_', fixed_point_f);
        load_filename  = strcat('Matlab_data\Optimised_fit\', file_name, '_grid_search.mat');
        load(load_filename)
        if f == 1
            grid_search_rmse_1 = grid_search_rmse;
        else
            grid_search_rmse_1 = grid_search_rmse_1 + grid_search_rmse;
        end
    end
    grid_search_rmse = grid_search_rmse_1;

    [minRowVals, rowIndices] = min(grid_search_rmse);        % min of each column, get row indices
    [minVal, col] = min(minRowVals)          % find overall min and its column
    row = rowIndices(col);                    % get corresponding row

    nreg = 10;
    ls = logspace(-4,5,nreg)'; % smoothing parameter is varied in the k-fold CV

    reg_weights = [ls(col), ls(col), ls(row), ls(row), ls(col)]
    % fprintf('Lowest value is %.2f at row %d, column %d\n', minVal, row, col);
end



function flow = get_tidal_model(V, constituents, mesh, bathy, xs, ef, reg_weights, save_parameters, model_name)

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

function flow_tracks = get_model_individual_track(V, mesh, bathy, xs, reg_weights, save_parameters, model_name)
    flow_model = VelocityModel;
    % no constituents as this is a very simple model (average in u, v, w)

    opts = SolverOptions(extrapolate_vert = 0, lat_weight_factor = 1);

    V_ef_tracks = struct();
    flow_tracks = struct();

    for i=1:length(unique(V.fileid))
        disp(strcat(string(i), '/', string(length(unique(V.fileid)))))
        track_name = strcat('Track_',string(i));

        V_ef_track_name = strcat('V_ef_Track_',string(i));
        V_ef_tracks.(V_ef_track_name) = EnsembleFilter(V, V.fileid ~= i);

        V_ef_tracks_i = V_ef_tracks.(V_ef_track_name);

        flow_regs = regularization.Velocity.get_all_regs(mesh, bathy, xs, flow_model, opts, 'NoExpand', V_ef_tracks_i);

        flow_regs(1).weight =  reg_weights(1);
        flow_regs(2).weight =  reg_weights(2);
        flow_regs(3).weight =  reg_weights(3);
        flow_regs(4).weight =  reg_weights(4);
        flow_regs(5).weight =  reg_weights(5);

        flow_solv = LocationBasedVelocitySolver(mesh, bathy, xs, V_ef_tracks_i, flow_model, opts, 'NoExpand', V, flow_regs);

        flow_solv.rotation = xs.angle;

        flow_tracks.(track_name) = flow_solv.get_solution();

        if save_parameters
            save_name = strcat("Matlab_data\Fitted_parameters\Individual_track_fit\",model_name,'_pars_track.mat');
            save(save_name, "flow_tracks");
        end
    end
end

function plot_ts_random_cells(mesh, lett, nr ,pars_U,t_plot,constituents,xs,V, channel, flow_tracks, save_figure, model_name,a,b,c, d, e, f, g, h)
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


    x_cell = 40;
    y_cell = 15;
    CellID = 157 ; %find(mesh.col_to_cell == x_cell & mesh.row_to_cell == y_cell);
    % random_cells = sort(randperm(mesh.ncells, 8), "ascend");
    random_cells = round(linspace(20,mesh.ncells-20,8));
   
    for i = 1:length(random_cells)
        cell = random_cells(i);
        if mesh.domains(cell,1) == 4 || mesh.domains(cell,1) == 3 || mesh.domains(cell,1) == 2
            cell = cell + 1;
        end
        random_cells(i) = cell;
    end

    random_cells = [a,b,c, d, e, f, g, h];

    screenSize = get(0, 'ScreenSize');
    figure('Position', [1, 1, screenSize(3)-100, screenSize(4)-100]);
    for k = 1:8
        CellID = random_cells(k);

        flows_single_in_cell = table();
        time_single = table();

        for i = 1:numel(fieldnames(flow_tracks))
            measurement_time = median(V.time(V.fileid == i));

            track_name = strcat("Track_", num2str(i));
            newRow = flow_tracks.(track_name).pars(CellID,:);

            if isempty(flows_single_in_cell)
                % If NewTable is empty, just initialize it with the first newRow
                flows_single_in_cell = newRow;
                time_single = measurement_time;
            else
                % Concatenate the newRow to the NewTable
                flows_single_in_cell = [flows_single_in_cell; newRow];
                time_single = [time_single; measurement_time];
            end
        end

        subplot(4, 4, k);

        calc_vel_U = pars_U(CellID,1);
        for constituent = 1:length(constituents)
            cur_constiturent = constituents{constituent};
            pars_col = constituent*2;
            period = eval(['o', cur_constiturent]);
            calc_vel_U = calc_vel_U + pars_U(CellID,pars_col)*cos(period.*t_plot) + pars_U(CellID,pars_col+1)*sin(period.*t_plot);
        end

        plot(datetime(t_plot/86400, 'ConvertFrom', 'datenum'), calc_vel_U);

        [time_in_column, vel_in_cell] = extract_measured_velocity_mesh_cell(CellID,V,xs,mesh,channel, lett);
        hold on

         plot(time_in_column,vel_in_cell,'k.','MarkerSize',4)

        hold on
         scatter(time_single(:,1), flows_single_in_cell(:,nr))

        title(strcat("CellID = ",string(CellID)))
    end
    subplot(2, 1, 2);
    patch(mesh.n_patch,mesh.z_patch,'white')
    hold on
    patch(mesh.n_patch(:,random_cells),mesh.z_patch(:,random_cells),'r')

    hbed = plot(mesh.nb_all,mesh.zb_all,'k','Linewidth',2);
    hwater = plot(mesh.nw,mesh.nw*0+mesh.water_level,'b','Linewidth',2);

    if save_figure
        fig_location = strcat("Matlab_data\Figures_models\Timeseries_indivual_cells\",model_name, "_timeseries.png");
        saveas(gcf,fig_location)
    end
end

function flow = load_tidal_model(model_name)
    load_name = strcat("Matlab_data\Fitted_parameters\",model_name,'_pars.mat');
    load(load_name)
end

function flow_tracks = load_model_individual_track(model_name)
    load_name = strcat("Matlab_data\Fitted_parameters\Individual_track_fit\",model_name,'_pars_track.mat');
    load(load_name)
end

function [time_in_column, vel_in_cell] = extract_measured_velocity_mesh_cell(cellID,V,xs,mesh,channel, direction)
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

    bounds_cellID_x = [min(mesh.x_patch(:,cellID))-1 max(mesh.x_patch(:,cellID))+1];
    bounds_cellID_y = [min(mesh.y_patch(:,cellID)) max(mesh.y_patch(:,cellID))];
    bounds_cellID_z = [min(mesh.z_patch(:,cellID))-0.25 max(mesh.z_patch(:,cellID))+0.25];

    if isequal(channel, 'Meghna')
        valid_column = V.horizontal_position(1,:) >= bounds_cellID_x(1) & V.horizontal_position(1,:) <= bounds_cellID_x(2);
    else
        valid_column = V.horizontal_position(2,:) >= bounds_cellID_y(1) & V.horizontal_position(2,:) <= bounds_cellID_y(2);
    end
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

function [velocity_prime, velocity_second, alfa] = rotate_sn2ps(velocity_U, velocity_V)
    mean_vel_U = mean(velocity_U, "omitnan");
    mean_vel_V = mean(velocity_V, "omitnan");
    alfa = atan2(mean_vel_V,mean_vel_U);


    theta = -alfa;
    R = [cos(theta), -sin(theta); sin(theta),  cos(theta)];
    velocity_prime = R(1,1) * velocity_U + R(1,2) * velocity_V;
    velocity_second = R(2,1) * velocity_U + R(2,2) * velocity_V;
end

function flow_pattern_video(flow, mesh, constituents, t_plot,  save_video, model_name)
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



    screenSize = get(0, 'ScreenSize');
    figure('Position', [1, 1, screenSize(3)-100, screenSize(4)-100]);
%     save_name = strcat("Matlab_data\Figures_models\Video_crosssection_primary\",model_name,"_flow_pattern");
    video = VideoWriter('C:\Users\arian\Documents\internship', 'MPEG-4');
    video.FrameRate = 5;  % Adjust frame rate as needed
    open(video);  % Open video file for writing



    [pars_U, pars_V, pars_W] = sort_flow_output(flow);

    for n = 1:length(t_plot)
        clf; % clear figure
        t_n = t_plot(n);
        velocity_U = pars_U(:,1);
        velocity_V = pars_V(:,1);
        velocity_W = pars_W(:,1);
        for constituent = 1:length(constituents)
            cur_constiturent = constituents{constituent};
            pars_col = constituent*2;
            period = eval(['o', cur_constiturent]);
            velocity_U = velocity_U + pars_U(:,pars_col)*cos(period.*t_n) + pars_U(:,pars_col+1)*sin(period.*t_n);
            velocity_V = velocity_V + pars_V(:,pars_col)*cos(period.*t_n) + pars_V(:,pars_col+1)*sin(period.*t_n);
            velocity_W = velocity_W + pars_W(:,pars_col)*cos(period.*t_n) + pars_W(:,pars_col+1)*sin(period.*t_n);
        end


        [velocity_prime, velocity_second, alfa] = rotate_sn2ps(velocity_U, velocity_V);
        % if ismember(n, [1,50,100])
        %     summary(velocity_second)
        %     summary(velocity_U)
        %     summary(velocity_V)
        %     summary(velocity_W)
        % end
        patch(mesh.n_patch,mesh.z_patch,velocity_prime, 'EdgeColor', 'none')
        hold on

        hbed = plot(mesh.nb_all,mesh.zb_all,'k','Linewidth',2);
        hwater = plot(mesh.nw,mesh.nw*0+mesh.water_level,'b','Linewidth',2);

        xlim([min(mesh.n_middle)-300, max(mesh.n_middle)+300])
        ylim([floor(min(mesh.zb_all)), ceil(mesh.water_level)])

        %add arrows
        x_limits = xlim;
        y_limits = ylim;
        dx = diff(x_limits);
        dy = diff(y_limits);

        aspect_ratio_fig = 1;% dy / dx;
        Fv_adjusted = velocity_second/ aspect_ratio_fig;
        Fw_adjusted = velocity_W;

        % Resulting velocity average
        average_sec_velocity = mean(sqrt(velocity_second.^2+velocity_W.^2),"omitnan");

        % threshold = prctile(sqrt(velocity_second.^2+velocity_W.^2),95);
        threshold = inf;
        outlier_indices = abs(Fv_adjusted*aspect_ratio_fig) > threshold | abs(Fw_adjusted) > threshold;

        Fv_adjusted(outlier_indices) = 0;
        Fw_adjusted(outlier_indices) = 0;

        scale_factor = 1 % 0.5/average_sec_velocity;
        q = quiver(mesh.n_middle(mesh.col_to_cell)',...
            mesh.z_center, Fv_adjusted*scale_factor, Fw_adjusted*scale_factor,...
            'Color','k');

        q.ShowArrowHead = 'off';
        q.Marker = '.';
        q.AutoScale = 'on';
        % q.AutoScaleFactor = 4;


        % x = mesh.n_middle(mesh.col_to_cell)';
        % y = mesh.z_center;
        % outliers = zeros(size(x(outlier_indices)));
        % q_outlier = quiver(x(outlier_indices), ...
        %     y(outlier_indices), outliers, outliers, ...
        %     'Color', 'r');
        % q_outlier.ShowArrowHead = 'off';
        % q_outlier.Marker = '.';
        % q_outlier.AutoScale = 'off';
        %
        % % Add legend
        % legend_value = ceil(average_sec_velocity/0.5)*0.5;
        % xpos = x_limits(2)-dx/8;
        % ypos = y_limits(1)+dy/10;
        % text_display = strcat(num2str(legend_value), ' m/s');
        % text(xpos,ypos,text_display,'HorizontalAlignment','left',...
        %     VerticalAlignment='top',BackgroundColor='w')
        % q_legend = quiver(xpos,ypos,(legend_value/ aspect_ratio_fig)*scale_factor,0, 'r');
        % q_legend.ShowArrowHead = 'off';
        % q_legend.Marker = '.';
        % q_legend.AutoScale = 'off';

        colorbar;
        cmap = cmocean('Balance');
        trim_amount = round(0.1 * size(cmap, 1));
        new_cmap = cmap(trim_amount+1:end-trim_amount, :);
        colormap(new_cmap)
        clim([-2 2]);
        title_str = strrep(model_name, '_', '\_');
        title(title_str)
        subtitle(strcat(string(datetime(t_n/86400, 'ConvertFrom', 'datenum')),' Angle with mesh: ', num2str(rad2deg(alfa)), '°'))

        angle_xs = rad2deg(atan2(mesh.xs.direction_orthogonal(1),mesh.xs.direction_orthogonal(2)));
        angle_xs_n = deg2rad(90 - rad2deg(atan2(mesh.xs.direction(1),mesh.xs.direction(2))));
        rotated_direction = rad2deg(alfa) + angle_xs;
        % Convert angle to radians (from x-axis)
        angle_rad = deg2rad(90 - rotated_direction);


        % Arrow start (normalized coordinates)
        x_start = 0.7;
        y_start = 0.95;
        arrow_length = 0.05;  % relative to figure size

        x_end = x_start + arrow_length * cos(angle_rad);
        y_end = y_start + arrow_length * sin(angle_rad);
        annotation('arrow', [x_start x_end], [y_start y_end], ...
            'Color', 'b', 'LineWidth', 4);

        line_length = 0.05;
        x_end = x_start + line_length * cos(angle_xs_n);
        y_end = y_start + line_length * sin(angle_xs_n);
        annotation('line', [x_start x_end], [y_start y_end], ...
            'Color', 'k', 'LineWidth', 2);
        x_end = x_start - line_length * cos(angle_xs_n);
        y_end = y_start - line_length * sin(angle_xs_n);

        annotation('line', [x_start x_end], [y_start y_end], ...
            'Color', 'k', 'LineWidth', 2);

        frame = getframe(gcf);

        writeVideo(video, frame);


    end
    if save_video
        close(video);
    end
end

function create_amplitude_bar_chart(constituents, pars_U, pars_V, pars_W, save_figure, model_name)
    A_phi_U = pars_U(:,1);
    A_phi_V = pars_V(:,1);
    A_phi_W = pars_W(:,1);
    for constituent = 1:length(constituents)
        cur_constiturent = constituents{constituent};
        pars_col = constituent*2;
        a_U = pars_U(:,pars_col);
        b_U = pars_U(:,pars_col+1);
        A_phi_U = [A_phi_U, sqrt(a_U.^2 + b_U.^2), atan(b_U ./ a_U)];
        a_V = pars_V(:,pars_col);
        b_V = pars_V(:,pars_col+1);
        A_phi_V = [A_phi_V, sqrt(a_V.^2 + b_V.^2), atan(b_V ./ a_V)];
        a_W = pars_W(:,pars_col);
        b_W = pars_W(:,pars_col+1);
        A_phi_W = [A_phi_W, sqrt(a_W.^2 + b_W.^2), atan(b_W ./ a_W)];
    end

    screenSize = get(0, 'ScreenSize');
    figure('Position', [1, 1, screenSize(3)-100, screenSize(4)-100]);
    hold on
    col =  2 * (1:length(constituents));
    subplot(3, 1, 1);
    plot(A_phi_U(:,col), '.')
    legend()
     bar(mean(A_phi_U(:,col),"omitnan"));
    title('Amplitudes U')
    xticklabels(constituents)

    subplot(3, 1, 2);
    bar(mean(A_phi_V(:,col),"omitnan"));
    title('Amplitudes V')
    xticklabels(constituents)

    subplot(3, 1, 3);
    bar(mean(A_phi_W(:,col),"omitnan"));
    title('Amplitudes W')
    xticklabels(constituents)

%     if save_figure
%         fig_location = strcat("Matlab_data\Figures_models\Constituent_amplitude\",model_name, "_ConstituentA.png");
%         saveas(gcf,fig_location)
%     end
end

function CV = CV_model_performance(flow, reg_weights, save_MSE, model_name)
    % Heavy function to run
    flow.opts.cv_iter = 100;
    if any(reg_weights == 0) % functions do not work with weight=0 if not all zero.
        smallValue = 1e-15;
        for i = 1:length(reg_weights)
            if reg_weights(i) == 0
                reg_weights(i) = smallValue;
            end
        end
    end

    CV = flow.cross_validate_single(reg_weights);

    gen_error = sqrt([CV{:,1}]); %generalization error / error on validation set
    train_error = sqrt([CV{:,2}]); %training error
    fprintf('Cross-validation generalization error RMSE: %d m/s\n', gen_error)
    fprintf('Cross-validation training error RMSE: %d m/s\n', train_error)

    if save_MSE
        save_name = strcat("Matlab_data\Fitted_parameters\Model_performance\",model_name,'_CV_errors.mat');
        save(save_name, "CV");
    end
end

function [water_velocity_fixed, northings_fixed, eastings_fixed] = get_fixed_point_excel(channel, tide, fixed_point, save_fixed_point)
    folder_adcp_fixed = 'C:\Users\joris\OneDrive - Wageningen University & Research\Documenten\Wageningen\Master\Thesis\Data\Fixed Boat Velocity_Sentinel V\';
    % Read data pre-processed .xlsx
    if isequal(channel,'Tetulia') & isequal(tide,'Neap') & isequal(fixed_point, 'V1')
        file_name = (strcat(folder_adcp_fixed, 'Process Data\Neap\Tentulia River\26092024\Inshore_monsoon_Tentu_V1_N_20240926.xlsx')); %path to .mat structs of processed data
        northings_fixed = 2524542;
        eastings_fixed = 254699;
    elseif isequal(channel,'Tetulia') & isequal(tide,'Neap') & isequal(fixed_point, 'V2')
        file_name = (strcat(folder_adcp_fixed, 'Process Data\Neap\Tentulia River\26092024\Inshore_monsoon_Tentu_V2_N_20240926.xlsx')); %path to .mat structs of processed data
        northings_fixed = 2522717;
        eastings_fixed = 256554;
    elseif isequal(channel,'Meghna') & isequal(tide,'Neap') & isequal(fixed_point, 'V1')
        file_name = (strcat(folder_adcp_fixed, 'Process Data\Neap\Meghna River\28092024\Inshore_monsoon_Meg_V1_N_20240928.xlsx')); %path to .mat structs of processed data
        northings_fixed = 2530028;
        eastings_fixed = 255347;
    elseif isequal(channel,'Meghna') & isequal(tide,'Neap') & isequal(fixed_point, 'V3')
        file_name = (strcat(folder_adcp_fixed, 'Process Data\Neap\Meghna River\28092024\Inshore_monsoon_Meg_V3_N_20240928.xlsx')); %path to .mat structs of processed data
        northings_fixed = 2532024;
        eastings_fixed = 259171;

    elseif isequal(channel,'Tetulia') & isequal(tide,'Spring') & isequal(fixed_point, 'V1')
        file_name = (strcat(folder_adcp_fixed, 'Process Data\Spring\Tentulia\03102024\Inshore_monsoon_Tentu_V1_S_20241003.xlsx')); %path to .mat structs of processed data
        northings_fixed = 2524562;
        eastings_fixed = 254654;
    elseif isequal(channel,'Tetulia') & isequal(tide,'Spring') & isequal(fixed_point, 'V2')
        file_name = (strcat(folder_adcp_fixed, 'Process Data\Spring\Tentulia\03102024\Inshore_monsoon_Tentu_V2_S_20241003.xlsx')); %path to .mat structs of processed data
        northings_fixed = 2522738;
        eastings_fixed = 256443;
    elseif isequal(channel,'Meghna') & isequal(tide,'Spring') & isequal(fixed_point, 'V1')
        file_name = (strcat(folder_adcp_fixed, 'Process Data\Spring\Meghna\02102024\Inshore_monsoon_Meg_V1_S_20241002.xlsx')); %path to .mat structs of processed data
        northings_fixed = 2530099;
        eastings_fixed = 255331;
    elseif isequal(channel,'Meghna') & isequal(tide,'Spring') & isequal(fixed_point, 'V3')
        file_name = (strcat(folder_adcp_fixed, 'Process Data\Spring\Meghna\02102024\Inshore_monsoon_Meg_V3_S_20241002.xlsx')); %path to .mat structs of processed data
        northings_fixed = 2531150;
        eastings_fixed = 259826;
    else
        error('The given combination of Channel, Tide and point does not exist')
    end

    sheetNames = sheetnames(file_name);
    numSheets = numel(sheetNames);

    adcp_data = struct();
    sheet_names_store = {'preprocessed', 'WaterVelocity', 'WaterDirection'};
    for i = 1:numSheets
        data = readtable(file_name, 'Sheet', sheetNames(i));

        if i~= 1
            data = data(~isnat(data.Var2), :);
            years = year(data.Var2);
            data = data(years == 2024, :);
            num_columns = width(data);
            column_names = {'num', 'Date', arrayfun(@(x) ['depth' num2str(x)], 1:num_columns-2, 'UniformOutput', false)};
            column_names = [column_names{:}];
            data.Properties.VariableNames = column_names;           % Assign custom column names to the data
        else
            data = data(~isnat(data.Date_Time), :);
            years = year(data.Date_Time);
            data = data(years == 2024, :);
        end
        adcp_data.(sheet_names_store{i}) = data;
    end

    % Excel depth avereged velocity
    water_velocity_table = adcp_data.WaterVelocity(~isnat(adcp_data.WaterVelocity.Date), :);
    water_velocity = water_velocity_table{:, 3:end};
    mean_velocity = mean(water_velocity, 2, 'omitnan');
    water_velocity_table.AvgVelocity = mean_velocity;

    % Excel depth avereged direction
    water_direction_table = adcp_data.WaterDirection(~isnat(adcp_data.WaterDirection.Date), :);
    water_direction = water_direction_table{:, 3:end};

    avg_directions = NaN(height(water_direction_table), 1);

    for i = 1:height(water_direction)
        % Get the directions for the current row, ignoring NaN values
        current_directions = water_direction(i, ~isnan(water_direction(i, :)));

        radians = deg2rad(current_directions);

        % Calculate x and y components
        x = cos(radians);
        y = sin(radians);

        % Calculate the average x and y
        avg_x = mean(x);
        avg_y = mean(y);

        % Calculate the average direction using atan2
        avg_angle_rad = atan2(avg_y, avg_x);

        % Convert back to degrees
        avg_directions(i) = mod(rad2deg(avg_angle_rad)+360,360);
    end
    water_direction_table.AvgDirection = avg_directions;

    [commonDates, idx1, idx2] = intersect(water_direction_table.Date, water_velocity_table.Date);

    % Filter the tables to keep only rows with common datetime values
    water_direction_table = water_direction_table(idx1, :);
    water_velocity_table = water_velocity_table(idx2, :);



    % Convert to U,V

    water_velocity_fixed = water_velocity_table(:,["Date","AvgVelocity"]);
    water_direction_fixed = water_direction_table(:,["Date","AvgDirection"]);
    vel_x = nan(height(water_velocity_fixed),1);
    vel_y = nan(height(water_velocity_fixed),1);
    AvgDirection = nan(height(water_velocity_fixed),1);
    for i = 1:height(water_direction_fixed)
        AvgDirection(i) = water_direction_fixed.AvgDirection(i);
        radians = deg2rad(water_direction_fixed.AvgDirection(i));
        velocity = water_velocity_fixed.AvgVelocity(i);
        vel_x(i) = sin(radians)*velocity;
        vel_y(i) = cos(radians)*velocity;
        % % Display each value on a separate line
        %
        % disp(['Velocity: ', num2str(velocity), ' m/s']);
        % disp(['Direction: ', num2str(water_direction_short.AvgDirection(i)), ' deg']);
        % disp(['U component: ', num2str(vel_U), ' m/s']);
        % disp(['V component: ', num2str(vel_V), ' m/s']);
    end
    water_velocity_fixed.AvgDirection = AvgDirection;
    water_velocity_fixed.vel_x = vel_x;
    water_velocity_fixed.vel_y = vel_y;

    if save_fixed_point
        save(strcat('Matlab_data\Fixed_stations\Fixed_point_',channel,'_', tide, '_', fixed_point, '.mat'), 'water_velocity_fixed')
        save(strcat('Matlab_data\Fixed_stations\Fixed_point_loc_N_',channel,'_', tide, '_', fixed_point, '.mat'), 'northings_fixed')
        save(strcat('Matlab_data\Fixed_stations\Fixed_point_loc_E_',channel,'_', tide, '_', fixed_point, '.mat'), 'eastings_fixed')
    end
end

function [F, fixed_vel_s, fixed_vel_n, northings_fixed, eastings_fixed] = get_fixed_point_pd0(channel, tide, constituents, fixed_point, xs)
    folder_adcp_fixed = 'C:\Users\joris\OneDrive - Wageningen University & Research\Documenten\Wageningen\Master\Thesis\Data\Fixed Boat Velocity_Sentinel V\';
    if isequal(channel,'Tetulia') & isequal(tide,'Neap') & isequal(fixed_point, 'V1')
        adcp = rdi.readADCP(strcat(folder_adcp_fixed, 'Raw Data\Neap\Tetlulia\Sentinel V\Tetlulia\V1\Tentulia 15m Depth 20240926T005223.pd0')); %path to .mat structs of processed data
        northings_fixed = 2524542;
        eastings_fixed = 254699;
    elseif isequal(channel,'Tetulia') & isequal(tide,'Neap') & isequal(fixed_point, 'V2')
        adcp = rdi.readADCP(strcat(folder_adcp_fixed, 'Raw Data\Neap\Tetlulia\Sentinel V\Tetlulia\V2\Tentulia_18m Velocity Profiling 20240926T003451.pd0')); %path to .mat structs of processed data
        northings_fixed = 2522717;
        eastings_fixed = 256554;
    elseif isequal(channel,'Meghna') & isequal(tide,'Neap') & isequal(fixed_point, 'V1')
        adcp = rdi.readADCP(strcat(folder_adcp_fixed, 'Raw Data\Neap\Meghna\Sentinel V\V1\28092024\Meghna 15m Depth 20240928T050035.pd0')); %path to .mat structs of processed data
        northings_fixed = 2530028;
        eastings_fixed = 255347;
    elseif isequal(channel,'Meghna') & isequal(tide,'Neap') & isequal(fixed_point, 'V3')
        adcp = rdi.readADCP(strcat(folder_adcp_fixed, 'Raw Data\Neap\Meghna\Sentinel V\V3\28092024\Meghna_18m Velocity Profiling 20240927T214934(1).pd0')); %path to .mat structs of processed data
        northings_fixed = 2532024;
        eastings_fixed = 259171;

    elseif isequal(channel,'Tetulia') & isequal(tide,'Spring') & isequal(fixed_point, 'V1')
        adcp = rdi.readADCP(strcat(folder_adcp_fixed, 'Raw Data\Spring\Tetulia\Sentinel V\V1\Tetulia_18m Velocity Profiling_Spring 20241002T205821.pd0')); %path to .mat structs of processed data
        northings_fixed = 2524562;
        eastings_fixed = 254654;
    elseif isequal(channel,'Tetulia') & isequal(tide,'Spring') & isequal(fixed_point, 'V2')
        adcp = rdi.readADCP(strcat(folder_adcp_fixed, 'Raw Data\Spring\Tetulia\Sentinel V\V2\Tetulia 15m Depth_Spring 20241002T214359.pd0')); %path to .mat structs of processed data
        northings_fixed = 2522738;
        eastings_fixed = 256443;
    elseif isequal(channel,'Meghna') & isequal(tide,'Spring') & isequal(fixed_point, 'V1')
        adcp = rdi.readADCP(strcat(folder_adcp_fixed, 'Raw Data\Spring\Meghna\Sentinel V\V1\02102024\Meghna 15m Depth_Spring 20241001T174606.pd0')); %path to .mat structs of processed data
        northings_fixed = 2530099;
        eastings_fixed = 255331;
    elseif isequal(channel,'Meghna') & isequal(tide,'Spring') & isequal(fixed_point, 'V3')
        adcp = rdi.readADCP(strcat(folder_adcp_fixed, 'Raw Data\Spring\Meghna\Sentinel V\V3\02102024\Meghna_18m Velocity Profiling_Spring 20241001T174029.pd0')); %path to .mat structs of processed data
        northings_fixed = 2531150;
        eastings_fixed = 259826;
    else
        error('The given combination of Channel, Tide and point does not exist')
    end


    % Read water level data
    [water_level_time, water_level_value] = load_water(channel);


    %%% Restructure ADCP-data
    % Does not seem to have any effect on plot_velocity. Does it work?
    % Add waterlevel
    F = rdi.ADCP(adcp);
    F.water_level_object = VaryingWaterLevel(water_level_time, water_level_value);
    F.water_level_object.model = TidalScalarModel(constituents = constituents);
    F.water_level_object.model.scalar_name = 'eta';
    F.water_level_object.get_parameters();
    F.horizontal_position_provider.position = [eastings_fixed, northings_fixed];

    % figure;
    % F.plot_velocity

    velocity = F.velocity(CoordinateSystem.Earth);
    velocity_E = velocity(:,:,1);
    velocity_N = velocity(:,:,2);

    % rotate fixed velocity to mesh coords
    fixed_vel_s = zeros(size(velocity_E,1),1);
    fixed_vel_n = zeros(size(velocity_E,1),1);
    for i = 1:size(velocity_E,1)
        for j = 1:size(velocity_E,2)
            [fixed_vel_s(i,j), fixed_vel_n(i,j)] = xs.xy2sn_vel(velocity_E(i,j), velocity_N(i,j));
        end
    end
    %
    % figure;
    % F.plot_velocity();
end

function [fixed_vel_s_smoothed_clean, fixed_vel_s_smoothed_time, fixed_vel_n_smoothed_clean, fixed_vel_n_smoothed_time] = smooth_filter_fixed_point_pd0(F,fixed_vel_s, fixed_vel_n, visual_check)

    vel_pos={F.depth_cell_position};
    vel_pos=cellfun(@(x) mean(x(:,:,:,3),3,'omitnan'), vel_pos,...
        'UniformOutput', false);

    mean_depth = mean(mean(vel_pos{1}));
    for i = 1:size(fixed_vel_s,2)
        if median(vel_pos{1}(:,i))>=mean_depth/5

            fixed_vel_s(:,i) = nan;
            fixed_vel_n(:,i) = nan;
        end
    end




    timeDiff = [diff(F.time),seconds(500)];

    vel_s_single = [];
    vel_s_time = [];
    vel_s_smoothed = [];
    fixed_vel_s_smoothed_time = [];
    for i = 1:length(F.time)
        if timeDiff(i) < seconds(300)
            vel_s_single = [vel_s_single, fixed_vel_s(:,i)];
            vel_s_time = [vel_s_time, F.time(i)];
        else
            vel_s_single = [vel_s_single, fixed_vel_s(:,i)];
            vel_s_time = [vel_s_time, F.time(i)];

            vel_s_smoothed = [vel_s_smoothed, mean(vel_s_single,2, 'omitmissing')];
            fixed_vel_s_smoothed_time = [fixed_vel_s_smoothed_time, mean(vel_s_time, 'omitmissing')];

            vel_s_single = [];
            vel_s_time = [];
        end
    end

    vel_n_single = [];
    vel_n_time = [];
    vel_n_smoothed = [];
    fixed_vel_n_smoothed_time = [];
    for i = 1:length(F.time)
        if timeDiff(i) < seconds(300)
            vel_n_single = [vel_n_single, fixed_vel_n(:,i)];
            vel_n_time = [vel_n_time, F.time(i)];
        else
            vel_n_single = [vel_n_single, fixed_vel_n(:,i)];
            vel_n_time = [vel_n_time, F.time(i)];

            vel_n_smoothed = [vel_n_smoothed, mean(vel_n_single,2, 'omitmissing')];
            fixed_vel_n_smoothed_time = [fixed_vel_n_smoothed_time, mean(vel_n_time, 'omitmissing')];

            vel_n_single = [];
            vel_n_time = [];
        end
    end





    vel_s_smoothed_mean = mean(vel_s_smoothed, "omitmissing");
    vel_n_smoothed_mean = mean(vel_n_smoothed, "omitmissing");

    % Initialize a matrix to store cleaned data
    fixed_vel_s_smoothed_clean = NaN(size(vel_s_smoothed));  % Matrix to store filtered data
    fixed_vel_n_smoothed_clean = NaN(size(vel_n_smoothed));  % Matrix to store filtered data
    outliers = isoutlier(vel_s_smoothed_mean, 'movmedian', 30, 'ThresholdFactor',1.5);

    % Loop through each row (time series) and remove outliers
    for i = 1:size(vel_s_smoothed, 1)
        % Remove the outliers from the current time series
        fixed_vel_s_smoothed_clean(i, ~outliers) = vel_s_smoothed(i,~outliers);
        fixed_vel_n_smoothed_clean(i, ~outliers) = vel_n_smoothed(i,~outliers);
    end


    % if visual_check == 0
    %     figure
    %     plot(fixed_vel_s_smoothed_time,vel_s_smoothed_mean)
    %     hold on
    %     plot(fixed_vel_s_smoothed_time,mean(fixed_vel_s_smoothed_clean))
    %     plot(fixed_vel_n_smoothed_time,vel_n_smoothed_mean);
    %     plot(fixed_vel_n_smoothed_time,mean(fixed_vel_n_smoothed_clean));
    % end
end

function [water_velocity_fixed, northings_fixed, eastings_fixed] = load_fixed_point(channel, tide, fixed_point)
    load(strcat('Matlab_data\Fixed_stations\Fixed_point_',channel,'_', tide, '_', fixed_point, '.mat'))
    load(strcat('Matlab_data\Fixed_stations\Fixed_point_loc_N',channel,'_', tide, '_', fixed_point, '.mat'))
    load(strcat('Matlab_data\Fixed_stations\Fixed_point_loc_E',channel,'_', tide, '_', fixed_point, '.mat'))
end

function [point_loc_x, point_loc_y, water_velocity_fixed] = compare_transect_to_fixed(channel, tide, transect, constituents, fit_new, V, mesh, xs, pars_U, pars_V, pars_W, t_plot, flow_tracks, save_figure, model_name)


    % select fixed point
    if isequal(channel,'Tetulia') & isequal(transect,'T2')
        fixed_point = {'V2'};
    elseif isequal(channel,'Tetulia') & isequal(transect,'T3')
        fixed_point = {'V1'};
    elseif isequal(channel,'Meghna')
        fixed_point = {'V1', 'V3'} ;
    end

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

    point_loc_x = zeros(1,length(fixed_point));
    point_loc_y = zeros(1,length(fixed_point));

    % Create figure
    screenSize = get(0, 'ScreenSize');
    figure('Position', [1, 1, screenSize(3)-100, screenSize(4)-100]);
    hold on

    for i = 1:length(fixed_point)
        % fixed point from excel
        if fit_new
            [water_velocity_fixed, northings_fixed, eastings_fixed] = get_fixed_point_excel(channel, tide, fixed_point{i}, 1);
        else
            [water_velocity_fixed, northings_fixed, eastings_fixed] = load_fixed_point(channel, tide, fixed_point{i});
        end
        point_loc_x(i) = eastings_fixed;
        point_loc_y(i) = northings_fixed;

        % rotate fixed velocity to mesh coords
        vel_s = zeros(size(water_velocity_fixed,1),1);
        vel_n = zeros(size(water_velocity_fixed,1),1);
        for j = 1:size(water_velocity_fixed,1)
            [vel_s(j), vel_n(j)] = xs.xy2sn_vel(water_velocity_fixed.vel_x(j), water_velocity_fixed.vel_y(j));
        end
        water_velocity_fixed.vel_us = vel_s;
        water_velocity_fixed.vel_vn = vel_n;

        % Fixed point from raw .pd0 file
        [F, fixed_vel_s_pd0, fixed_vel_n_pd0, ~, ~] = get_fixed_point_pd0(channel, tide, constituents, fixed_point{i}, xs);

        [fixed_vel_s_smoothed_clean, fixed_vel_s_smoothed_time, fixed_vel_n_smoothed_clean, fixed_vel_n_smoothed_time] = smooth_filter_fixed_point_pd0(F,fixed_vel_s_pd0, fixed_vel_n_pd0, 0);


        % Find closest points on mesh
        A = [eastings_fixed, northings_fixed];
        mesh_coords = [mesh.x_left', mesh.y_left'];
        % Calculate the Euclidean distances from point A to each coordinate
        distances = sqrt(sum((mesh_coords - A).^2, 2));
        % Find the indices of the two smallest distances
        [~, sortedIndices] = sort(distances);
        % Get the indices of the two closest coordinates
        column_number = sortedIndices(1); % first indices is the column number of the closest cell
        distance = distances(column_number);

        %Calculate depth averaged velocity in mesh column
        [~, ~, idx] = unique(mesh.col_to_cell);
        pars_U_col_mean  = [];
        pars_V_col_mean  = [];
        pars_W_col_mean  = [];
        for a = 1:size(pars_U,2)
            pars_U_col = pars_U(:,a);
            pars_V_col = pars_V(:,a);
            pars_W_col = pars_W(:,a);

            U_col_mean = accumarray(idx, pars_U_col, [], @(x) mean(x, 'omitnan'));
            V_col_mean = accumarray(idx, pars_V_col, [], @(x) mean(x, 'omitnan'));
            W_col_mean = accumarray(idx, pars_W_col, [], @(x) mean(x, 'omitnan'));

            pars_U_col_mean  = [pars_U_col_mean, U_col_mean ];
            pars_V_col_mean  = [pars_V_col_mean, V_col_mean ];
            pars_W_col_mean  = [pars_W_col_mean, W_col_mean ];
        end

        vel_U_boat_near_fixed = pars_U_col_mean(column_number,1);
        vel_V_boat_near_fixed = pars_V_col_mean(column_number,1);
        vel_W_boat_near_fixed = pars_W_col_mean(column_number,1);

        for constituent = 1:length(constituents)
            cur_constiturent = constituents{constituent};
            pars_col = constituent*2;
            period = eval(['o', cur_constiturent]);

            vel_U_boat_near_fixed = vel_U_boat_near_fixed + pars_U_col_mean(column_number,pars_col)*cos(period.*t_plot) + pars_U_col_mean(column_number,pars_col+1)*sin(period.*t_plot);
            vel_V_boat_near_fixed = vel_V_boat_near_fixed + pars_V_col_mean(column_number,pars_col)*cos(period.*t_plot) + pars_V_col_mean(column_number,pars_col+1)*sin(period.*t_plot);
            vel_W_boat_near_fixed = vel_W_boat_near_fixed + pars_W_col_mean(column_number,pars_col)*cos(period.*t_plot) + pars_W_col_mean(column_number,pars_col+1)*sin(period.*t_plot);
        end

        % Single track column average
        cells = (1:mesh.ncells);
        included_cells = (mesh.col_to_cell == column_number);
        CellID = cells(included_cells);

        flows_single_in_cell = table();
        time_single = table();
        for k = 1:numel(fieldnames(flow_tracks))
            measurement_time = median(V.time(V.fileid == k));

            track_name = strcat("Track_", num2str(k));
            newRow = flow_tracks.(track_name).pars(CellID,:);
            newRow = nanmean(newRow,1); % mesh column mean
            if isempty(flows_single_in_cell)
                % If NewTable is empty, just initialize it with the first newRow
                flows_single_in_cell = newRow;
                time_single = measurement_time;
            else
                % Concatenate the newRow to the NewTable
                flows_single_in_cell = [flows_single_in_cell; newRow];
                time_single = [time_single; measurement_time];
            end
        end


        subplot(length(fixed_point), 2, i*2-1)
        plot(datetime(t_plot/86400, 'ConvertFrom', 'datenum'), vel_U_boat_near_fixed, "LineWidth",2,"Color","k","DisplayName",'Depth averaged flow');
        hold on

        measure_time = [];
        measure_vel = [];
        for c = 1:length(CellID)
            Cellid = CellID(c);
            calc_vel_U = pars_U(Cellid,1);
            for constituent = 1:length(constituents)
                cur_constiturent = constituents{constituent};
                pars_col = constituent*2;
                period = eval(['o', cur_constiturent]);
                calc_vel_U = calc_vel_U + pars_U(Cellid,pars_col)*cos(period.*t_plot) + pars_U(Cellid,pars_col+1)*sin(period.*t_plot);
            end
            % plot(datetime(t_plot/86400, 'ConvertFrom', 'datenum'), calc_vel_U, "DisplayName",  ['Cell ' num2str(Cellid)]);

            % Add raw measurements
            [time_in_column, vel_in_cell] = extract_measured_velocity_mesh_cell(Cellid,V,xs,mesh,channel, 'U');
            measure_time = [measure_time, time_in_column];

            measure_vel = [measure_vel, vel_in_cell];
            size(measure_vel)
            % measurements_plot = plot(time_in_column,vel_in_cell,'k.','MarkerSize',4);

            % set(measurements_plot, 'HandleVisibility', 'off');
        end


        % 1. Sort time and data
        [timeSorted, sortIdx] = sort(measure_time);
        dataSorted = measure_vel(:, sortIdx);
        % --- Find breaks where time gap > 30 minutes ---
        timeDiff = minutes(diff(timeSorted));
        splitIdx = find(timeDiff > 30);
        splitIdx = [0, splitIdx, length(timeSorted)];

        % --- Prepare data for boxplot ---
        plotData = [];       % All data values
        groupIDs = [];       % One group ID per value
        groupLabels = {};    % Label for each group ID
        groupCount = 0;
        xSeparators = [];

        for d = 1:length(splitIdx)-1
            startIdx = splitIdx(d)+1;
            endIdx = splitIdx(d+1);

            segment = dataSorted(:, startIdx:endIdx);
            values = segment(~isnan(segment));
            t = timeSorted(startIdx);

            % Compute boxplot stats
            q1 = quantile(values, 0.25);
            q2 = quantile(values, 0.5);
            q3 = quantile(values, 0.75);
            iqr = q3 - q1;
            lowerWhisker = min(values(values >= q1 - 1.5*iqr));
            upperWhisker = max(values(values <= q3 + 1.5*iqr));
            outliers = values(values < lowerWhisker | values > upperWhisker);

            % Convert datetime to numeric for plotting (datetime axis still works)
            x = t; %datenum(t);

            % Width of the box (in days)
            w = 0.25/24;

            % Plot manually
            hold on;
            % Box

            h = fill([x-w x+w x+w x-w], [q1 q1 q3 q3], [0.2 0.5 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Median line
            h = plot([x-w x+w], [q2 q2], 'k', 'LineWidth', 1.5);
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Whiskers
            h = plot([x x], [q3 upperWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            h = plot([x x], [q1 lowerWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Whisker caps
            h = plot([x-w/2 x+w/2], [upperWhisker upperWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            h = plot([x-w/2 x+w/2], [lowerWhisker lowerWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Outliers (only plot if non-empty)
            if ~isempty(outliers)
                h = plot(repmat(x, size(outliers)), outliers, 'ko');
                if ~isempty(h)
                    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
            end

        end


        % boxplot(plotData, groupLabels, 'LabelOrientation', 'inline');
        % scatter(time_single(:,1), flows_single_in_cell(:,1), "DisplayName", 'Track average','MarkerFaceColor','b')

        plot(fixed_vel_s_smoothed_time,mean(fixed_vel_s_smoothed_clean,"omitmissing"), "LineWidth",2,"Color", "r", "DisplayName", 'Moored ADCP (pd0)')
        % plot(water_velocity_fixed.Date,water_velocity_fixed.vel_us, "LineWidth",2,"Color",[1 0.7 0.7], "DisplayName", 'Moored ADCP (IWM)')

        title_str = strrep(model_name, '_', '\_');
        title(['Comparison Vel U: ', strcat(title_str, ' vs. fixed', fixed_point{i})])
        ylabel('Velocity [m/s]')


        subplot(length(fixed_point), 2, i*2)
        plot(datetime(t_plot/86400, 'ConvertFrom', 'datenum'), vel_V_boat_near_fixed, "LineWidth",2,"Color","k","DisplayName",'Depth averaged flow');
        hold on

        measure_time = [];
        measure_vel = [];
        for c = 1:length(CellID)
            Cellid = CellID(c);
            calc_vel_V = pars_V(Cellid,1);
            for constituent = 1:length(constituents)
                cur_constiturent = constituents{constituent};
                pars_col = constituent*2;
                period = eval(['o', cur_constiturent]);
                calc_vel_V = calc_vel_V + pars_V(Cellid,pars_col)*cos(period.*t_plot) + pars_V(Cellid,pars_col+1)*sin(period.*t_plot);
            end

            % plot(datetime(t_plot/86400, 'ConvertFrom', 'datenum'), calc_vel_V, "DisplayName",  ['Cell ' num2str(Cellid)]);

            % Add raw measurements
            [time_in_column, vel_in_cell] = extract_measured_velocity_mesh_cell(Cellid,V,xs,mesh,channel, 'V');
            measure_time = [measure_time, time_in_column];

            measure_vel = [measure_vel, vel_in_cell];
            % measurements_plot = plot(time_in_column,vel_in_cell,'k.','MarkerSize',4);

            % set(measurements_plot, 'HandleVisibility', 'off');
        end

        % scatter(time_single(:,1), flows_single_in_cell(:,2), "DisplayName", 'Track average','MarkerFaceColor','b')


        % 1. Sort time and data
        [timeSorted, sortIdx] = sort(measure_time);
        dataSorted = measure_vel(:, sortIdx);
        % --- Find breaks where time gap > 30 minutes ---
        timeDiff = minutes(diff(timeSorted));
        splitIdx = find(timeDiff > 30);
        splitIdx = [0, splitIdx, length(timeSorted)];

        % --- Prepare data for boxplot ---
        plotData = [];       % All data values
        groupIDs = [];       % One group ID per value
        groupLabels = {};    % Label for each group ID
        groupCount = 0;
        xSeparators = [];

        for d = 1:length(splitIdx)-1
            startIdx = splitIdx(d)+1;
            endIdx = splitIdx(d+1);

            segment = dataSorted(:, startIdx:endIdx);
            values = segment(~isnan(segment));
            t = timeSorted(startIdx);

            % Compute boxplot stats
            q1 = quantile(values, 0.25);
            q2 = quantile(values, 0.5);
            q3 = quantile(values, 0.75);
            iqr = q3 - q1;
            lowerWhisker = min(values(values >= q1 - 1.5*iqr));
            upperWhisker = max(values(values <= q3 + 1.5*iqr));
            outliers = values(values < lowerWhisker | values > upperWhisker);

            % Convert datetime to numeric for plotting (datetime axis still works)
            x = t; %datenum(t);

            % Width of the box (in days)
            w = 0.25/24;

            % Plot manually
            hold on;
            % Box
            h = fill([x-w x+w x+w x-w], [q1 q1 q3 q3], [0.2 0.5 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Median line
            h = plot([x-w x+w], [q2 q2], 'k', 'LineWidth', 1.5);
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Whiskers
            h = plot([x x], [q3 upperWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            h = plot([x x], [q1 lowerWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Whisker caps
            h = plot([x-w/2 x+w/2], [upperWhisker upperWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            h = plot([x-w/2 x+w/2], [lowerWhisker lowerWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Outliers (only plot if non-empty)
            if ~isempty(outliers)
                h = plot(repmat(x, size(outliers)), outliers, 'ko');
                if ~isempty(h)
                    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
            end

        end
        plot(fixed_vel_n_smoothed_time,mean(fixed_vel_n_smoothed_clean,"omitmissing"), "LineWidth",2,"Color","r", "DisplayName", 'Moored ADCP-(raw)')
        % plot(water_velocity_fixed.Date,water_velocity_fixed.vel_vn, "LineWidth",2,"Color",[1 0.7 0.7], "DisplayName", 'Moored ADCP')

        title_str = strrep(model_name, '_', '\_');
        title(['Comparison Vel V: ', strcat(title_str, ' vs. fixed', fixed_point{i})])
        % xlabel('velocity [m/s]')
        legend


    end

    if save_figure
        fig_location = strcat("Matlab_data\Figures_models\Fixed_point_comparison\",model_name, "_Fixed_point3.png");
        saveas(gcf,fig_location)
    end

end

function [point_loc_x, point_loc_y, water_velocity_fixed] = compare_transect_to_fixed_simp(channel, tide, transect, fixed_point, constituents, fit_new, V, mesh, xs, pars_U, pars_V, pars_W, t_plot, flow_tracks, save_figure, model_name)




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

    point_loc_x = zeros(1,length(fixed_point));
    point_loc_y = zeros(1,length(fixed_point));

    % Create figure
    fig = figure;
    fig.Units = 'centimeters';
    fig.Position = [0.5, 0.5, 10, 10];
    hold on

    for i = 1:length(fixed_point)
        % fixed point from excel
        if fit_new
            [water_velocity_fixed, northings_fixed, eastings_fixed] = get_fixed_point_excel(channel, tide, fixed_point{i}, 1);
        else
            [water_velocity_fixed, northings_fixed, eastings_fixed] = load_fixed_point(channel, tide, fixed_point{i});
        end
        point_loc_x(i) = eastings_fixed;
        point_loc_y(i) = northings_fixed;

        % rotate fixed velocity to mesh coords
        vel_s = zeros(size(water_velocity_fixed,1),1);
        vel_n = zeros(size(water_velocity_fixed,1),1);
        for j = 1:size(water_velocity_fixed,1)
            [vel_s(j), vel_n(j)] = xs.xy2sn_vel(water_velocity_fixed.vel_x(j), water_velocity_fixed.vel_y(j));
        end
        water_velocity_fixed.vel_us = vel_s;
        water_velocity_fixed.vel_vn = vel_n;

        % Fixed point from raw .pd0 file
        [F, fixed_vel_s_pd0, fixed_vel_n_pd0, ~, ~] = get_fixed_point_pd0(channel, tide, constituents, fixed_point{i}, xs);

        [fixed_vel_s_smoothed_clean, fixed_vel_s_smoothed_time, fixed_vel_n_smoothed_clean, fixed_vel_n_smoothed_time] = smooth_filter_fixed_point_pd0(F,fixed_vel_s_pd0, fixed_vel_n_pd0, 0);


        % Find closest points on mesh
        A = [eastings_fixed, northings_fixed];
        mesh_coords = [mesh.x_left', mesh.y_left'];
        % Calculate the Euclidean distances from point A to each coordinate
        distances = sqrt(sum((mesh_coords - A).^2, 2));
        % Find the indices of the two smallest distances
        [~, sortedIndices] = sort(distances);
        % Get the indices of the two closest coordinates
        column_number = sortedIndices(1); % first indices is the column number of the closest cell
        distance = distances(column_number);

        %Calculate depth averaged velocity in mesh column
        [~, ~, idx] = unique(mesh.col_to_cell);
        pars_U_col_mean  = [];
        pars_V_col_mean  = [];
        pars_W_col_mean  = [];
        for a = 1:size(pars_U,2)
            pars_U_col = pars_U(:,a);
            pars_V_col = pars_V(:,a);
            pars_W_col = pars_W(:,a);

            U_col_mean = accumarray(idx, pars_U_col, [], @(x) mean(x, 'omitnan'));
            V_col_mean = accumarray(idx, pars_V_col, [], @(x) mean(x, 'omitnan'));
            W_col_mean = accumarray(idx, pars_W_col, [], @(x) mean(x, 'omitnan'));

            pars_U_col_mean  = [pars_U_col_mean, U_col_mean ];
            pars_V_col_mean  = [pars_V_col_mean, V_col_mean ];
            pars_W_col_mean  = [pars_W_col_mean, W_col_mean ];
        end

        vel_U_boat_near_fixed = pars_U_col_mean(column_number,1);
        vel_V_boat_near_fixed = pars_V_col_mean(column_number,1);
        vel_W_boat_near_fixed = pars_W_col_mean(column_number,1);

        for constituent = 1:length(constituents)
            cur_constiturent = constituents{constituent};
            pars_col = constituent*2;
            period = eval(['o', cur_constiturent]);

            vel_U_boat_near_fixed = vel_U_boat_near_fixed + pars_U_col_mean(column_number,pars_col)*cos(period.*t_plot) + pars_U_col_mean(column_number,pars_col+1)*sin(period.*t_plot);
            vel_V_boat_near_fixed = vel_V_boat_near_fixed + pars_V_col_mean(column_number,pars_col)*cos(period.*t_plot) + pars_V_col_mean(column_number,pars_col+1)*sin(period.*t_plot);
            vel_W_boat_near_fixed = vel_W_boat_near_fixed + pars_W_col_mean(column_number,pars_col)*cos(period.*t_plot) + pars_W_col_mean(column_number,pars_col+1)*sin(period.*t_plot);
        end

        % Single track column average
        cells = (1:mesh.ncells);
        included_cells = (mesh.col_to_cell == column_number);
        CellID = cells(included_cells);

        flows_single_in_cell = table();
        time_single = table();
        for k = 1:numel(fieldnames(flow_tracks))
            measurement_time = median(V.time(V.fileid == k));

            track_name = strcat("Track_", num2str(k));
            newRow = flow_tracks.(track_name).pars(CellID,:);
            newRow = nanmean(newRow,1); % mesh column mean
            if isempty(flows_single_in_cell)
                % If NewTable is empty, just initialize it with the first newRow
                flows_single_in_cell = newRow;
                time_single = measurement_time;
            else
                % Concatenate the newRow to the NewTable
                flows_single_in_cell = [flows_single_in_cell; newRow];
                time_single = [time_single; measurement_time];
            end
        end


        % subplot(length(fixed_point), 2, i*2-1)
        plot(datetime(t_plot/86400, 'ConvertFrom', 'datenum'), vel_U_boat_near_fixed, "LineWidth",2,"Color","k","DisplayName",'Modelled flow');
        hold on

        measure_time = [];
        measure_vel = [];
        for c = 1:length(CellID)
            Cellid = CellID(c);
            calc_vel_U = pars_U(Cellid,1);
            for constituent = 1:length(constituents)
                cur_constiturent = constituents{constituent};
                pars_col = constituent*2;
                period = eval(['o', cur_constiturent]);
                calc_vel_U = calc_vel_U + pars_U(Cellid,pars_col)*cos(period.*t_plot) + pars_U(Cellid,pars_col+1)*sin(period.*t_plot);
            end
            % plot(datetime(t_plot/86400, 'ConvertFrom', 'datenum'), calc_vel_U, "DisplayName",  ['Cell ' num2str(Cellid)]);

            % Add raw measurements
            [time_in_column, vel_in_cell] = extract_measured_velocity_mesh_cell(Cellid,V,xs,mesh,channel, 'U');
            measure_time = [measure_time, time_in_column];

            measure_vel = [measure_vel, vel_in_cell];
            size(measure_vel)
            % measurements_plot = plot(time_in_column,vel_in_cell,'k.','MarkerSize',4);

            % set(measurements_plot, 'HandleVisibility', 'off');
        end


        % 1. Sort time and data
        [timeSorted, sortIdx] = sort(measure_time);
        dataSorted = measure_vel(:, sortIdx);
        % --- Find breaks where time gap > 30 minutes ---
        timeDiff = minutes(diff(timeSorted));
        splitIdx = find(timeDiff > 30);
        splitIdx = [0, splitIdx, length(timeSorted)];

        % --- Prepare data for boxplot ---
        plotData = [];       % All data values
        groupIDs = [];       % One group ID per value
        groupLabels = {};    % Label for each group ID
        groupCount = 0;
        xSeparators = [];

        for d = 1:length(splitIdx)-1
            startIdx = splitIdx(d)+1;
            endIdx = splitIdx(d+1);

            segment = dataSorted(:, startIdx:endIdx);
            values = segment(~isnan(segment));
            t = timeSorted(startIdx);

            % Compute boxplot stats
            q1 = quantile(values, 0.25);
            q2 = quantile(values, 0.5);
            q3 = quantile(values, 0.75);
            iqr = q3 - q1;
            lowerWhisker = min(values(values >= q1 - 1.5*iqr));
            upperWhisker = max(values(values <= q3 + 1.5*iqr));
            outliers = values(values < lowerWhisker | values > upperWhisker);

            % Convert datetime to numeric for plotting (datetime axis still works)
            x = t; %datenum(t);

            % Width of the box (in days)
            w = 0.25/24;

            % Plot manually
            hold on;
            % Box

            h = fill([x-w x+w x+w x-w], [q1 q1 q3 q3], [0.2 0.5 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Median line
            h = plot([x-w x+w], [q2 q2], 'k', 'LineWidth', 1.5);
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Whiskers
            h = plot([x x], [q3 upperWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            h = plot([x x], [q1 lowerWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Whisker caps
            h = plot([x-w/2 x+w/2], [upperWhisker upperWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            h = plot([x-w/2 x+w/2], [lowerWhisker lowerWhisker], 'k');
            if ~isempty(h)
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            % Outliers (only plot if non-empty)
            if ~isempty(outliers)
                h = plot(repmat(x, size(outliers)), outliers, 'ko');
                if ~isempty(h)
                    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
            end

        end


        % boxplot(plotData, groupLabels, 'LabelOrientation', 'inline');
        % scatter(time_single(:,1), flows_single_in_cell(:,1), "DisplayName", 'Track average','MarkerFaceColor','b')

        plot(fixed_vel_s_smoothed_time,mean(fixed_vel_s_smoothed_clean,"omitmissing"), "LineWidth",2,"Color", "r", "DisplayName", 'Moored ADCP')
        % plot(water_velocity_fixed.Date,water_velocity_fixed.vel_us, "LineWidth",2,"Color",[1 0.7 0.7], "DisplayName", 'Moored ADCP (IWM)')

        title_str = strrep(model_name, '_', '\_');
        % title(['Comparison Vel U: ', strcat(title_str, ' vs. fixed', fixed_point{i})])
        ylabel('Velocity [m/s]')


        
        legend('Location','southwest')


    end

    if save_figure
        fig_location = strcat("Matlab_data\Figures_models\Fixed_point_comparison\",model_name, "_Fixed_point4.png");
        saveas(gcf,fig_location)
    end

end

function plot_areal_image(V, xs, channel, point_loc_x, point_loc_y, save_figure, model_name)
    figure;
    [lat,lon] = utm2ll(V.horizontal_position(1,:),V.horizontal_position(2,:),46);
    geoplot(lat,lon, 'b.')
    hold on
    [lat,lon] = utm2ll(xs.origin(1),xs.origin(2),46);
    geoplot(lat,lon, 'r.')
    [lat_end,lon_end] = utm2ll(xs.origin(1)+xs.direction(1)*xs.scale,xs.origin(2)+xs.direction(2)*xs.scale,46);
    geoplot([lat lat_end],[lon lon_end],'LineWidth',1,'Color','g');
    [lat_end,lon_end] = utm2ll(xs.origin(1)+xs.direction_orthogonal(1)*xs.scale,xs.origin(2)+xs.direction_orthogonal(2)*xs.scale,46);
    geoplot([lat lat_end],[lon lon_end],'LineWidth',1,'Color','r');
    [lat,lon] = utm2ll(point_loc_x,point_loc_y,46);
    geoplot(lat,lon, 'g.', MarkerSize=10)

    if isequal(channel,'Meghna')
        geolimits([22.8331 22.9027],[90.5865 90.6892])
    else
        geolimits([22.7816 22.8260],[90.5931 90.6407])
    end
    geobasemap satellite
    hold off
    if save_figure
        fig_location = strcat("Matlab_data\Figures_models\areal_figure\",model_name, "_areal.png");
        saveas(gcf,fig_location)
    end
end



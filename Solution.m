classdef Solution < handle & helpers.ArraySupport
    %Solution Class for capturing model parameters as output of
    %get_parameters

    %   Detailed explanation goes here

    properties(SetAccess = ?Solver)
        M (:,:) double % matrix M such that b = Mp

        b (:,1) double % rhs of system of eqs (b = Mp)

        p (:,:) double % model parameters (b = Mp)

        cell_idx (:,1) double

        solver (1,1)

        name

        time (:,1) datetime % measurement time data
        fileid (:,1) double % fileid of the measuerment data

    end

    %  Solver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   ensemble_filter - defines repeat transects to be processed
    %   data_model - defines the velocity model to fit

    properties(Dependent)
        pars %p arranged in matrix of size Nc x np.
        ns   %number of measurements per cell.
    end
    methods
        function obj = Solution(varargin)
            obj = obj@helpers.ArraySupport(varargin{:})
        end
        function pars = get.pars(obj)
            pars = obj.p2pars(obj.p);
        end
        function ns = get.ns(obj)
            ns = accumarray(obj.cell_idx, ones(size(obj.cell_idx)),[obj.solver.mesh.ncells, 1]);
        end
        function varargout = get_data(obj)
            varargout = cell(1,nargout);
            [varargout{:}] = obj.solver.model.get_data(obj.pars);
        end

        function gof = validate(obj)
            if ~isscalar(obj)
                obj.run_method('validate');
                return
            end
            gof = obj.get_residuals();
        end

        function res = evaluate(obj, sol_idx, X)
            % Wrapper for model.get_data, with arbitrary number of query
            % points
            % input: sol_idx: index of the solution to be evaluated
            % X.T : Nq x 1 vector of time values
            % X.N : Nq x 1 vector of lateral coordinate values
            % X.Sig : Nq x 1 vector of sigma values

            % evaluate(obj, sol_idx, T = ..., N = ..., Sig = ...)
            % Used for interpolation or plotting
            arguments
                obj
                sol_idx (:,:) double  = 1
                X.T (:,1) double = 0
                X.dX (1,1) double = 0 % scalar used to investigate along-channel gradients.
                X.N (:,1) double = 0
                X.Sig (:,1) double = 0
                X.extrapolate (1,1) logical = false
            end
            % if statements: numel XT, XN, XSig must be equal
            % npars must have equal elements
            nx = numel(X.T);
            np = obj.solver.model.npars(1);

            % find mesh cell query points belong to
            ev_cidx = obj.solver.mesh.index(X.N, X.Sig);

            % Query points too close to bottom or surface.

            if X.extrapolate
                ev_cidx = obj.solver.mesh.index_extrapolate(ev_cidx, X.N, X.Sig);
            else
                error("Query points outside mesh domain. " + ...
                    "Enter extrapolate = true to linearly " + ...
                    "extrapolate the solution.")
            end
            n_center = reshape(...
                obj.solver.mesh.n_middle(obj.solver.mesh.col_to_cell), [], 1);
            dS = X.dX*ones(size(ev_cidx));
            dN = X.N - n_center(ev_cidx); % delta_n
            dZ = zeros(size(ev_cidx)); % delta_z - not used, dSig takes precedence
            dSig = X.Sig - obj.solver.mesh.sig_center(ev_cidx); % delta_sig

            dT = datetime(X.T, 'ConvertFrom', 'datenum');
            M0 = obj.solver.model.get_model(dT, dS, dN, dZ,dSig);

            pars_cell = obj.pars(ev_cidx, :, sol_idx);


            ip = repmat(1:nx, [np, 1]); ip = ip(:);
            jp = 1:nx*np;
            res = nan([nx, obj.solver.model.ncomponents, numel(sol_idx)]);
            v = cell([obj.solver.model.ncomponents,1]);
            P = cell([obj.solver.model.ncomponents,numel(sol_idx)]);
            Mp = cell([obj.solver.model.ncomponents,1]);
            pidx = [0, cumsum(obj.solver.model.npars)];

            %The following could be coded more elegantly but does not
            %affect computation times significantly (probably
            %matmult/pagemtimes/permute)
            for dim = 1:obj.solver.model.ncomponents
                v{dim} = M0(:,:,dim)'; v{dim} = v{dim}(:);
                Mp{dim} = sparse(ip,jp,v{dim});
                for ri = 1:numel(sol_idx)
                    P{dim, ri} = squeeze(pars_cell(:,(pidx(dim)+1):pidx(dim+1), ri))'; P{dim, ri} = P{dim, ri}(:);
                    res(:,dim, ri) = Mp{dim}*P{dim, ri}; % Compute u = Mx*px;
                end
            end
        end

        function [H, Wl, Zb] = get_depth(obj)
            [xvec, yvec] = obj.solver.mesh.xs.sn2xy(zeros(size(obj.decomp.X.y)), obj.decomp.X.y); % convert to x,y
            zvec = obj.bathy.get_bed_elev(xvec, yvec);  % get (time-indep) bathy

            wl = obj.solver.adcp.water_level_object.get_water_level_model(datetime(obj.decomp.X.t, 'ConvertFrom', 'datenum'));


            [Wl, Zb] = ndgrid(wl, zvec);

            H = Wl - Zb;
        end


        function plot_solution(obj, plot_names, opts)
            % Wrapper for mesh.plot()
            % Does not support quiver plots
            arguments
                obj
                plot_names (1,:) cell = obj.solver.model.all_names
                opts.sol_idx (1,:) double  = 1:size(obj.p,2)
                opts.representation = "ab"
            end

            if ~isscalar(obj)
                obj.run_method('plot_solution');
                return
            end

            nn = length(plot_names);
            np = sum(obj.solver.model.npars); % Number of parameters in each cell
            m = figure;%makefigure(9, 3*nn);

            if numel(opts.sol_idx)>1 % Compare different vectors
                t = tiledlayout(m, nn, numel(opts.sol_idx), TileSpacing = "tight", Padding = "tight", TileIndexing = "columnmajor");
            else
                t = tiledlayout(m, 'flow', TileSpacing="tight", Padding="tight");
            end

            t.XLabel.String = 'y';
            t.YLabel.String = 'z';
            t.XLabel.Interpreter = 'latex';
            t.YLabel.Interpreter = 'latex';

            par_idx = obj.get_par_idx(plot_names);
            if strcmp(opts.representation, "Aphi")
                %pa = obj.p2pars();
                [tid_pars, tid_names] = ab2Aphi(obj.pars, obj.solver.model.all_names);
                P = obj.pars2p(tid_pars);
                P = P(:, opts.sol_idx);
                plot_names = tid_names(par_idx);
            else
                P = obj.p(:, opts.sol_idx);
            end

            titles = obj.modify_names(plot_names);

            for col = 1:numel(opts.sol_idx)
                for row = 1:nn
                    ax = nexttile(t);
                    amax = max(abs(P(par_idx(row):np:end, numel(opts.sol_idx))), [], 'omitnan') + 1e-5; % For paper, we assume nreg = 3
                    var = P(par_idx(row):np:end, col);
                    hold on
                    [hbed, hwater, hmesh] = obj.solver.mesh.plot(ax,'var', var, 'FixAspectRatio', false);
                    hmesh.LineStyle = 'none';

                    % Colormap
                    % Default: velmap
                    colormap(ax, helpers.cmaps("velmap"))
                    clim([-amax, amax])
                    if contains(titles{row}, 'A')    % linear variable
                        clim([-amax, amax])
                        colormap(ax, helpers.cmaps("Amap"))
                    elseif contains(titles{row}, 'phi')
                        clim([-pi, pi])              % cyclic variable
                        colormap(ax, helpers.cmaps("phimap"))
                    end

                    % Titles
                    lam = {'\mathbf{\lambda}_N', '\mathbf{\lambda}_L', '\mathbf{\lambda}_H'};
                    if 0
                        title(['$', titles{row}, ', \hspace{.1cm}  \lambda = ', lam{col}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                    else
                        title(['$', titles{row}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                    end

                    % Colorbar
                    if col == numel(opts.sol_idx)
                        c = colorbar;
                        c.TickLabelInterpreter = 'latex';
                        c.FontSize = 10;
                        if false
                            if ~contains(titles{row}, '\partial') % Velocities
                                unit = '[m/s]';
                            elseif contains(titles{row}, '\sigma') % Velocities
                                unit = '[m/s]';
                            else % derivatives of velocities in x,y,z directions: m/s/m = 1/s
                                unit = '[1/s]';
                            end



                            pos = get(c,'Position');
                            if row == 1
                                pos1 = pos;
                            end
                            % disp(pos)
                            c.Label.String = unit;
                            c.Label.Interpreter = 'latex';
                            %c.Label.Position(1) = pos1(1) + 3; % to change its position
                            %c.Label.Position(2) = c.Label.Position(2) + .3; % to change its position
                            c.Label.VerticalAlignment = 'middle'; % to change its position
                            c.Label.HorizontalAlignment = 'right';

                            c.Label.Rotation = 270; % to rotate the text

                        end
                    end
                    axis tight
                    set(gca, 'XDir','reverse') % Very important

                    set(gca, 'XTick',[])
                    set(gca, 'YTick',[])

                end
            end
            % Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
            % %[Left Bottom Right Top] spacing
            % NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            % set(gca, 'Position', NewPos);
        end


        function pars = p2pars(obj,p)
            [np, nsols] = size(p);
            ncells = obj.solver.mesh.ncells;
            npars = np/ncells;
            pars = reshape(obj.p, npars, ncells, nsols);
            pars = permute(pars, [2 1 3]);
        end

        function p = pars2p(obj, pars)
            p = nan(size(obj.p));
            for i=1:size(pars,3)
                p(:,i) = reshape(pars(:,:,i)', 1, [])';
            end
        end

    function [CV, rpc, rps] = cross_validate_2D_track(obj, min, max, N, k, track)
            % 2D analysis: two-element inputs min, max, N. (i.e. min =
            % [0,0])
            reg_pars_cont = obj.reg_pars_symlog(min(1), max(1), N(1));
            reg_pars_smooth = obj.reg_pars_symlog(min(2), max(2), N(2));

            [rpc, rps] = meshgrid(reg_pars_cont, reg_pars_smooth);

            reg_pars = helpers.vecs2mat_full(reg_pars_cont, reg_pars_smooth);

            % paste in because of ordering of constraints.
            reg_pars_mat = [reg_pars(:,1), reg_pars(:,1), reg_pars(:,2), reg_pars(:,2), reg_pars(:,1)];

            CV = obj.cross_validate(reg_pars_mat,track);
    end

        function CV = cross_validate_0D(obj, k)
            % simple cross-validation of the solution obj.p.
            CV = obj.kfold_cross_validate([obj.solver.regularization.weight], k);
        end

        function CV = cross_validate_1D(obj, min, max, N, k)
            % 1D analysis: scalar min, max, N.
            %reg_pars_mat = reg_pars_symlog(obj, min, max)
            rp  = reg_pars_symlog(obj, min, max, N);
            reg_pars_mat = repmat(rp, 1, 5);
            CV = obj.kfold_cross_validate(reg_pars_mat, k);

            figure;
            plot(rp, [CV{:,1}])
            xlabel('reg pars')
            ylabel('generalization error')
            title('lambda vs scaled generalization error')
        end

 function [rho, reg_pars_mat, eta] = cross_validate_1D_track(obj, min, max, N, track)
            % 1D analysis: scalar min, max, N.
            %reg_pars_mat = reg_pars_symlog(obj, min, max)
            rp  = reg_pars_symlog(obj, min, max, N);
            reg_pars_mat = repmat(rp, 1, 5);
            [CV, eta] = obj.cross_validate(reg_pars_mat,track);
            rho = CV;
            figure
            plot(rp, [CV{:,1}])
            xlabel('reg pars')
            ylabel('generalization error')
            title('lambda vs scaled generalization error')
            
        end

        function [CV, rpc, rps] = depth_average_validation_track(obj, min, max, N, track)
            reg_pars_cont = obj.reg_pars_symlog(min(1), max(1), N(1));
            reg_pars_smooth = obj.reg_pars_symlog(min(2), max(2), N(2));

            [rpc, rps] = meshgrid(reg_pars_cont, reg_pars_smooth);

            reg_pars = helpers.vecs2mat_full(reg_pars_cont, reg_pars_smooth);

            % paste in because of ordering of constraints.
            reg_pars_mat = [reg_pars(:,1), reg_pars(:,1), reg_pars(:,2), reg_pars(:,2), reg_pars(:,1)];


            CV = obj.cross_validate_depth_average(reg_pars_mat,track);
        end

        function [CV, rpc, rps] = cross_validate_2D(obj, min, max, N, k)
            % 2D analysis: two-element inputs min, max, N. (i.e. min =
            % [0,0])
            reg_pars_cont = obj.reg_pars_symlog(min(1), max(1), N(1));
            reg_pars_smooth = obj.reg_pars_symlog(min(2), max(2), N(2));

            [rpc, rps] = meshgrid(reg_pars_cont, reg_pars_smooth);

            reg_pars = helpers.vecs2mat_full(reg_pars_cont, reg_pars_smooth);

            % paste in because of ordering of constraints.
            reg_pars_mat = [reg_pars(:,1), reg_pars(:,1), reg_pars(:,2), reg_pars(:,2), reg_pars(:,1)];

            
            CV = obj.kfold_cross_validate(reg_pars_mat, k);

        end
    end

    

    methods(Access=protected)

        function CV = kfold_cross_validate(obj, reg_pars_mat, k)
            % reg_pars_mat is a matrix of size nreg x 5, with the 5 known
            % regularization constraints.

            % This function is an ad hoc function and has not been
            % optimized using matrix algebra.
            CV = cell([size(reg_pars_mat, 1), 2]);
            i = 0;
            tot = k*size(reg_pars_mat, 1);
            fprintf('K-fold cross-validation percentage: %2.2f percent \n', 0)
            kfold_idx = obj.kfold_split_dataset(k); % k fold partition - random indices
            
            for rp = 1:size(reg_pars_mat, 1)
                E = zeros([size(reg_pars_mat, 1), k]);
                for test_fold = 1:k
                    test_idx = kfold_idx{test_fold};
                    p_train = zeros([size(obj.p,1), k]);
                    for train_fold = setdiff(1:k, test_fold) % loop over all subsets except test set
                        train_idx = kfold_idx{train_fold};
                        % Construct training matrix and data
                        M0 = obj.M(train_idx, :);
                        b0 = obj.b(train_idx);

                        Mp = M0'*M0;

                        p_train(:,train_fold) = obj.assemble_solve_single(M0, b0, Mp, reg_pars_mat(rp,:));
                    end
                    
                    p_mean = sum(p_train,2)/(k-1);
                    % having obtained an estimate (p_mean), apply it to the
                    % test set and compare to the real measurement data.
                    M1 = obj.M(test_idx,:);
                    b1 = obj.b(test_idx);
                    E(rp, test_fold) = mean((M1*p_mean - b1).^2); % Generalization error
                    i = i+1;
                    fprintf('K-fold cross-validation percentage: %2.2f percent \n', 100*i/tot)
                end
                CV{rp,1} = mean(E(rp,:));
            end
        end


        function kfold_idx = kfold_split_dataset(obj, k)
            nb = numel(obj.cell_idx);
            shuffled_idx = randperm(nb);
            limits = round(linspace(1, nb+1, k+1));
            for i = 1:k
                kfold_idx{i} = shuffled_idx(limits(i):(limits(i+1)-1));
            end
        end


        function [CV, S] = cross_validate(obj, reg_pars_mat,track)
            % reg_pars_mat is a matrix of size nreg x 5, with the 5 known
            % regularization constraints.

            % This function is an ad hoc function and has not been
            % optimized using matrix algebra.
            p_train = cell([size(reg_pars_mat, 1), 1]);
            p_size = cell([size(reg_pars_mat, 1), 1]);

            if strcmp(obj.solver.opts.cv_mode, 'random')
                nepochs = obj.solver.opts.cv_iter;
            else
                nepochs = 1;
            end
            niter = size(reg_pars_mat, 1)*nepochs;
            E = cell([size(reg_pars_mat, 1), 2]);
            S = cell([size(reg_pars_mat, 1), 1]);
            CV = cell([size(reg_pars_mat, 1), 2]);
            i = 0;
            fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
            for ep = 1:obj.solver.opts.cv_iter % randomized iterations: Only meaningful if cv_mode == 'random'
                train_idx = logical(obj.split_dataset(track));
                test_idx = ~train_idx;

                % Construct training matrix and data
                Mt = obj.M;
                bt = obj.b;

                M0 = obj.M(train_idx, :);
                b0 = obj.b(train_idx);

                M1 = obj.M(test_idx,:);
                b1 = obj.b(test_idx);

                Mp = M0'*M0;
                Mtp = Mt' *Mt;

                for rp = 1:size(reg_pars_mat, 1)
                    p_train{rp}(:, ep) = obj.assemble_solve_single(M0, b0, Mp, reg_pars_mat(rp,:));
                    i = i+1;
                    fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
                    E{rp,1}(1, ep) = mean((M1*p_train{rp}(:, ep) - b1).^2); % Generalization error
                    E{rp,2}(1, ep) = mean((M0*p_train{rp}(:, ep) - b0).^2); % Training error
                    p_size{rp}(:, ep) = obj.assemble_solve_single(Mt,bt,Mtp, reg_pars_mat(rp,:));
                    if rp == 1
                    S0 = Mt*p_size{rp}(:,ep);
                    end
                    S{rp,1}(1,ep) = mean(abs(Mt*p_size{rp}(:,ep)-S0));
                end
            end
            for rp = 1:size(reg_pars_mat, 1)
                %                 %p_avg{rp} = mean(p_train{rp}, 2); % k-fold cross-validated estimate
                CV{rp, 1} = mean(E{rp,1}); % Ensemble average
                CV{rp, 2} = mean(E{rp,2}); % Ensemble average
            end
        end

        function CV = cross_validate_depth_average(obj, reg_pars_mat,track, V,xs,mesh)
            % reg_pars_mat is a matrix of size nreg x 5, with the 5 known
            % regularization constraints.

            % This function is an ad hoc function and has not been
            % optimized using matrix algebra.
            % depth_av data
            cellID = [];
            vel_cell_size = size(1,length(mesh.z_top_right));
            for i =1:length(mesh.cell_to_mat)
            cellID = i;
            [time_in_column, vel_in_cell_tot] = depth_av_measuements(V, xs, mesh,cellID);
            vel_cell_size(:,cellID) = sqrt(vel_in_cell_tot(:,1)^2 + vel_in_cell_tot(:,2)^2 + vel_in_cell_tot(:,3)^2);
            end
            % depth_av model
            
            timemin = min(V.time(V.fileid == track));
            [minimum,minidx] = min(time_in_column-timemin);

            timemax = max(V.time(V.fileid == track));
            [maximum,maxidx] = min(time_in_column-timemax);

            Av_vel_in_cell = mean(vel_in_cell_tot(minidx:maxidx,:));

            p_train = cell([size(reg_pars_mat, 1), 1]);

            nepochs = 1;

            niter = size(reg_pars_mat, 1)*nepochs;
            E = cell([size(reg_pars_mat, 1), 2]);
            CV = cell([size(reg_pars_mat, 1), 2]);
            i = 0;
            fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
            for ep = 1:obj.solver.opts.cv_iter % randomized iterations: Only meaningful if cv_mode == 'random'
                train_idx = logical(obj.split_dataset(track));
                test_idx = ~train_idx;

                % Construct training matrix and data
                M = obj.M(:, :);
                b = obj.b(:);

                Mp = obj.M'*obj.M;

                for rp = 1:size(reg_pars_mat, 1)
                    p_train{rp}(:, ep) = obj.assemble_solve_single(M, b, Mp, reg_pars_mat(rp,:));
                    i = i+1;
                    fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
                    E{rp,1}(1, ep) = mean((M1*p_train{rp}(:, ep) - b1).^2); % Generalization error
                    E{rp,2}(1, ep) = mean((M0*p_train{rp}(:, ep) - b0).^2); % Training error
                end
            end
            for rp = 1:size(reg_pars_mat, 1)
                %                 %p_avg{rp} = mean(p_train{rp}, 2); % k-fold cross-validated estimate
                CV{rp, 1} = mean(E{rp,1}); % Ensemble average
                CV{rp, 2} = mean(E{rp,2}); % Ensemble average
            end
        end



        function training_idx = split_dataset(obj, track)
            %ci = vertcat(obj.cell_idx{:});
            ci = obj.cell_idx;
            training_idx = ones(size(ci));

            if strcmp(obj.solver.opts.cv_mode, 'none')
                training_idx = ones(size(ci));
            elseif strcmp(obj.solver.opts.cv_mode, 'random')
                tp = obj.solver.opts.training_perc;
                rand0 = rand(size(training_idx));
                training_idx = (rand0 <= tp);
            elseif strcmp(obj.solver.opts.cv_mode, 'omit_cells')
                oc = obj.solver.opts.omit_cells;
                for occ = 1:length(oc)
                    training_idx(ci == oc(occ)) = 0;
                end
            elseif strcmp(obj.solver.opts.cv_mode, 'kfold') % Disjoint partitions (in contrast to 'random'). see

            elseif strcmp(obj.solver.opts.cv_mode, 'omit_time') % to be implemented
            elseif strcmp(obj.solver.opts.cv_mode, 'track') % to be implemented
                if track == 14
                training_idx = obj.fileid ~= track & track-1;
                else
                training_idx = obj.fileid ~= track & track+1;
                end
            end
        end

        function rp = reg_pars_symlog(obj, min, max, N)%, res_near_zero = 1.0)
            %Variation upon logspace.
            ymin = helpers.symlog(min);
            ymax = helpers.symlog(max);
            rp = helpers.symexp(linspace(ymin, ymax, N))';
        end


        function reg_pars_sens_vec = vectorize_reg_pars(obj, cont, smooth)
            if strcmp(obj.solver.opts.reg_vary, 'coupled')
                L = obj.reg_pars_sens(1, [1,3]); % Only vary two reg. parameters
                n = length(L);
                [L{:}] = ndgrid(L{end:-1:1});
                L = cat(n+1,L{:});
                L = fliplr(reshape(L,[],n));
                RP = [L(:, 1), L(:,1), L(:,2), L(:,2), L(:,1)]; % Coupling of parameters
            elseif strcmp(obj.reg_vary, 'full')
                error('Full analysis is not supported. Change opts.reg_vary to coupled.')
                L = obj.reg_pars_sens; % Vary all reg. parameters
                n = length(L);
                [L{:}] = ndgrid(L{end:-1:1});
                L = cat(n+1,L{:});
                RP = fliplr(reshape(L,[],n));
            else
                error('Invalid reg_vary option')
            end

            for i = 1:size(RP,2)
                RP(:,i) = obj.reg_relative_weights(i);
                if obj.force_zero(i)
                    RP(:,i) = 0;
                end
            end
            reg_pars_sens_vec = RP;
        end



        function [A, rhs, L] = assemble_single(obj, M, b, Mp, regp)
            A = Mp + regp(1)*obj.solver.regularization(1).Cg+ regp(2)*obj.solver.regularization(2).Cg +...
                regp(3)*obj.solver.regularization(3).Cg + regp(4)*obj.solver.regularization(4).Cg + regp(5)*obj.solver.regularization(5).Cg;
            obj.solver.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
            L = ichol(A, obj.solver.opts.preconditioner_opts);
            rhs = M'*b + regp(5)*obj.solver.regularization(5).C'*obj.solver.regularization(5).rhs;
        end

        function p = assemble_solve_single(obj, M, b, Mp, regp)
            % Solve system of eqs
            [A, rhs, L] = assemble_single(obj, M, b, Mp, regp);
            p = obj.solve_single(A, rhs, L); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function p = assemble_solve_single_guess(obj, M, b, Mp, regp, p0)
            % Solve system of eqs
            [A, rhs, L] = assemble_single(obj, M, b, Mp, regp);
            p = obj.solve_single_guess(A, rhs, L, p0); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function p = solve_single(obj, A, rhs, L)
            % Solve system of eqs
            [p, ~, ~, ~] = pcg(A, rhs, obj.solver.opts.pcg_tol, obj.solver.opts.pcg_iter, L, L'); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function p = solve_single_guess(obj, A, rhs, L, p0)
            % Solve system of eqs
            [p, ~, ~, ~] = pcg(A, rhs, obj.solver.opts.pcg_tol, obj.solver.opts.pcg_iter, L, L', p0); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function rhs = b2rhs(obj, M, b, regp)
            rhs =  M'*b + regp(5)*obj.regularization.C{5}'*obj.regularization.rhs;
        end

        function par_idx = get_par_idx(obj, names_selection)
            par_idx = nan([1,numel(names_selection)]);
            for name_idx = 1:numel(names_selection)
                par_idx(name_idx) = find(strcmp(obj.solver.model.all_names, names_selection{name_idx}));
            end
        end

        function [time_in_column, vel_in_cell_tot] = depth_av_measuements(V, xs, mesh,cellID)% [time_in_column, vel_in_cell] = depth_av(cellID,V,xs,mesh,channel, direction)
         for dir = 1:3
                direct = ['U', 'V', 'W'];
                direction = direct(dir);
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

            bounds_cellID_x = [min(mesh.x_patch(:,cellID)) max(mesh.x_patch(:,cellID))];
            bounds_cellID_y = [min(mesh.y_patch(:,cellID)) max(mesh.y_patch(:,cellID))];
            bounds_cellID_z = [min(mesh.z_patch(:,cellID)) max(mesh.z_patch(:,cellID))];
        
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
            vel_in_cell_tot = NaN(size(vel_in_column, 3));
            for col = 1:size(vel_in_column, 2)
                valid_rows = (depth_in_column(:, col) >= bounds_cellID_z(1)) & (depth_in_column(:, col) <= bounds_cellID_z(2));
                % Extract the velocities for the valid rows in this column
                vel_in_cell(valid_rows, col) = vel_in_column(valid_rows, col);
                vel_in_cell_tot(:,direct) =  vel_in_cell(valid_rows, col);
            end
         end
        end

        function mod_names = modify_names(obj, names_selection)
            mod_names = names_selection;
            for idx = 1:numel(names_selection)
                mod_names{idx} = strrep(mod_names{idx}, 'sig', '\sigma');
                mod_names{idx} = strrep(mod_names{idx}, '^1', '');
                mod_names{idx} = strrep(mod_names{idx}, 'u0', 'u_0');
                mod_names{idx} = strrep(mod_names{idx}, 'v0', 'v_0');
                mod_names{idx} = strrep(mod_names{idx}, 'w0', 'w_0');

                mod_names{idx} = strrep(mod_names{idx}, 'ds', '\partial x');
                mod_names{idx} = strrep(mod_names{idx}, 'dn', '\partial y');
                mod_names{idx} = strrep(mod_names{idx}, 'd\sigma', '\partial \sigma');

                mod_names{idx} = strrep(mod_names{idx}, 'M0', 'M_0');
                mod_names{idx} = strrep(mod_names{idx}, 'M2', 'M_2');
                mod_names{idx} = strrep(mod_names{idx}, 'M4', 'M_4');

                mod_names{idx} = strrep(mod_names{idx}, 'A', '-A');
                mod_names{idx} = strrep(mod_names{idx}, '\phi', '-\phi');
                %mod_names{idx} = strrep(mod_names{idx}, 'M4', 'M_4');
            end
        end

        end
end
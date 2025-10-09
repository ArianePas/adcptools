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

        opts (1,1) SolverOptions

        name

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



        function CV = cross_validate_0D(obj)
            % simple cross-validation of the solution obj.p.
            CV = obj.cross_validate([obj.solver.regularization.weight]);
        end

        function CV = cross_validate_1D(obj, min, max, N)
            % 1D analysis: scalar min, max, N.
            %reg_pars_mat = reg_pars_symlog(obj, min, max)
            rp  = reg_pars_symlog(obj, min, max, N);
            reg_pars_mat = repmat(rp, 1, 5);
            CV = obj.cross_validate(reg_pars_mat);

            figure;
            plot(rp, [CV{:,1}])
            xlabel('reg pars')
            ylabel('generalization error')
            title('lambda vs scaled generalization error')
        end

        function CV = cross_validate_2D(obj, min, max, N)
            % 2D analysis: two-element inputs min, max, N. (i.e. min =
            % [0,0])
            reg_pars_cont = obj.reg_pars_symlog(min(1), max(1), N(1));
            reg_pars_smooth = obj.reg_pars_symlog(min(2), max(2), N(2));

            [rpc, rps] = meshgrid(reg_pars_cont, reg_pars_smooth);

            reg_pars = helpers.vecs2mat_full(reg_pars_cont, reg_pars_smooth);
            
            % paste in because of ordering of constraints.
            reg_pars_mat = [reg_pars(:,1), reg_pars(:,1), reg_pars(:,2), reg_pars(:,2), reg_pars(:,1)];
            
            CV = obj.cross_validate(reg_pars_mat);

            figure;
            % make use of the ordering of the constraints: first continuity
            
            contourf(helpers.symlog(rpc), helpers.symlog(rps), reshape([CV{:,1}]./CV{1,1}, N(2), N(1)), 100)
            colorbar
            % colormap()
            xlabel('cont (symlog10)')
            ylabel('smoothness (symlog10)')
            title('2D cross-validation')
            


        end

        


    end

    methods(Access=protected)



        function CV = cross_validate(obj, reg_pars_mat)
            % reg_pars_mat is a matrix of size nreg x 5, with the 5 known
            % regularization constraints.

            % This function is an ad hoc function and has not been
            % optimized using matrix algebra.
            p_train = cell([size(reg_pars_mat, 1), 1]);

            if strcmp(obj.opts.cv_mode, 'random')
                nepochs = obj.opts.cv_iter;
            else
                nepochs = 1;
            end
            niter = size(reg_pars_mat, 1)*nepochs;
            E = cell([size(reg_pars_mat, 1), 2]);
            CV = cell([size(reg_pars_mat, 1), 2]);
            i = 0;
            fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
            for ep = 1:nepochs % randomized iterations: Only meaningful if cv_mode == 'random'
                train_idx = logical(obj.split_dataset());
                test_idx = ~train_idx;

                % Construct training matrix and data
                M0 = obj.M(train_idx, :);
                b0 = obj.b(train_idx);

                M1 = obj.M(test_idx,:);
                b1 = obj.b(test_idx);

                Mp = M0'*M0;

                for rp = 1:size(reg_pars_mat, 1)
                    p_train{rp}(:, ep) = obj.assemble_solve_single(M0, b0, Mp, reg_pars_mat(rp,:));
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

        function training_idx = split_dataset(obj)
            %ci = vertcat(obj.cell_idx{:});
            ci = obj.cell_idx;
            training_idx = ones(size(ci));

            if strcmp(obj.opts.cv_mode, 'none')
                training_idx = ones(size(ci));
            elseif strcmp(obj.opts.cv_mode, 'random')
                tp = obj.opts.training_perc;
                rand0 = rand(size(training_idx));
                training_idx = (rand0 <= tp);
            elseif strcmp(obj.opts.cv_mode, 'omit_cells')
                oc = obj.opts.omit_cells;
                for occ = 1:length(oc)
                    training_idx(ci == oc(occ)) = 0;
                end
            elseif strcmp(obj.opts.cv_mode, 'kfold') % Disjoint partitions (in contrast to 'random')

            elseif strcmp(obj.opts.cv_mode, 'omit_time') % to be implemented
            end
        end

        function rp = reg_pars_symlog(obj, min, max, N)%, res_near_zero = 1.0)
            %Variation upon logspace.
            ymin = helpers.symlog(min);
            ymax = helpers.symlog(max);
            rp = helpers.symexp(linspace(ymin, ymax, N))';
        end


        function reg_pars_sens_vec = vectorize_reg_pars(obj, cont, smooth)
            if strcmp(obj.opts.reg_vary, 'coupled')
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
            obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
            L = ichol(A, obj.opts.preconditioner_opts);
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
            [p, ~, ~, ~] = pcg(A, rhs, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L'); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function p = solve_single_guess(obj, A, rhs, L, p0)
            % Solve system of eqs
            [p, ~, ~, ~] = pcg(A, rhs, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L', p0); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
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
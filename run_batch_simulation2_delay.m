clc;
clear;
close all;
addpath functions;

rng shuffle;
%% 1. DEFINE MODEL LIBRARIES
% Creation: Poisson process
creation_meantime_values = 0.01;

% Degradation: two-channel-delay
rate1_values = logspace(-2,0,3);          % channel 1 decay rate (during delay)
rate2_values = 1;          % channel 2 decay rate (after delay)
delay_meantime_value = 20; % mean delay time
delay_relvar_values = [0, 0.01, 0.1, 0.5, 1];   % delay CV^2 (0 = deterministic, >0 = Gamma)

[creation_models, degradation_models] = define_models_("Poisson_Delay", ...
    creation_meantime_values = creation_meantime_values, ...
    rate1_values = rate1_values, ...
    rate2_values = rate2_values, ...
    delay_meantime_values = delay_meantime_value, ...
    delay_relvar_values = delay_relvar_values);
%% 2. SET SIMULATION PARAMETERS
simulation_parameters = define_parameters_(creation_models, degradation_models,...
    numTrajectories = 1e2, ...
    T0_over_meantime_c = [-10, -10], ...
    dt = delay_meantime_value*1e-2, ...
    Tmax = delay_meantime_value*250, ...
    ss_time_over_meantime_d = 50);
%% set directory to save data
[results_dir, stt_time] = set_directory_();
%% 3. BATCH EXECUTION (NESTED LOOPS)
clc
end_time = batch_execution_(simulation_parameters, results_dir, stt_time = stt_time);
%% functions
function [creation_models, degradation_models] = define_models_(ver, options)
    arguments
        ver
        options.creation_relvar_values = 1;
        options.creation_meantime_values = 1;
        options.degradation_relvar_values = 1;
        options.degradation_meantime_values = 1;
        options.rate1_values = 1;
        options.rate2_values = 1;
        options.delay_meantime_values = 1;
        options.delay_relvar_values = 0;
    end
    creation_models = {};
    degradation_models = {};

    switch ver
        case {"Poisson_Delay"}
            % Creation: Poisson models
            creation_meantime_values = options.creation_meantime_values;
            for ii = 1:numel(creation_meantime_values)
                meantc_ = creation_meantime_values(ii);
                parameters = struct('rate', 1/meantc_);
                new_model = struct("name", "Poisson", ...
                    'parameters', parameters, ...
                    'meantime', meantc_, ...
                    'detail_name', sprintf("E(r=%g)", 1/meantc_));
                creation_models = [creation_models, new_model]; %#ok<*AGROW>
            end
            % Degradation: two-channel-delay models
            % 4D grid: rate1 x rate2 x delay_meantime x delay_relvar
            [r1g, r2g, dmg, drg] = ndgrid(options.rate1_values, options.rate2_values, ...
                options.delay_meantime_values, options.delay_relvar_values);
            for ii = 1:numel(r1g)
                rate1 = r1g(ii);
                rate2 = r2g(ii);
                delay_mean = dmg(ii);
                delay_rv = drg(ii);
                % Construct delay model based on relvar
                if delay_rv == 0
                    delay_model = struct('name', 'deterministic', ...
                        'parameters', struct('tau', delay_mean));
                    delay_str = sprintf("Det(%g)", delay_mean);
                else
                    delay_model = struct('name', 'Gamma', ...
                        'parameters', struct('alpha', 1/delay_rv, 'beta', delay_mean*delay_rv));
                    delay_str = sprintf("G(a=%g,b=%g)", 1/delay_rv, delay_mean*delay_rv);
                end
                params = struct('rate1', rate1, ...
                                'rate2', rate2, ...
                                'delay_model', delay_model);
                new_models = struct('name', 'two-channel-delay', ...
                                    'parameters', params, ...
                                    'detail_name', sprintf("TCD(r1=%g,r2=%g,d=%s)", rate1, rate2, delay_str));
                new_models.meantime = get_meantime(new_models);
                degradation_models = [degradation_models, new_models];
            end
        case "Models"
            % define creation models
            [relvarc, meantc] = ndgrid(options.creation_relvar_values, options.creation_meantime_values);
            [relvard, meantd] = ndgrid(options.degradation_relvar_values, options.degradation_meantime_values);
            for ii = 1:numel(meantc)
                meantc_ = meantc(ii);
                relvarc_ = relvarc(ii);
                if relvarc_ > 1
                    % define 3 channel creation process
                    % let k3 = 0, k13 = 0, k31 = 0, k12 = k23 = kd, k21 = k32 = ku for simplicity
                    % set k2 = 1/meantc_ and r = kd/ku = 1 for additional simplilcity.
                    k2 = 1/meantc_; % do not change this; If k2 is not given by 1/meantc_, k1 and ku should be defined differently.
                    r = 1; % You can change the value of r
                    k1 = (1 + r^2)/meantc_;
                    ku = 2*r^2*(1 + r)/(meantc_*(1 + r + r^2)*(relvarc_ - 1));
                    kd = r.*ku;
                    parameters = struct(transition_rates = [0, kd, 0; ku, 0, kd; 0, ku, 0], ...
                        reaction_rates = [k1, k2, 0], ...
                        initial_probability = [1/3, 1/3, 1/3]);
                    new_models = struct(name = "multi-channel", ...
                        parameters = parameters, ...
                        meantime = meantc_, ...
                        detail_name = sprintf("3C(%g,%g,%g,%g,%g)",kd,ku,k1,k2));
                elseif relvarc_ > 0
                    gamma_params = struct('alpha',1/relvarc_,'beta',meantc_*relvarc_);
                    new_models = define_models("Gamma", gamma_params);
                elseif relvarc_ == 0
                    new_models = struct(name = "deterministic", ...
                        parameters = struct('tau', meantc_), ...
                        meantime = meantc_, ...
                        detail_name = sprintf("Det(%g)",meantc_));
                end
                creation_models = [creation_models, new_models]; %#ok<*AGROW>
            end
            for ii = 1:numel(meantd)
                meantd_ = meantd(ii);
                relvard_ = relvard(ii);
                if relvard_ > 1
                    k21 = 0; k12 = 1/meantd_;
                    lambda = sqrt(-7 + 10*relvard_ + relvard_^2);
                    gamma1 = (relvard_ + 1 + lambda)/(4*meantd_);
                    gamma2 = k12*(1 + relvard_ - lambda)/(2 - 2*relvard_);
                    parameters = struct(transition_rates = [0, k12; k21, 0], ...
                        reaction_rates = [gamma1, gamma2], ...
                        initial_probability = [1, 0]);
                    new_models = struct(name = "multi-state", ...
                        parameters = parameters, ...
                        meantime = meantd_, ...
                        detail_name = sprintf("2S(%g,%g,%g,%g)",k12,k21,gamma1,gamma2));
                elseif relvard_ > 0
                    gamma_params = struct(alpha = 1/relvard_,beta = meantd_*relvard_);
                    new_models = define_models("Gamma", gamma_params);
                elseif relvard_ == 0
                    new_models = struct(name = "deterministic", ...
                        parameters = struct('tau', meantd_), ...
                        meantime = meantd_, ...
                        detail_name = sprintf("Det(%g)",meantd_));
                end
                degradation_models = [degradation_models, new_models];
            end
        case {"G", "Gamma"}
            % Generate Gamma models
            [relvarc, meantc] = ndgrid(options.creation_relvar_values, options.creation_meantime_values);
            [relvard, meantd] = ndgrid(options.degradation_relvar_values, options.degradation_meantime_values);
            alpha_vals = 1./relvarc(:);
            beta_vals = meantc(:) ./ alpha_vals;
            gamma_params = arrayfun(@(a,b) struct('alpha',a,'beta',b), alpha_vals, beta_vals);
            creation_models = [creation_models, define_models("Gamma", gamma_params)];
            
            alpha_vals = 1./relvard(:);
            beta_vals = meantd(:) ./ alpha_vals;
            gamma_params_d = arrayfun(@(a,b) struct('alpha',a,'beta',b), alpha_vals, beta_vals);
            degradation_models = [degradation_models, define_models("Gamma", gamma_params_d)];
    end
end
function simulation_parameters = define_parameters_(creation_models, degradation_models, options)
    arguments
        creation_models
        degradation_models
        options.numTrajectories = 1;
        options.T0_over_meantime_c = 0;
        options.dt_over_meantime_d;
        options.dt
        options.Tmax_over_meantime_d;
        options.Tmax
        options.ss_time_over_meantime_d = 10;
    end
    simulation_parameters = options;
    simulation_parameters.creation_models = creation_models;
    simulation_parameters.degradation_models = degradation_models;
    fprintf('Simulation parameters are set.\n');
end
function [results_dir, stt_time] = set_directory_()
    stt_time = datetime("now");
    data_dir = uigetdir(pwd, 'Select a directory to save the data');
    if data_dir == 0
        disp('Folder selection cancelled by user. Aborting.');
        return;
    end
    results_dir = sprintf("simulations_%s", datetime("now", Format = "yyyyMMddHHmmss"));
    results_dir = fullfile(data_dir, results_dir);
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
end
function end_time = batch_execution_(simulation_parameters, results_dir, options)
    arguments
        simulation_parameters
        results_dir
        options.stt_time = datetime("now")
    end
    stt_time = options.stt_time;
    creation_models = simulation_parameters.creation_models;
    degradation_models = simulation_parameters.degradation_models;
    
    parameters_filename = fullfile(results_dir, 'simulation_parameters.mat');
    fprintf('\nSaving simulation parameters to %s\n', parameters_filename);
    save(parameters_filename, 'simulation_parameters');
    
    fprintf('Total simulations to run: %d (%d creation x %d degradation models)\n', ...
        length(creation_models) * length(degradation_models), length(creation_models), length(degradation_models));
    p = gcp("nocreate");
    if isempty(p)
        gcp
    end
    run_index = 0;
    for ii = 1:length(creation_models)
        for jj = 1:length(degradation_models)
            run_index = run_index + 1;
            % --- Get models for the current run ---
            creation_model = creation_models{ii};
            degradation_model = degradation_models{jj};
            
            % --- Create a unique, model-aware name for this run ---
            create_str = creation_model.detail_name;
            degrade_str = degradation_model.detail_name;
            run_name = sprintf("Run%04d_Create-%s_Degrade-%s", run_index, create_str, degrade_str);
            fprintf('\n--- Starting: %s ---\n', run_name);
            
            % --- Run Simulation ---
            fprintf('Running simulation...\n');
            tic;
            [moments, ~, psd_avg, omega, acf, lags] = run_simulation2( ...
                creation_model, degradation_model, simulation_parameters);
            sim_runtime = toc;
            % --- Run Analysis ---
            fprintf('Analyzing results...\n');
            
            tot_runtime = toc;
            % --- Save Results ---
            results_filename = fullfile(results_dir,run_name + ".mat");
            fprintf('Saving results to %s\n', results_filename);
            P_avg_rect = psd_avg(1,:);
            P_avg_hann = psd_avg(2,:);
            P_avg_hamm = psd_avg(3,:);
            P_avg_bkmn = psd_avg(4,:);
            P_avg_kisr = psd_avg(5,:);

            save(results_filename, 'P_avg_rect', 'P_avg_hann', 'P_avg_hamm', 'P_avg_bkmn', 'P_avg_kisr', 'moments', ...
                                    'omega', 'lags', 'acf', 'creation_model', 'degradation_model','sim_runtime','tot_runtime');
        end
    end
    
    fprintf('\nBatch simulation complete. Results and plots saved in %s folder.\n', results_dir);
    for ii = 1:100
        fprintf("=")
    end
    
    end_time = datetime("now");
    simulation_parameters.elapsed_time = end_time - stt_time;
    simulation_parameters.cpu_info = getenv('PROCESSOR_IDENTIFIER');
    simulation_parameters.matlab_version = version;
    
    parameters_filename = fullfile(results_dir, 'simulation_parameters.mat');
    fprintf('\nSaving simulation parameters to %s\n', parameters_filename);
    save(parameters_filename, 'simulation_parameters');
end
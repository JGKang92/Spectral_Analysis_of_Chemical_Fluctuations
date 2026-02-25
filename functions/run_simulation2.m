function [moments, time, psd_avg, omega, acf_avg, lags] = run_simulation2(creation_model, degradation_model, simulation_parameters)
    arguments
        creation_model
        degradation_model
        simulation_parameters
    end
    numTrajectories = simulation_parameters.numTrajectories;
    meantime_c = get_meantime(creation_model);
    meantime_d = get_meantime(degradation_model);
    T0 = simulation_parameters.T0_over_meantime_c * meantime_c;
    if isfield(simulation_parameters, "dt")
        dt = simulation_parameters.dt;
    else
        dt = simulation_parameters.dt_over_meantime_d * meantime_d;
    end
    if isfield(simulation_parameters, "Tmax")
        Tmax = simulation_parameters.Tmax;
    else
        Tmax = simulation_parameters.Tmax_over_meantime_d * meantime_d;
    end
    ss_time = simulation_parameters.ss_time_over_meantime_d * meantime_d;
    ss_idx = ceil(ss_time / dt);
    fs = 1/dt;
    
    lent = ceil(Tmax./dt);
    edges = -0.5 + (1:1:lent+1);
    moments = zeros(4, lent);
    acf_avg = zeros(1, lent-ss_idx+1);
    lenw = floor((lent-ss_idx+1)/2)+1;
    psd_avg = zeros(5, lenw);

    parfor ii = 1:numTrajectories
        tcList = Generate_tcList(T0, Tmax, creation_model);
        if isempty(tcList)
            lifetimes = inf(size(tcList));
        else
            lifetimes = generate_random_variates(degradation_model, numel(tcList));
        end
        tdList = tcList + lifetimes;
        
        icList = ceil(tcList/dt);
        idList = ceil(tdList/dt);

        n_Birth = histcounts(icList, edges);
        n_Death = histcounts(idList, edges);

        number = cumsum(n_Birth - n_Death);
        number_ss = number(ss_idx:end);
        P_avg_rect = calculate_psd(number_ss, fs, "rectangular");
        % P_avg_rect = nan(1, lenw);
        P_avg_hann = calculate_psd(number_ss, fs, "hann");
        % P_avg_hann = nan(1, lenw);
        P_avg_hamm = calculate_psd(number_ss, fs, "hamming");
        % P_avg_hamm = nan(1, lenw);
        P_avg_bkmn = calculate_psd(number_ss, fs, "blackman");
        % P_avg_bkmn = nan(1, lenw);
        P_avg_kisr = calculate_psd(number_ss, fs, "kaiser", 5);
        % P_avg_kisr = nan(1, lenw);
        psd_avg = psd_avg + [P_avg_rect; P_avg_hann; P_avg_hamm; P_avg_bkmn; P_avg_kisr]./numTrajectories;
        
        acf = number_ss.*number_ss(1);
        acf_avg = acf_avg + acf./numTrajectories;
        
        moments = moments + [number; number.^2; number.^3; number.^4]./numTrajectories;
    end
    
    time = (1:lent) * dt;
    acf_avg = acf_avg - moments(1,ss_idx)*moments(1,ss_idx:end);
    omega = 2*pi*fs*(0:((lent-ss_idx+1)/2))/(lent-ss_idx+1);
    lags = (0:lent-ss_idx) * dt;
    fprintf("simulation finished numTrajectories = %d\n", numTrajectories)
end
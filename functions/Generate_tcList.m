function tcList = Generate_tcList(T0, T_max, model)
    arguments
        T0
        T_max
        model
    end
    rng shuffle
    if isscalar(T0)
        T = T0;
    elseif isnumeric(T0)
        T = T0(1) + T0(2)*rand();
        % fprintf("T0 is set to %f\n",T)
    end
    tcList = [];

    if strcmp(model.name, "multi-channel")
        parameters = model.parameters;
        K_transition = parameters.transition_rates;
        k_reaction = parameters.reaction_rates;
        initial_probability = parameters.initial_probability;
        channel = randsample(1:length(initial_probability), 1, true, initial_probability);
        while T < T_max
            % Determine the time to the next channel transition
            rates_from_current_state = K_transition(channel, :);
            [channel_time, next_channel] = min(arrayfun(@(x)random("Exponential", 1/x),rates_from_current_state));
            T_prev = T;
            while true
                reaction_time = random("Exponential", 1/k_reaction(channel));
                T = T + reaction_time;
                if T > T_prev + channel_time
                    T = T_prev + channel_time;
                    channel = next_channel;
                    break
                elseif T < T_max && T > 0
                    tcList = [tcList, T];
                end
            end
        end
    else
        while T < T_max
            T = T + generate_random_variates(model, 1);
            if T < T_max && T > 0
                tcList = [tcList, T];
            end
        end
    end
    tcList = tcList';
end
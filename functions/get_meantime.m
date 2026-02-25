function [meantime, model] = get_meantime(model)
    parameters = model.parameters;
    if isfield(model, "meantime")
        meantime = model.meantime;
    else
        switch model.name
            case {'deterministic'}
                meantime = parameters.tau;
            case {'Poisson', 'Exponential', '1-state'}
                meantime = 1./parameters.rate;
            case 'Gamma'
                meantime = parameters.alpha * parameters.beta;
            case "multi-state"
                K_transition = parameters.transition_rates;
                k_reaction = parameters.reaction_rates;
                initial_probability = parameters.initial_probability;
                switch numel(initial_probability)
                    case 2
                        k12 = K_transition(1, 2);
                        k21 = K_transition(2, 1);
                        gamma1 = k_reaction(1);
                        gamma2 = k_reaction(2);
                        P1init = initial_probability(1);
                        P2init = initial_probability(2);
                        numer = k12 + k21 + P1init*gamma1 + P2init*gamma2;
                        denom = k12*gamma2 + k21*gamma1 + gamma1*gamma2;
                        meantime = numer / denom;
                    case 3
                        error("Not implemented yet.")
                end
            
            case "multi-channel"
                switch numel(parameters.reaction_rates)
                    case 2
                        K_transition = parameters.transition_rates;
                        k_reaction = parameters.reaction_rates;
                        time_scale = 1./(K_transition(1, 2) + K_transition(2, 1));
                        P1ss = K_transition(2, 1)*time_scale;
                        P2ss = K_transition(1, 2)*time_scale;
                        if P1ss+P2ss ~= 1
                            error('Steady-state probabilities do not sum to 1.');
                        end
                        meantime = 1 / (k_reaction(1)*P1ss + k_reaction(2)*P2ss);
                    case 3
                        K_transition = parameters.transition_rates;
                        k_reaction = parameters.reaction_rates;
                        k1 = k_reaction(1);
                        k2 = k_reaction(2);
                        k3 = k_reaction(3);
                        k12 = K_transition(1, 2);
                        k13 = K_transition(1, 3);
                        k21 = K_transition(2, 1);
                        k23 = K_transition(2, 3);
                        k31 = K_transition(3, 1);
                        k32 = K_transition(3, 2);

                        numer = k21*k31 + k23*k31 + k21*k32 + k13*(k21 + k23 + k32) + k12*(k23 + k31 + k32);
                        denom = k13*(k21*k3 + k23*k3 + k2*k32) + k12*(k23*k3 + k2*k31 + k2*k32) + k1*(k21*k31 + k23*k31 + k21*k32);
                        meantime = numer / denom;
                end

            case 'two-channel-delay'
                rate1 = parameters.rate1;
                rate2 = parameters.rate2;
                delay_model = parameters.delay_model;
                P_D_g1 = calculate_characteristic_function(delay_model, rate1);
                meantime = 1/rate1 - P_D_g1 * (rate2 - rate1) / (rate1 * rate2);
        end
    end
    model.meantime = meantime;
end
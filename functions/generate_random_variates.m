function variates = generate_random_variates(model, numSamples)
%GENERATE_RANDOM_VARIATES Generates random variates from a specified model.
%
%   variates = GENERATE_RANDOM_VARIATES(model, numSamples) generates
%   numSamples random variates based on the statistical model defined in
%   the 'model' struct.

    arguments
        model
        numSamples (1,1) {mustBeInteger, mustBePositive} = 1
    end
    parameters = model.parameters;
    switch model.name
        case {'Poisson', 'Exponential', '1-state'}
            rate = parameters.rate;
            variates = random('Exponential', 1/rate, numSamples, 1);
        case "multi-state"
            K_transition = parameters.transition_rates;
            k_reaction = parameters.reaction_rates;
            initial_probability = parameters.initial_probability;
            variates = nan(numSamples, 1);
            switch numel(initial_probability)
                case 2
                    for ii = 1:numSamples
                        variates(ii) = generate_random_variates_2state(K_transition, k_reaction, initial_probability);
                    end
                case 3
                    error("Not implemented yet.")
            end
        case 'multi-step'
            k = parameters.rate;
            n = parameters.numSteps;
            variates = random('Gamma', n, 1/k, numSamples, 1);

        case 'Gamma'
            a = parameters.alpha;
            b = parameters.beta;
            variates = random('Gamma', a, b, numSamples, 1);
        case 'deterministic'
            tau = parameters.tau;
            variates = tau * ones(numSamples, 1);

        case 'two-channel-delay'
            rate1 = parameters.rate1;
            rate2 = parameters.rate2;
            delay_model = parameters.delay_model;
            T1 = random('Exponential', 1/rate1, numSamples, 1);
            D = generate_random_variates(delay_model, numSamples);
            T2 = D + random('Exponential', 1/rate2, numSamples, 1);
            variates = T1;
            idx = T1 >= D;
            variates(idx) = T2(idx);

        otherwise
            error('Unknown model name: %s', model.name);
    end
end
function time = generate_random_variates_2state(K_transition, k_reaction, initial_probability)
    state = randsample(1:length(initial_probability), 1, true, initial_probability);
    time = 0;
    while true
        state_transition_time = random('Exponential', 1./K_transition(state, :));
        reaction_time = random('Exponential', 1./k_reaction(state));
        [val, idx] = min([reaction_time, state_transition_time]);
        time = time + val;
        if idx == 1
            break
        else
            state = idx - 1; % switch between 1 and 2
        end
    end
end
function P_s = calculate_characteristic_function(model, s)
%CALCULATE_CHARACTERISTIC_FUNCTION Computes the characteristic function P(s).
%
%   P_s = CALCULATE_CHARACTERISTIC_FUNCTION(model, s) calculates the
%   characteristic function (Laplace transform of the PDF) for a given
%   lifetime distribution model, evaluated at the Laplace variable values s.
%
%   Inputs:
%   model - A struct defining the lifetime distribution model.
%           It must have a .name and a .parameters field.
%   s     - A vector or scalar of Laplace variable values (can be complex).
%
%   Outputs:
%   P_s   - The calculated characteristic function, P(s) = L{PDF}(s).

    % Ensure s is a column vector for consistent indexing
    s = s(:).';
    
    % Initialize output
    P_s = zeros(size(s));
    
    % Find zero and non-zero s indices
    zero_idx = (s == 0);
    nonzero_idx = (s ~= 0);
    s_nonzero = s(nonzero_idx);
    parameters = model.parameters;

    % Calculate P(s) for non-zero s
    switch model.name
        case {'deterministic'}
            tau = parameters.tau;
            P_s(nonzero_idx) = exp(-s_nonzero*tau);

        case {'Poisson', 'Exponential', '1-state'}
            r = parameters.rate;
            P_s(nonzero_idx) = r ./ (s_nonzero + r);

        case 'multi-step'
            k = parameters.rate;
            n = parameters.numSteps;
            P_s(nonzero_idx) = (k ./ (s_nonzero + k)).^n;

        case 'Gamma'
            a = parameters.alpha; % shape
            b = parameters.beta;  % scale
            P_s(nonzero_idx) = (1 + b.*s_nonzero).^(-a);

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
                    numer = k21*gamma1 + k12*gamma2 + gamma1*gamma2 + s_nonzero.*(P1init*gamma1 + P2init*gamma2);
                    denom = k12*(s_nonzero + gamma2) + (s_nonzero + gamma1).*(s_nonzero + k21 + gamma2);
                otherwise
                    error("Not implemented yet.")
            end
            P_s(nonzero_idx) = numer ./ denom;
        case 'two-channel-delay'
            rate1 = parameters.rate1;
            rate2 = parameters.rate2;
            delay_model = parameters.delay_model;
            % Claude's code
            P_D_shifted = calculate_characteristic_function(delay_model, s_nonzero + rate1);
            P_s(nonzero_idx) = rate1 ./ (s_nonzero + rate1) ...
                + P_D_shifted .* s_nonzero .* (rate2 - rate1) ...
                ./ ((s_nonzero + rate1) .* (s_nonzero + rate2));

        case {'multi-channel'}
            meantime = get_meantime(model);
            transition_rates = parameters.transition_rates;
            reaction_rates = parameters.reaction_rates;
            switch numel(parameters.reaction_rates)
                case 2
                    k1 = reaction_rates(1);
                    k2 = reaction_rates(2);
                    k12 = transition_rates(1, 2);
                    k21 = transition_rates(2, 1);
                    P1ss = k21 / (k12 + k21);
                    P2ss = k12 / (k12 + k21);
                    if P1ss+P2ss ~= 1
                        warning('Steady-state probabilities do not sum to 1.');
                    end
                    CR_Laplace = P1ss*P2ss*(k1 - k2)^2 / (s + k12 + k21);
                case 3
                    k1 = reaction_rates(1);
                    k2 = reaction_rates(2);
                    k3 = reaction_rates(3);
                    k12 = transition_rates(1, 2);
                    k13 = transition_rates(1, 3);
                    k21 = transition_rates(2, 1);
                    k23 = transition_rates(2, 3);
                    k31 = transition_rates(3, 1);
                    k32 = transition_rates(3, 2);
                    A1 = k12 + k13 + k21 + k23 + k31 + k32;
                    A2 = k13*k21 + k12*k23 + k13*k23 + k12*k31 + k21*k31 + k23*k31 + k12*k32 + k13*k32 + k21*k32;
                    P1ss = (k21*k31 + k23*k31 + k21*k32)/A2;
                    P2ss = (k12*k31 + k12*k32 + k13*k32)/A2;
                    P3ss = (k13*k21 + k12*k23 + k13*k23)/A2;
                    if P1ss+P2ss+P3ss ~= 1
                        warning('Steady-state probabilities do not sum to 1.');
                    end
                    a12 = -(k12^2*k23*(2*k21*(k31 + k32) - k31*(-k23 + k31 + k32))) - ...
                            k13*k32*(k13*(k21 + k23)*(k21 - 2*k31) - k13*k21*k32 - (k21*k31 + k23*k31 + k21*k32)*(k21 + 2*(k23 + k31 + k32))) - ...
                            k12*(k13*(k23*k31*(k23 - 2*k31 - 3*k32) + 2*k21^2*(k31 + k32) + k21*(3*k23 - 2*k31 - k32)*(k31 + k32)) - ...
                            (k21*k31 + k23*k31 + k21*k32)*(2*(k31 + k32)^2 + k23*(k31 + 2*k32)));
                    a13 = k13*(2*(k21 + k23)^3*k31 + (k21 + k23)*(2*k21*(k21 + k23) + k13*(k21 - 2*k31) + (k21 + 2*k23)*k31)*k32 + k21*(-k13 + k21 + 2*k23)*k32^2) + ...
                            k12^2*k23*(2*k21*(k31 + k32) - k31*(-k23 + k31 + k32)) + ...
                            k12*(k23*(2*k21 + 2*k23 + k31 + 2*k32)*(k23*k31 + k21*(k31 + k32)) + ...
                            k13*(k23*k31*(k23 - 2*k31 - 3*k32) + 2*k21^2*(k31 + k32) + k21*(3*k23 - 2*k31 - k32)*(k31 + k32)));
                    a23 = 2*k12^3*k23*(k31 + k32) + k13*k32*(2*k13^2*(k21 + k23) + k13*(k21 + k23)*(k21 + 2*k31) + k13*k21*k32 - (k21 + 2*k23)*(k23*k31 + k21*(k31 + k32))) + ...
                            k12^2*(2*k21*k23*(k31 + k32) + k23*k31*(k23 + k31 + k32) + 2*k13*(k21*(k31 + k32) + k23*(2*k31 + 3*k32))) + ...
                            k12*(-(k23*(k31 + 2*k32)*(k23*k31 + k21*(k31 + k32))) + 2*k13^2*(k21*(k31 + 2*k32) + k23*(k31 + 3*k32)) + ...
                            k13*(2*k21^2*(k31 + k32) + k21*(k31 + k32)*(3*k23 + 2*k31 + k32) + k23*k31*(k23 + 2*k31 + 3*k32)));
                    numer0 = (a12*(k1 - k2)^2 + a13*(k1 - k3)^2 + a23*(k2 - k3)^2)/A2^2;
                    numer1 = P1ss*P2ss*(k1 - k2)^2 + P1ss*P3ss*(k1 - k3)^2 + P2ss*P3ss*(k2 - k3)^2;
                    denom = A2 + A1*s + s.^2;
                    CR_Laplace = (numer0 + numer1*s)/denom;
                otherwise
                    error('Multi-channel models with more than 3 states are not supported.');
            end
            P_s(nonzero_idx) = CR_Laplace*meantime/(1 + CR_Laplace*meantime);
        otherwise
            error('Unknown model name: %s', model.name);
    end

    % For s=0, P(0) is the integral of the PDF, which is always 1.
    if any(zero_idx)
        P_s(zero_idx) = 1;
    end
    
end
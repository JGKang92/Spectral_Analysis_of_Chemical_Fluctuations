function S_rate = calculate_rate_psd(model, omega)
%CALCULATE_RATE_PSD Calculates the PSD of rate fluctuations.
%
%   S_rate = CALCULATE_RATE_PSD(model, omega) calculates the power
%   spectral density of rate fluctuations for a given lifetime distribution
%   model. This is defined as the real part of the characteristic function
%   (Laplace transform of the PDF) evaluated at s = i*omega.
%
%   Inputs:
%   model - A struct defining the lifetime distribution model.
%           It must have a .name and a .parameters field.
%   omega - A vector of angular frequencies (rad/s) at which to
%           calculate the spectrum.
%
%   Outputs:
%   S_rate - The calculated power spectral density of rate fluctuations.

    % Define Laplace variable
    s = 1i * omega(:).';
    meantime = get_meantime(model);

    % Calculate P(s) using the dedicated function
    if model.name == "multi-channel"
        transition_rates = model.parameters.transition_rates;
        reaction_rates = model.parameters.reaction_rates;
        switch numel(model.parameters.reaction_rates)
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
                CR_Laplace = P1ss*P2ss*(k1 - k2)^2 ./ (s + k12 + k21);
                S_rate = 2*real(CR_Laplace);
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
                a12 = -(k12^2*k23*(2*k21*(k31 + k32) - k31*(-k23 + k31 + k32))) ...
                    - k13*k32*(k13*(k21 + k23)*(k21 - 2*k31) - k13*k21*k32 - (k21*k31 + k23*k31 + k21*k32)*(k21 + 2*(k23 + k31 + k32))) ...
                    - k12*(k13*(k23*k31*(k23 - 2*k31 - 3*k32) + 2*k21^2*(k31 + k32) + k21*(3*k23 - 2*k31 - k32)*(k31 + k32)) ...
                    - (k21*k31 + k23*k31 + k21*k32)*(2*(k31 + k32)^2 + k23*(k31 + 2*k32)));
                a13 = k13*(2*(k21 + k23)^3*k31 + (k21 + k23)*(2*k21*(k21 + k23) + k13*(k21 - 2*k31) + (k21 + 2*k23)*k31)*k32 + k21*(-k13 + k21 + 2*k23)*k32^2) ...
                    + k12^2*k23*(2*k21*(k31 + k32) - k31*(-k23 + k31 + k32)) ...
                    + k12*(k23*(2*k21 + 2*k23 + k31 + 2*k32)*(k23*k31 + k21*(k31 + k32))...
                    + k13*(k23*k31*(k23 - 2*k31 - 3*k32) + 2*k21^2*(k31 + k32) + k21*(3*k23 - 2*k31 - k32)*(k31 + k32)));
                a23 = 2*k12^3*k23*(k31 + k32) + k13*k32*(2*k13^2*(k21 + k23) + k13*(k21 + k23)*(k21 + 2*k31) + k13*k21*k32 - (k21 + 2*k23)*(k23*k31 + k21*(k31 + k32))) + ...
                        k12^2*(2*k21*k23*(k31 + k32) + k23*k31*(k23 + k31 + k32) + 2*k13*(k21*(k31 + k32) + k23*(2*k31 + 3*k32))) + ...
                        k12*(-(k23*(k31 + 2*k32)*(k23*k31 + k21*(k31 + k32))) + 2*k13^2*(k21*(k31 + 2*k32) + k23*(k31 + 3*k32)) + ...
                        k13*(2*k21^2*(k31 + k32) + k21*(k31 + k32)*(3*k23 + 2*k31 + k32) + k23*k31*(k23 + 2*k31 + 3*k32)));
                numer0 = (1/2)*(a12*(k1 - k2)^2 + a13*(k1 - k3)^2 + a23*(k2 - k3)^2)/A2^2;
                numer1 = (P1ss*P2ss*(k1 - k2)^2 + P1ss*P3ss*(k1 - k3)^2 + P2ss*P3ss*(k2 - k3)^2);
                denom = A2 + A1.*s + s.^2;
                CR_Laplace = (numer0 + numer1.*s)./denom;
                S_rate = 2*real(CR_Laplace);
            otherwise
                error('Multi-channel models with more than 3 states are not supported.');
        end
    elseif ismember(string(model.name), ["Poisson", "Exponential","1-state"])
        S_rate = zeros(size(s));
    else
        P_s = calculate_characteristic_function(model, s);
        
        % The PSD of rate fluctuations is the real part of P(s)
        S_rate = 2*real(P_s./(1 - P_s))./meantime;
    end
end
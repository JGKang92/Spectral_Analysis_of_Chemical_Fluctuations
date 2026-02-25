function chi = calculate_susceptibility(degradation_model, omega)
%CALCULATE_SUSCEPTIBILITY Calculates the susceptibility of a degradation model.
%
%   chi = CALCULATE_SUSCEPTIBILITY(degradation_model, omega) calculates the
%   susceptibility, defined as the squared absolute value of the Laplace
%   transform of the lifetime SURVIVAL PROBABILITY, evaluated at s = i*omega.
%
%   Inputs:
%   degradation_model - A struct defining the lifetime distribution model.
%                       It must have a .name and a .parameters field.
%   omega             - A vector of angular frequencies (rad/s) at which to
%                       calculate the susceptibility.
%
%   Outputs:
%   chi               - The calculated susceptibility, |L{S(t)}(i*omega)|^2.

    % Ensure omega is a column vector for indexing
    omega = omega(:)';
    s = 1i * omega;
    
    % Initialize output
    chi = zeros(size(omega));
    
    % Find zero and non-zero omega indices
    zero_idx = (omega == 0);
    nonzero_idx = (omega ~= 0);
    s_nonzero = s(nonzero_idx);

    % Calculate P(s) using the dedicated function
    P_s = calculate_characteristic_function(degradation_model, s);

    % Calculate mean lifetime for the s=0 case
    if isfield(degradation_model, "meantime")
        meantime = degradation_model.meantime;
    else
        meantime = get_meantime(degradation_model);
    end

    % For omega=0, s=0. L{S(t)}(0) is the mean lifetime.
    if any(zero_idx)
        chi(zero_idx) = meantime^2;
    end

    % For omega~=0, s~=0. L{S(t)}(s) = (1 - P(s))/s
    if any(nonzero_idx)
        S_L_s = (1 - P_s(nonzero_idx)) ./ s_nonzero;
        chi(nonzero_idx) = abs(S_L_s).^2;
    end

end
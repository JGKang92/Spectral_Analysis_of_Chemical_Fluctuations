function Sz0 = calculate_basal_spectrum(mean_creation_rate, degradation_model, omega)
%CALCULATE_BASAL_SPECTRUM Calculates the basal spectrum of a degradation model.
%
%   Sz0 = CALCULATE_BASAL_SPECTRUM(degradation_model, omega) calculates the
%   basal spectrum, defined as the imaginary part of the Laplace
%   transform of the lifetime SURVIVAL PROBABILITY, evaluated at s = i*omega.
%
%   Inputs:
%   degradation_model - A struct defining the lifetime distribution model.
%                       It must have a .name and a .parameters field.
%   omega             - A vector of angular frequencies (rad/s) at which to
%                       calculate the spectrum.
%
%   Outputs:
%   Sz0             - The calculated basal spectrum, Im(L{S(t)}(i*omega)).

    % Ensure omega is a column vector for indexing
    omega = omega(:)';
    s = 1i * omega;
    
    % Initialize output
    Sz0 = zeros(size(omega));
    
    % Find zero and non-zero omega indices
    zero_idx = (omega == 0);
    nonzero_idx = (omega ~= 0);
    s_nonzero = s(nonzero_idx);

    % Calculate P(s), the Laplace transform of the PDF
    P_s = calculate_characteristic_function(degradation_model, s);

    % For omega=0, s=0. L{S(t)}(0) is the mean lifetime (real), so Im is 0.
    if any(zero_idx)
        Sz0(zero_idx) = 0;
    end

    % For omega~=0, s~=0. L{S(t)}(s) = (1 - P(s))/s
    if any(nonzero_idx)
        S_L_s = (1 - P_s(nonzero_idx)) ./ s_nonzero;
        Sz0(nonzero_idx) = -2*(1/mean_creation_rate).*imag(S_L_s)./omega(nonzero_idx);
    end

end
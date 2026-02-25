function [P_avg, omega] = calculate_psd(x, fs, window_name, window_param)
%CALCULATE_PSD Calculates the power spectral density of a signal or a set of signals.
%
%   [P_avg, omega] = CALCULATE_PSD(x, fs, window_name, window_param) calculates the
%   power spectral density (PSD) of a time-domain signal x using the Fast
%   Fourier Transform (FFT) and a specified window.
%
%   Inputs:
%   x            - The input signal. Can be a vector (single trajectory) or a
%                  matrix where each row is a trajectory.
%   fs           - The sampling frequency of the signal.
%   window_name  - (Optional) The name of the window function to use.
%                  Options: 'rectangular', 'hann', 'hamming', 'blackman', 'kaiser'.
%                  Default is 'rectangular'.
%   window_param - (Optional) Parameter for the window function. Currently only
%                  used for the 'kaiser' window (beta value). Default is 5.
%
%   Outputs:
%   P_avg - The single-sided power spectral density. If x is a matrix,
%           this is the average PSD over all trajectories.
%   omega- The corresponding frequency vector in rad/s.

    if nargin < 3
        window_name = 'rectangular';
    end
    if nargin < 4 && strcmpi(window_name, 'kaiser')
        window_param = 5; % Default beta for Kaiser window
    end

    if isvector(x)
        x = x(:)'; % Ensure x is a row vector
    end

    [num_trajectories, N] = size(x);

    % Create the window
    switch lower(window_name)
        case "rectangular"
            w = rectwin(N)';
        case "hann"
            w = hann(N)';
        case "hamming"
            w = hamming(N)';
        case "blackman"
            w = blackman(N)';
        case "kaiser"
            if nargin < 4
                window_param = 5; % Default beta for Kaiser window
            end
            w = kaiser(N, window_param)';
        otherwise
            error("Unknown window name: %s. Use 'rectangular', 'hann', 'hamming', 'blackman', or 'kaiser'.", window_name);
    end

    % Initialize matrix to store PSDs for each trajectory
    if num_trajectories > 1
        x_mean = mean(x, 1); % Ensemble average
    else
        x_mean = mean(x); % Time average for a single trajectory
    end
    x = x - x_mean;
    
    % Apply window to all trajectories
    x = x .* w; % w is broadcasted
    
    % Calculate FFT for all trajectories (row-wise)
    Y = fft(x, [], 2);
    
    % Calculate the power spectrum for all trajectories
    P2 = abs(Y).^2 / (fs * sum(w.^2)); % Power spectral density
    % P1 = abs(Y).^2 / (N*sum(w.^2));% Power spectrum, 

    % Extract single-sided spectrum
    P_matrix = P2(:, 1:floor(N/2)+1);

    % Average the PSDs
    P_avg = mean(P_matrix, 1);

    % Define the frequency domain (angular frequency)
    omega = 2*pi*fs*(0:(N/2))/N;
end

function phase_factor = factor(pulse_shape,ref_volt,pulse_volt,pulse_dur,pulse_freq)
% Compute the constant factor between the relative B1 map and the phase 
% shift introduced by the Bloch-Siegert pulse: Dphi = +/- K*b1^2
%
% FORMAT phase_factor = b1p.bs.factor(pulse_shape, pulse_freq, pulse_dur,
%                                     pulse_volt, ref_volt) 
% pulse_shape  - Pulse shape (gauss|fermi)
% pulse_freq   - Pulse off-resonance frequency [in Hz]
% pulse_dur    - Pulse duration [in s]
% pulse_volt   - Pulse voltage [in ?]
% ref_volt     - Reference voltage [in ?]
% phase_factor - Phase shift to squared relative B1 ratio [in rad]
%
% REFERENCES
% . Sacolick, L.I., Wiesinger, F., Hancu, I., Vogel, M.W., 2010. 
%   B1 mapping by Bloch-Siegert shift. Magnetic Resonance in Medicine 63, 
%   1315–1322. https://doi.org/10.1002/mrm.22357
% . Corbin, N., Acosta‐Cabronero, J., Malik, S.J., Callaghan, M.F., 2019. 
%   Robust 3D Bloch-Siegert based mapping using multi-echo general linear 
%   modeling. Magnetic Resonance in Medicine 82, 2003–2015. 
%   https://doi.org/10.1002/mrm.27851

% Constant parameters
gamma  = 2 * pi * 42.57e6;        % H1 resonance frequency (rad/s/T)
ref_b1 = pi / (gamma * 1e-3);     % Nominal B1 for 180 deg/ms

% Pulse centred at time t = 0
pulse_freq = 2 * pi * pulse_freq; % Convert pulse frequency (rad/s)
switch lower(pulse_shape)
    case 'gauss'
        pulse_integral = pulse_integral_gauss(pulse_dur);
    case 'fermi'
        pulse_integral = pulse_integral_fermi(pulse_dur);
end
phase_factor = 0.5 * gamma^2 *  pulse_integral / pulse_freq;

% Normalise using reference and BS voltage
phase_factor = phase_factor * ref_b1 * pulse_volt / ref_volt;

function pulse_integral = pulse_integral_gauss(pulse_dur)
% Compute the integral of the square of a Gaussian pulse.
%
%   The Gaussian pulse has shape exp(-t2/(2*s2)), where the
%   standard-deviation s depends on the pulse duration through s=d/7.4338.
%   It would be good to know where this magic number comes from.

% B1(t)^2 = exp(-t2/(2*s^2))^2 = exp(-2*t^2/(2*s^2)) = exp(-t^2/s^2)
% => int B1(t)2 = sqrt(pi) * erf(t/s)
pulse_width    = pulse_dur / 7.4338; % Convert duration to Gaussian s.d.
pulse_integral = sqrt(pi) * erf(pulse_dur/pulse_width);

function pulse_integral = pulse_integral_fermi(pulse_dur)
% Compute the integral of the square of a Fermi pulse.
%
%   The Fermi pulse shape is 1/(1 + exp(abs(t)-t0)/a), where a and t0 are
%   hard-coded pulse parameters. They should be probably be made parameters
%   of the function in the future.
corr_factor = pulse_dur*1e6/8192;  % based on Sequence code
t0 = 3e-3*corr_factor;             % pulse width, s, as per sequence
a  = 0.16e-3*corr_factor;          % timing parameter s.t. as a -> 0 Fermi -> Hard pulse
dt = 0.01e-3;
t  = -(pulse_dur/2):dt:(pulse_dur/2);
pulse_integral = 1 ./ (1 + exp( (abs(t) - t0) ./ a));
pulse_integral = dt * sum(pulse_integral.^2);
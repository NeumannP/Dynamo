function [out, desc] = export(self, in_polar, keep_TU)
% Export a control sequence from Dynamo into an array.
%
%  All control fields are divided by 2*pi, to convert angular
%  frequencies to normal frequencies.
%
%  If in_polar is true, converts (x,y) control fields to polar
%  coordinates (amplitude, phase). Assumes that all control fields
%  are of the (x,y) type.
%
%  NOTE that the amplitude corresponds to a Rabi frequency only if
%  the control Hamiltonians are properly normalized.
%
%  The output array rows correspond to pulses/bins, each row consists of 
%    [tau (s), amp_1 (Hz), amp_2 (Hz), ...]
%  or, if in_polar is true,
%    [tau (s), f_rabi_1 (Hz), \phi_1, f_rabi_2 (Hz), \phi_2, ...]

% Ville Bergholm 2013-2016


if nargin < 3
    keep_TU = false;
    if nargin < 2
        % default: export the controls as they are
        in_polar = false;
    end
end

TU = self.system.TU;
tu_str = '';
if isempty(TU)
    % not defined
    if keep_TU
        error('Time unit not set.');
    end
    TU = 1; % implicit default
    t_unit = '';
    f_unit = '';
else
    if keep_TU
        tu_str = sprintf('TU = %g s\n', TU);
        switch TU
          case 1e-3
            t_unit = 'ms';
            f_unit = 'kHz';
          case 1e-6
            t_unit = '\mu s';
            f_unit = 'MHz';
          case 1e-9
            t_unit = 'ns';
            f_unit = 'GHz';
          otherwise
            t_unit = 'TU';
            f_unit = '1/TU';
        end
    else
        t_unit = 's';
        f_unit = 'Hz';
    end
end

% from angular to normal frequency
f   = self.seq.fields / (2*pi);
tau = self.seq.tau;
if ~keep_TU
    f = f / TU;
    tau = tau * TU;
end

% polar transformation?
if in_polar
    %n_bins = size(f, 1);
    n_controls = size(f, 2) / 2;
    if ceil(n_controls) ~= n_controls
        error('Number of controls must be even for polar transformation.')
    end
    % use polar coordinates, assumes that all controls are x/y pairs
    legend = sprintf('[tau (%s), f_rabi_1 (%s), phi_1, f_rabi_2 (%s), phi_2, ...]', t_unit, f_unit, f_unit);

    for k=1:2:2*n_controls
        x = f(:, k);
        y = f(:, k+1);
        f(:, k) = sqrt(x.^2 +y.^2);
        f(:, k+1) = atan2(y, x);
    end
else
    % cartesian coordinates, or simply abstract unrelated controls
    legend = sprintf('[tau (%s), amp_1 (%s), amp_2 (%s), ...]', t_unit, f_unit, f_unit);
end

desc = sprintf(['%s\nGenerated: %s UTC\n'...
               'Rows correspond to pulses, each row consists of\n%s\n%s\n'],...
               self.system.description, self.config.date_UTC, legend, tu_str);
out = [tau, f];
end

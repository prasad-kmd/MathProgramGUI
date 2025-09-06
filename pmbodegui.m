function pmbodegui(sys)
    % PMBODEGUI - Creates a GUI to display Bode plots of a transfer function.
    %
    % Syntax: pmbodegui(sys)
    %
    % Description:
    %   pmbodegui(sys) creates a GUI that displays the Bode plot for the
    %   transfer function 'sys'. The GUI shows the individual term plots and
    %   the final asymptotic Bode plot.
    %
    %   'sys' is a transfer function object created using the tf() function.
    %

    % Create the main figure
    fig = figure('Name', 'Bode Plot GUI', 'Position', [100, 100, 1200, 600]);

    % --- UI Layout ---
    % Left Panel for Individual Terms
    left_panel = uipanel('Title', 'Individual Terms', 'FontSize', 12, ...
                         'Position', [0.05 0.2 0.43 0.75]);
    ax_mag_ind = axes('Parent', left_panel, 'Position', [0.1 0.55 0.85 0.35]);
    title(ax_mag_ind, 'Magnitude');
    ylabel(ax_mag_ind, 'Magnitude (dB)');
    grid(ax_mag_ind, 'on');
    grid(ax_mag_ind, 'minor');
    ax_phase_ind = axes('Parent', left_panel, 'Position', [0.1 0.1 0.85 0.35]);
    title(ax_phase_ind, 'Phase');
    ylabel(ax_phase_ind, 'Phase (deg)');
    xlabel(ax_phase_ind, 'Frequency (rad/s)');
    grid(ax_phase_ind, 'on');
    grid(ax_phase_ind, 'minor');

    % Right Panel for Asymptotic Bode Plot
    right_panel = uipanel('Title', 'Asymptotic Bode Plot', 'FontSize', 12, ...
                          'Position', [0.52 0.2 0.43 0.75]);
    ax_mag_final = axes('Parent', right_panel, 'Position', [0.1 0.55 0.85 0.35]);
    title(ax_mag_final, 'Magnitude');
    ylabel(ax_mag_final, 'Magnitude (dB)');
    grid(ax_mag_final, 'on');
    grid(ax_mag_final, 'minor');
    ax_phase_final = axes('Parent', right_panel, 'Position', [0.1 0.1 0.85 0.35]);
    title(ax_phase_final, 'Phase');
    ylabel(ax_phase_final, 'Phase (deg)');
    xlabel(ax_phase_final, 'Frequency (rad/s)');
    grid(ax_phase_final, 'on');
    grid(ax_phase_final, 'minor');

    % Info Text Area
    info_panel = uipanel('Title', 'Transfer Function Info', 'FontSize', 12, ...
                         'Position', [0.05 0.02 0.9 0.15]);
    info_text = uicontrol('Parent', info_panel, 'Style', 'text', 'FontSize', 10, ...
                          'HorizontalAlignment', 'left', 'Units', 'normalized', ...
                          'Position', [0.01 0.05 0.98 0.9]);

    % --- Transfer function parsing ---
    [num, den] = tfdata(sys, 'v');
    [z, p, k] = zpkdata(sys, 'v');

    terms = {};

    % Gain
    if k ~= 1
        terms{end+1} = sprintf('Gain: %g', k);
    end

    % Zeros
    for i = 1:length(z)
        if z(i) == 0
            terms{end+1} = 'Differentiator: s';
        elseif isreal(z(i))
            terms{end+1} = sprintf('1st Order Zero: (1 + s/%g)', -z(i));
        else % Complex
            if imag(z(i)) > 0 % Process only one of the conjugate pair
                wn = abs(z(i));
                zeta = -real(z(i))/wn;
                terms{end+1} = sprintf('2nd Order Zero: (1 + 2*%g*s/%g + (s/%g)^2)', zeta, wn, wn);
            end
        end
    end

    % Poles
    for i = 1:length(p)
        if p(i) == 0
            terms{end+1} = 'Integrator: 1/s';
        elseif isreal(p(i))
            terms{end+1} = sprintf('1st Order Pole: 1/(1 + s/%g)', -p(i));
        else % Complex
            if imag(p(i)) > 0 % Process only one of the conjugate pair
                wn = abs(p(i));
                zeta = -real(p(i))/wn;
                terms{end+1} = sprintf('2nd Order Pole: 1/(1 + 2*%g*s/%g + (s/%g)^2)', zeta, wn, wn);
            end
        end
    end

    % --- Plotting ---

    % Define frequency range
    all_freqs = abs([p.', z.']);
    all_freqs = all_freqs(all_freqs > 0);
    if isempty(all_freqs)
        w = logspace(-2, 2, 400);
    else
        w = logspace(floor(log10(min(all_freqs)))-1, ceil(log10(max(all_freqs)))+1, 400);
    end

    hold(ax_mag_ind, 'on');
    hold(ax_phase_ind, 'on');

    legend_entries = {};
    mag_final = zeros(size(w));
    phase_final = zeros(size(w));

    % Plot individual terms
    % Gain
    if k ~= 0
        mag = 20*log10(abs(k)) * ones(size(w));
        phase = (k<0)*-180 * ones(size(w));
        semilogx(ax_mag_ind, w, mag, 'LineWidth', 1.5);
        semilogx(ax_phase_ind, w, phase, 'LineWidth', 1.5);
        legend_entries{end+1} = 'Gain';
        mag_final = mag_final + mag;
        phase_final = phase_final + phase;
    end

    % Zeros
    for i = 1:length(z)
        if z(i) == 0 % Differentiator
            mag = 20*log10(w);
            phase = 90 * ones(size(w));
            legend_entries{end+1} = 'Differentiator';
        elseif isreal(z(i)) % 1st order zero
            w_z = -z(i);
            mag = zeros(size(w));
            mag(w >= w_z) = 20*log10(w(w >= w_z)/w_z);
            phase = zeros(size(w));
            phase(w > w_z*10) = 90;
            sel = w >= w_z/10 & w <= w_z*10;
            phase(sel) = 45 * (log10(w(sel)) - log10(w_z/10));
            legend_entries{end+1} = sprintf('Zero at %g', w_z);
        else % 2nd order zero
            if imag(z(i)) > 0
                wn = abs(z(i));
                mag = zeros(size(w));
                mag(w >= wn) = 40*log10(w(w >= wn)/wn);
                phase = zeros(size(w));
                phase(w > wn*10) = 180;
                sel = w >= wn/10 & w <= wn*10;
                phase(sel) = 90 * (log10(w(sel)) - log10(wn/10));
                legend_entries{end+1} = sprintf('Complex Zero at %g rad/s', wn);
            else
                continue; % Skip conjugate pair
            end
        end
        semilogx(ax_mag_ind, w, mag, '--', 'LineWidth', 1.5);
        semilogx(ax_phase_ind, w, phase, '--', 'LineWidth', 1.5);
        mag_final = mag_final + mag;
        phase_final = phase_final + phase;
    end

    % Poles
    for i = 1:length(p)
        if p(i) == 0 % Integrator
            mag = -20*log10(w);
            phase = -90 * ones(size(w));
            legend_entries{end+1} = 'Integrator';
        elseif isreal(p(i)) % 1st order pole
            w_p = -p(i);
            mag = zeros(size(w));
            mag(w >= w_p) = -20*log10(w(w >= w_p)/w_p);
            phase = zeros(size(w));
            phase(w > w_p*10) = -90;
            sel = w >= w_p/10 & w <= w_p*10;
            phase(sel) = -45 * (log10(w(sel)) - log10(w_p/10));
            legend_entries{end+1} = sprintf('Pole at %g', w_p);
        else % 2nd order pole
             if imag(p(i)) > 0
                wn = abs(p(i));
                mag = zeros(size(w));
                mag(w >= wn) = -40*log10(w(w >= wn)/wn);
                phase = zeros(size(w));
                phase(w > wn*10) = -180;
                sel = w >= wn/10 & w <= wn*10;
                phase(sel) = -90 * (log10(w(sel)) - log10(wn/10));
                legend_entries{end+1} = sprintf('Complex Pole at %g rad/s', wn);
             else
                 continue; % Skip conjugate
             end
        end
        semilogx(ax_mag_ind, w, mag, ':', 'LineWidth', 1.5);
        semilogx(ax_phase_ind, w, phase, ':', 'LineWidth', 1.5);
        mag_final = mag_final + mag;
        phase_final = phase_final + phase;
    end

    legend(ax_mag_ind, legend_entries);
    legend(ax_phase_ind, legend_entries);

    % Plot final asymptotic Bode plot
    semilogx(ax_mag_final, w, mag_final, 'LineWidth', 2, 'Color', 'r');
    semilogx(ax_phase_final, w, phase_final, 'LineWidth', 2, 'Color', 'r');

    hold(ax_mag_ind, 'off');
    hold(ax_phase_ind, 'off');

    % Display info
    tf_str_raw = evalc('disp(sys)');
    tf_lines = strsplit(strtrim(tf_str_raw), '\n');
    terms_str = ['Terms: ' strjoin(terms, ', ')];
    display_str = [tf_lines, {terms_str}];
    set(info_text, 'String', display_str);

end

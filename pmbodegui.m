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
                         'Position', [0.05 0.1 0.4 0.85]);
    ax_mag_ind = axes('Parent', left_panel, 'Position', [0.1 0.55 0.8 0.35]);
    title(ax_mag_ind, 'Magnitude');
    ylabel(ax_mag_ind, 'Magnitude (dB)');
    grid(ax_mag_ind, 'on');
    ax_phase_ind = axes('Parent', left_panel, 'Position', [0.1 0.1 0.8 0.35]);
    title(ax_phase_ind, 'Phase');
    ylabel(ax_phase_ind, 'Phase (deg)');
    xlabel(ax_phase_ind, 'Frequency (rad/s)');
    grid(ax_phase_ind, 'on');

    % Right Panel for Asymptotic Bode Plot
    right_panel = uipanel('Title', 'Asymptotic Bode Plot', 'FontSize', 12, ...
                          'Position', [0.55 0.1 0.4 0.85]);
    ax_mag_final = axes('Parent', right_panel, 'Position', [0.1 0.55 0.8 0.35]);
    title(ax_mag_final, 'Magnitude');
    ylabel(ax_mag_final, 'Magnitude (dB)');
    grid(ax_mag_final, 'on');
    ax_phase_final = axes('Parent', right_panel, 'Position', [0.1 0.1 0.8 0.35]);
    title(ax_phase_final, 'Phase');
    ylabel(ax_phase_final, 'Phase (deg)');
    xlabel(ax_phase_final, 'Frequency (rad/s)');
    grid(ax_phase_final, 'on');

    % Info Text Area
    info_panel = uipanel('Title', 'Transfer Function Info', 'FontSize', 12, ...
                         'Position', [0.05 0.01 0.9 0.08]);
    info_text = uicontrol('Parent', info_panel, 'Style', 'text', 'FontSize', 10, ...
                          'HorizontalAlignment', 'left', ...
                          'Position', [10 5 800 20]);

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
            mag = 20*log10(sqrt(1 + (w/w_z).^2));
            phase = atand(w/w_z);
            legend_entries{end+1} = sprintf('Zero at %g', w_z);
        else % 2nd order zero
            if imag(z(i)) > 0
                wn = abs(z(i));
                zeta = -real(z(i))/wn;
                u = w/wn;
                mag = 20*log10(sqrt((1-u.^2).^2 + (2*zeta*u).^2));
                phase = atand((2*zeta*u)./(1-u.^2));
                phase(u > 1) = phase(u > 1) + 180;
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
            mag = -20*log10(sqrt(1 + (w/w_p).^2));
            phase = -atand(w/w_p);
            legend_entries{end+1} = sprintf('Pole at %g', w_p);
        else % 2nd order pole
             if imag(p(i)) > 0
                wn = abs(p(i));
                zeta = -real(p(i))/wn;
                u = w/wn;
                mag = -20*log10(sqrt((1-u.^2).^2 + (2*zeta*u).^2));
                phase = -atand((2*zeta*u)./(1-u.^2));
                phase(u > 1) = phase(u > 1) - 180;
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
    tf_str = evalc('disp(sys)');
    info_str = sprintf('G(s) = %s\nTerms: %s', tf_str, strjoin(terms, ', '));
    set(info_text, 'String', info_str);

end

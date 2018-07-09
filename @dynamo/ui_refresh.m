function ui_refresh(self, full_plot, err)
% Refreshes the UI figure, plots stuff.
    
    % Do a full refresh or just redraw the graphical objects?
    if nargin < 2
        full_plot = true;
    end
    
    h = self.config.ui_fig;
    % It's incredible how much work it takes just to make
    % MATLAB not steal window focus when it plots something.
    set(0, 'CurrentFigure', h);

    if 0
        ax = subplot(2, 1, 1);
        self.plot_seq(ax, full_plot);

        ax = subplot(2, 1, 2);
        self.plot_X([], [], ax, full_plot);
        %self.plot_stats(ax);
    else
        ax = get(h, 'CurrentAxes');
        self.plot_seq(ax, full_plot);
        title(ax, self.system.description);
    end
    if nargin == 3
        if self.config.frobenius_like_error
            % plot consistent Frobenius norm error
            err = sqrt(2 * self.system.norm2 * err);
            str = 'Frobenius norm';
        else
            % currently overlap errors are the only ones that are
            % incommensurate with Frobenius norm errors
            str = 'Overlap';
        end
        text(0.05, 0.9, sprintf('%s error: %6.6g', str, err), 'Units','normalized',...
             'FontSize',18, 'BackgroundColor',[.8 .8 1])
    end
    drawnow();
end

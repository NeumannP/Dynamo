function ui_open(self)
% Create the UI figure.

    self.config.ui_fig = figure('Name', 'Dynamo UI', 'CloseRequestFcn', {@close_req, self});
    ax = axes();
    self.ui_refresh(true); % initial draw
end


function close_req(src, event, self)
% Callback for closing the UI figure window.
    delete(src); % close the figure
    self.config.ui_fig = [];
    self.config.stop = true; % signal monitor_func that we should stop
end

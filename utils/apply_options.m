function [out, unused] = apply_options(defaults, opts, only_existing)
% MATLAB-style options struct processing.
% Applies the options in the struct 'opts' to the struct 'defaults'.
% Returns the updated struct, and the struct of options that could not be parsed.

% Ville Bergholm 2015-2016

    % fields in opts
    names = fieldnames(opts);

    if only_existing
        % logical array: are the corresponding fields present in defaults?
        present = isfield(defaults, names);
    else
        present = true(size(names));
    end

    % could not find a function for doing this
    out = defaults;
    for f = names(present).'
        out = setfield(out, f{1}, getfield(opts, f{1}));
    end
    % remove the fields we just used from opts
    unused = rmfield(opts, names(present));
end

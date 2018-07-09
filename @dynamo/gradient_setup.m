function gradient_setup(self, control_mask)
% Set up all the required _needed_now fields for computing the gradient.
% TODO kind of a hack, the gradient funcs themselves should maybe keep this information


% Find out which Hs, Us, & Ls we need to compute the gradient.
% It's usually more efficient to do all the required calculations at once, and not piece-meal.
slot_mask = any(control_mask, 2);
self.cache.H_needed_now(slot_mask) = true;           % H_{slot}
self.cache.L_needed_now([false; slot_mask]) = true;  % L_{slot+1}

switch self.config.dP
  case 'series'
    self.cache.U_needed_now([false; slot_mask]) = true;  % U_{slot+1}

  otherwise
    % taus and other controls require different things
    tau_slot_mask = control_mask(:, end);
    c_slot_mask   = any(control_mask(:, 1:end-1), 2);
    u_slot_mask = [c_slot_mask; false] | [false; tau_slot_mask];
    
    % P_{c_slot} is required by eig, series_ss and fd_dp
    % TODO not by aux, fd?
    self.cache.P_needed_now(c_slot_mask) = true;  % P_{c_slot}, for H_v and H_eig_factor
    self.cache.U_needed_now(u_slot_mask) = true;  % U_{c_slot},  U_{tau_slot+1}
end

if self.config.UL_mixed
    self.cache.U_needed_now([false; slot_mask]) = true;
end

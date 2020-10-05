function LAST_scheduler
%



% get TargetList including history


% run night simulation in order to estimate the footprints of fast cadence

% prep TargetList for hig cadence

% prep TargetList for low cadence

% prep Dynamic TargetList (ToO) - set to empty

% start with most limited telescopes (1,2)
% for each telescope
    % For each telescope mount call target_selection

    % telescopes are pre-assigned for low/high cadence  (????)
    % 1 3 5 7 9  11
    % 2 4 6 8 10 12
    % by default 5,6,9,10 are assigned to fast cadence

    % For each telescope mount call validate_coo.m 
% end
% output is list of possible targets, target category (fast/slow), and
% which telescope can observe the target


% Select target for telescope
% and mark the selected targets
% AD9361 Filter Wizard requires MATLAB 2015b or later, older versions are unsupported.

% Returns 1 if running a supported version, 0 otherwise.
function result = check_compat()
result = ~verLessThan('matlab', '8.6')

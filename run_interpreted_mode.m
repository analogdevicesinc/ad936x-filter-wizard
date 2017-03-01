%% Generate example coefficients in interpreted mode

% Get example config
a = load('ad9361_settings.mat');
in = a.ad9361_settings.tx.LTE5;
% Fill out all fields
input = process_input(in);
% Get filter coefficients
callMexVersion = false;
taps = call_filter_designer_cg(input,callMexVersion); % Taps are int16
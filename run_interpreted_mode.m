% Generate coefficients in interpreted mode

% Get example config
a = load('ad9361_settings.mat');
in = a.ad9361_settings.tx.LTE5;
% Fill out all fields
input = process_input(in);
% Manipulate input struct
input.int_FIR = 0;
% Get filter coefficients
out = filter_designer_cg(input,false);
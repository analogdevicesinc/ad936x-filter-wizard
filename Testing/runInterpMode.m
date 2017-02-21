

a = load('ad9361_settings.mat');

in = a.ad9361_settings.tx.LTE5;

input = process_input(in);
input.int_FIR = 0;
out = filter_designer_cg(input,false);
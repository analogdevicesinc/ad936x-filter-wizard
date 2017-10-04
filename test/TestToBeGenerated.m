

function [output,r] = TestToBeGenerated()

% Get example inputs and reference outputs
a = load('ad9361_settings_processed_test.mat');
input = a.input;
firtaps = a.firtaps;

% Return preallocation and type def
output = zeros(1,128,'int16');
numOutputTaps = 1;

% Initialize generated designer
coder.ceval('internal_design_filter_cg_initialize');

% Call generated version
coder.ceval('internal_design_filter_cg',...
    input.Rdata,input.Fpass, input.Fstop, input.caldiv, input.FIR,...
    input.HB1, input.PLL_mult, input.Apass, input.Astop, input.phEQ,...
    input.HB2, input.HB3, input.Type,input.RxTx, input.RFbw, ...
    input.DAC_div, input.converter_rate, input.PLL_rate, input.Fcenter,...
    input.wnom, input.FIRdBmin, input.int_FIR, coder.wref(output), ...
    coder.wref(numOutputTaps));

% Teardown generated designer
coder.ceval('internal_design_filter_cg_terminate');

% Check outputs
r = isequal(output(1:numOutputTaps),int16(firtaps));


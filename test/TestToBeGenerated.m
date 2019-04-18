

function [output,r,Apass,Astop,maxInputFS,maxInputDB] = TestToBeGenerated()

% Get example inputs and reference outputs
a = load('ad9361_settings_processed_test.mat');
input = a.input;
firtaps = a.firtaps;

if input.HB3 == 2
    hb3 = 2;
elseif input.HB3 == 3
    hb3 = 1;
else
    hb3 = 1;
end
if strcmp(input.RxTx, 'Rx')
    if hb3 == 1
        N = min(16*floor(input.converter_rate/(input.Rdata)),128);
    else
        N = min(16*floor(input.converter_rate/(2*input.Rdata)),128);
    end
else
    switch input.FIR
        case 1
            Nmax = 64;
        case 2
            Nmax = 128;
        case 4
            Nmax = 128;
        otherwise
            error('Wrong FIR Type');
    end
    N = min(16*floor(input.converter_rate*input.DAC_div/(2*input.Rdata)),Nmax);
end

% Return preallocation and type def
output = zeros(1,128,'int16');
numOutputTaps = 1;
filterGain = 0;
Apass = 0;
Astop = 0;
maxInputFS = 0;
maxInputDB = 0;

% Initialize generated designer
coder.ceval('internal_design_filter_cg_initialize');

% Call generated version
coder.ceval('internal_design_filter_cg',...
    input.Rdata,input.Fpass, input.Fstop, input.caldiv, input.FIR,...
    input.HB1, input.PLL_mult, input.Apass, input.Astop, input.phEQ,...
    input.HB2, input.HB3, input.Type,input.RxTx, input.RFbw, ...
    input.DAC_div, input.converter_rate, input.PLL_rate, input.Fcenter,...
    input.wnom, input.FIRdBmin, input.int_FIR, N, coder.wref(output), ...
    coder.wref(numOutputTaps),coder.wref(filterGain),...
    coder.wref(Apass),coder.wref(Astop),...
    coder.wref(maxInputFS),coder.wref(maxInputDB));

% Teardown generated designer
coder.ceval('internal_design_filter_cg_terminate');

% Check outputs
r = isequal(output(1:numOutputTaps),int16(firtaps));


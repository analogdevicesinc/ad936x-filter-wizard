% Example generation of code for internal_filter_designer_cg.m

addpath(genpath('Testing')); % Grab helper functions to manage structs
a = load('ad9361_settings.mat');
inputVar = a.ad9361_settings.tx.LTE5;
% Fill out necessary fields
input = process_input(inputVar);
args = '{input.Rdata, input.Fpass, input.Fstop, input.caldiv, input.FIR, input.HB1, input.PLL_mult, input.Apass, input.Astop, input.phEQ, input.HB2, input.HB3, input.Type, input.RxTx, input.RFbw, input.DAC_div, input.converter_rate, input.PLL_rate, input.Fcenter, input.wnom, input.FIRdBmin, input.int_FIR}';
functionName = 'internal_design_filter_cg';

%% Build library
% This will generate code in the codegen folder and compile a dll/so/dylib
% with the necessary headers
cfg = coder.config('dll');
cfg.TargetLang = 'C++';
cfg.FilePartitionMethod = 'SingleFile';
outputLIBName = 'libinternal_filter_designer';
result = codegen('-config','cfg',functionName,'-O ','disable:openmp','-args', args,'-o',outputLIBName);

%% Build executable
% This will generate code in the codegen folder and compile an exe with the
% ex_main.cpp file provided in the 'cpp' folder
cfg = coder.config('exe');
cfg.TargetLang = 'C++';
cfg.FilePartitionMethod = 'SingleFile';
cfg.CustomSource = 'ex_main.cpp';
cfg.CustomInclude = 'cpp';
outputEXEName = 'designer';
result = codegen('-config','cfg',functionName,'-O ','disable:openmp','-args', args,'-o',outputEXEName);

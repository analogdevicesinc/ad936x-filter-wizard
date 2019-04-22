function output = internal_design_filter_opt_ripple(input)

% support a simple data rate input otherwise it must be a structure
if isfloat(input)
    input = struct('Rdata', input);
else
    error('Input must be the required samplerate only');
end

input = cook_input(input);

% use the internal FIR if unspecified
if ~isfield(input, 'int_FIR')
    input.int_FIR = 1;
end

% nominal frequency can't be zero
if ~input.wnom
    input.wnom = calculate_rfbw(input.PLL_rate, input.caldiv, input.RxTx, true);
    input.wnom = double(input.wnom);
end

FIR = [1,2,4];
outs = [];
ripples = [];
for i = 1:3
    % Recalc clocks
    [input0, configIndx] = auto_fast_rates(input,FIR(i));
    if configIndx == 0
        continue;
    end
    % Update position of wnom since PLL can change
    if strcmp(input.RxTx, 'Rx')
        input0.wnom = 1.4 * input0.Fstop; % Rx
    else
        input0.wnom = 1.6 * input0.Fstop; % Tx
    end
    div = ceil((input0.PLL_rate/input0.wnom)*(log(2)/(2*pi)));
    input0.caldiv = min(max(div,1),511);
    input0.wnom = double(calculate_rfbw(input0.PLL_rate, input0.caldiv,...
        input0.RxTx, true));
    % Design filter
    out = internal_design_filter(input0);
    outs = [outs; out]; %#ok<AGROW>
    ripples = [ripples;out.Apass_actual]; %#ok<AGROW>
end

if isempty(ripples)
    error('No valid design found');
end

[~,i] = min(ripples);
output = outs(i);


end
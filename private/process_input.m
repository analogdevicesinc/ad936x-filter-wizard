function input = process_input(input)

% support a simple data rate input otherwise it must be a structure
if isfloat(input)
    input = struct('Rdata', input);
end
input = cook_input(input);

% use the internal FIR if unspecified
if ~isfield(input, 'int_FIR')
    input.int_FIR = 1;
end

% nominal frequency can't be zero
if ~input.wnom
    input.wnom = (input.PLL_rate/input.caldiv)*(log(2)/(2*pi));
end

end
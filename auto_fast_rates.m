function [input,savedFilterConfig] = auto_fast_rates(input)
%auto_fast_rates Configures halfband filters, PLL clock, and
%converter rates in an optimal configuration based on the
%current 'Rdata' property

bounds9361 = rate_bounds;

if input.Rdata>=bounds9361.MAX_DATA_RATE
    input.Rdata = (bounds9361.MAX_DATA_RATE);
elseif input.Rdata<=bounds9361.MIN_DATA_RATE
    input.Rdata = (bounds9361.MIN_DATA_RATE);
end

currentMaxADCRate = 0;
savedDACDivider = 0;
rate_gov = false;
savedFIR = 0;
savedPLL = 0;
savedFilterConfig = 0;

DR = input.Rdata;

for FIR = [4,2,1]
    
    FilterConfig = [...
        3,2,2,FIR;...
        2,2,2,FIR;...
        3,2,1,FIR;...
        2,2,1,FIR;...
        2,1,1,FIR;...
        3,1,1,FIR;...
        1,1,1,FIR];
    
    for k=1:7
        % Update configuration
        divs = struct;
        divs.RxFIR = FilterConfig(k,4);
        divs.RxHB1 = FilterConfig(k,3);
        divs.RxHB2 = FilterConfig(k,2);
        divs.RxHB3 = FilterConfig(k,1);
        divs.TxFIR = FilterConfig(k,4);
        divs.TxHB1 = FilterConfig(k,3);
        divs.TxHB2 = FilterConfig(k,2);
        divs.TxHB3 = FilterConfig(k,1);
        
        rates = gen_rates(divs,DR);
        
        % HB3 cannot be 3 if rate_gov enabled
        if rate_gov && (divs.RxHB3==3)
            continue;
        end
        
        % Check valid rates (ignore some section pre-dacdiv)
        valid = validatePathRates(rates);
        v = valid.RxRateHB3.pass && valid.RxRateHB2.pass && valid.RxRateHB1.pass && valid.RxRateFIR.pass;
        v = v && valid.TxRateHB2.pass && valid.TxRateHB1.pass && valid.TxRateFIR.pass;
        if v
            % Determine PLL divider
            pll = determine_pll_div(rates.ADCRate);
            if pll>0
                rates.PLLRate = rates.ADCRate * pll;
                % Determine DAC divider setting and check ADC/DAC settings
                dac_div = check_dac_adc_config(pll,k,rates.PLLRate);
                if dac_div>0
                    rates.DACRate = rates.ADCRate/dac_div;
                    if rates.ADCRate>currentMaxADCRate
                        % Final check
                        [~,v] = validatePathRates(rates);
                        if ~v
                            continue
                        end
                        currentMaxADCRate = rates.ADCRate;
                        savedFilterConfig = k;
                        savedDACDivider = dac_div;
                        savedFIR = FIR;
                        savedPLL = pll;
                    end
                end
            end
        end
    end
end


% Set HBs based on best config found
PLLDivider = savedPLL;
RxFIR = savedFIR;
RxHB1 = FilterConfig(savedFilterConfig,3);
RxHB2 = FilterConfig(savedFilterConfig,2);
RxHB3 = FilterConfig(savedFilterConfig,1);
TxFIR = savedFIR;
TxHB1 = FilterConfig(savedFilterConfig,3);
TxHB2 = FilterConfig(savedFilterConfig,2);
TxHB3 = FilterConfig(savedFilterConfig,1);
% Apply dac div
if (savedFilterConfig<3) && (savedDACDivider>1)
    TxHB1 = 1;
elseif (savedFilterConfig<5) && (savedDACDivider>1)
    TxHB2 = 1;
end

% Save
input.DAC_div = savedDACDivider;
input.PLL_mult = PLLDivider;
if strcmp(input.RxTx, 'Rx')
    input.FIR = RxFIR;
    input.HB1 = RxHB1;
    input.HB2 = RxHB2;
    input.HB3 = RxHB3;
    input.converter_rate = currentMaxADCRate;
else
    input.FIR = TxFIR;
    input.HB1 = TxHB1;
    input.HB2 = TxHB2;
    input.HB3 = TxHB3;
    input.converter_rate = currentMaxADCRate/savedDACDivider;
end
input.PLL_rate = input.converter_rate * input.DAC_div * input.PLL_mult;

% Final check
if strcmp(input.RxTx, 'Rx')
    if input.converter_rate > bounds9361.MAX_ADC_CLK || input.converter_rate < bounds9361.MIN_ADC_CLK
        error('Invalid ADC Rate');
    end
else
    if input.converter_rate > bounds9361.MAX_DAC_CLK || input.converter_rate < bounds9361.MIN_DAC_CLK
        error('Invalid DAC Rate');
    end
end


end


function pll = determine_pll_div(ADCRate)

% Determine necessary PLL multiplier
PLL_mult = 64; %MAX_BBPLL_DIV

bounds9361 = rate_bounds;

while (PLL_mult>1)
    rate0 = ADCRate*PLL_mult;
    v = (rate0 >= bounds9361.MIN_BBPLL_FREQ) && (rate0 <= bounds9361.MAX_BBPLL_FREQ);
    if v
        pll = PLL_mult;
        return
    end
    PLL_mult = PLL_mult/2;
end
pll = -1;

end

function r = check_dac_adc_config(PLL_mult,dec_table_index,PLLRate)

with_dd = PLLRate/PLL_mult/2;
without_dd = PLLRate/PLL_mult/1;

bounds9361 = rate_bounds;

a = rangeCheck(with_dd,'DAC Rate',bounds9361.MIN_DAC_CLK,bounds9361.MAX_DAC_CLK);
b = rangeCheck(without_dd,'ADC Rate',bounds9361.MIN_ADC_CLK,bounds9361.MAX_ADC_CLK );
c = rangeCheck(without_dd,'DAC Rate',bounds9361.MIN_DAC_CLK,bounds9361.MAX_DAC_CLK);

if (c.pass && b.pass)
    r = 1; % Run without dac div
elseif (a.pass && b.pass && (dec_table_index<6))
    r = 2; % Run with dac div
else
    r = -1; % All rates invalid
end

end

function e = rangeCheck(val,name,min,max)

if ~isreal(val) || isnan(val) || ~isscalar(val)
    e = struct;
    e.msg = sprintf('%s must be integer, scalar, and finite',name);
    e.pass = false;
elseif val<min || val>max
    e = struct;
    e.msg = sprintf('%s must be in range [%d %d]',name,int32(min),int32(max));
    e.pass = false;
else
    e = struct;
    e.msg = '';
    e.pass = true;
end

end

% Path Rate Validations
function [valid,v] = validatePathRates(rates)

bounds9361 = rate_bounds;

valid = struct;
valid.ADCRate = rangeCheck(rates.ADCRate,'ADC Rate',bounds9361.MIN_ADC_CLK,bounds9361.MAX_ADC_CLK);
% Rates into stage
valid.RxRateHB3 = rangeCheck(rates.RxRateHB3,'Rx HB3',0,bounds9361.MAX_RX.HB3);
valid.RxRateHB2 = rangeCheck(rates.RxRateHB2,'Rx HB2',0,bounds9361.MAX_RX.HB2);
valid.RxRateHB1 = rangeCheck(rates.RxRateHB1,'Rx HB1',0,bounds9361.MAX_RX.HB1);
valid.RxRateFIR = rangeCheck(rates.RxRateFIR,'Rx FIR',0,bounds9361.MAX_FIR);
valid.DACRate = rangeCheck(rates.DACRate,'DAC Rate',bounds9361.MIN_DAC_CLK,bounds9361.MAX_DAC_CLK);
% Rates out of stage
valid.TxRateHB3 = rangeCheck(rates.TxRateHB3,'Tx HB3',0,bounds9361.MAX_TX.HB3);
valid.TxRateHB2 = rangeCheck(rates.TxRateHB2,'Tx HB2',0,bounds9361.MAX_TX.HB2);
valid.TxRateHB1 = rangeCheck(rates.TxRateHB1,'Tx HB1',0,bounds9361.MAX_TX.HB1);
valid.TxRateFIR = rangeCheck(rates.TxRateFIR,'Tx FIR',0,bounds9361.MAX_FIR);
valid.PLLRate = rangeCheck(rates.PLLRate,'PLL Rate',bounds9361.MIN_BBPLL_FREQ,bounds9361.MAX_BBPLL_FREQ);

v = valid.ADCRate.pass && valid.RxRateHB3.pass && valid.RxRateHB2.pass && valid.RxRateHB1.pass && valid.RxRateFIR.pass;
v = v && valid.DACRate.pass && valid.TxRateHB3.pass && valid.TxRateHB2.pass && valid.TxRateHB1.pass && valid.TxRateFIR.pass;
v = v && valid.PLLRate.pass;

end

function rates = gen_rates(divs,DR)

rates.DataRate = DR;

rates.RxRateFIR = rates.DataRate * divs.RxFIR;
rates.RxRateHB1 = rates.RxRateFIR * divs.RxHB1;
rates.RxRateHB2 = rates.RxRateHB1 * divs.RxHB2;
rates.RxRateHB3 = rates.RxRateHB2 * divs.RxHB3;

rates.TxRateFIR = rates.DataRate * divs.TxFIR;
rates.TxRateHB1 = rates.TxRateFIR * divs.TxHB1;
rates.TxRateHB2 = rates.TxRateHB1 * divs.TxHB2;
rates.TxRateHB3 = rates.TxRateHB2 * divs.TxHB3;

rates.ADCRate = rates.RxRateHB3;

rates.PLLRate = -1;
rates.DACRate = -1;

end

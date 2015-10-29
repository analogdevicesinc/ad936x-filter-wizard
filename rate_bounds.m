% AD9361/AD9364 specific max/min clock rates
bounds.MAX_BBPLL_FREQ = 1430000000;                         % 1430.0 MHz
bounds.MIN_BBPLL_FREQ =  715000000;                         %  715.0 MHz
bounds.MAX_ADC_CLK    =  640000000;                         %  640.0 MHz
bounds.MIN_ADC_CLK    =  bounds.MIN_BBPLL_FREQ / (2 ^ 6);  %   11.2 MHz
bounds.MAX_DAC_CLK    =  bounds.MAX_ADC_CLK / 2;           % (MAX_ADC_CLK / 2)
bounds.MAX_DATA_RATE  =   61440000;                         %   61.44 MSPS
bounds.MIN_DATA_RATE  =  bounds.MIN_BBPLL_FREQ / (48 * (2 ^ 6));
bounds.MAX_FIR        =  bounds.MAX_DATA_RATE * 2;
bounds.MAX_RX.HB1     =  245760000;
bounds.MAX_RX.HB2     =  320000000;
bounds.MAX_RX.HB3     =  640000000;
bounds.MAX_TX.HB1     =  160000000;
bounds.MAX_TX.HB2     =  320000000;
bounds.MAX_TX.HB3     =  320000000;

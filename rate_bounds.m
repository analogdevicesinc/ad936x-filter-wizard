function bounds9361 = rate_bounds
% AD9361/AD9364 specific max/min clock rates
bounds9361.MAX_BBPLL_FREQ = 1430000000;                         % 1430.0 MHz
bounds9361.MIN_BBPLL_FREQ =  715000000;                         %  715.0 MHz
bounds9361.MAX_ADC_CLK    =  640000000;                         %  640.0 MHz
bounds9361.MIN_ADC_CLK    =   25000000;                         %   25.0 MHz
bounds9361.MAX_DAC_CLK    =  bounds9361.MAX_ADC_CLK / 2;
bounds9361.MIN_DAC_CLK    =   25000000;                         %   25.0 MHz
bounds9361.MAX_DATA_RATE  =   61440000;                         %   61.44 MSPS
bounds9361.MIN_DATA_RATE  =  bounds9361.MIN_ADC_CLK / 48;           %  520.83 kSPS
bounds9361.MAX_FIR        =  bounds9361.MAX_DATA_RATE * 2;
bounds9361.MAX_RX.HB1     =  245760000;
bounds9361.MAX_RX.HB2     =  320000000;
bounds9361.MAX_RX.HB3     =  640000000;
bounds9361.MAX_TX.HB1     =  160000000;
bounds9361.MAX_TX.HB2     =  320000000;
bounds9361.MAX_TX.HB3     =  320000000;

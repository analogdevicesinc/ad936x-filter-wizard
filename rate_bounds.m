% AD9361/AD9364 specific max/min clock rates
boundsAD.MAX_BBPLL_FREQ = 1430000000;                         % 1430.0 MHz
boundsAD.MIN_BBPLL_FREQ =  715000000;                         %  715.0 MHz
boundsAD.MAX_ADC_CLK    =  640000000;                         %  640.0 MHz
boundsAD.MIN_ADC_CLK    =  boundsAD.MIN_BBPLL_FREQ / (2 ^ 6);  %   11.2 MHz
boundsAD.MAX_DAC_CLK    =  boundsAD.MAX_ADC_CLK / 2;           % (MAX_ADC_CLK / 2)
boundsAD.MAX_DATA_RATE  =   61440000;                         %   61.44 MSPS
boundsAD.MIN_DATA_RATE  =  boundsAD.MIN_BBPLL_FREQ / (48 * (2 ^ 6));
boundsAD.MAX_FIR        =  boundsAD.MAX_DATA_RATE * 2;
boundsAD.MAX_RX.HB1     =  245760000;
boundsAD.MAX_RX.HB2     =  320000000;
boundsAD.MAX_RX.HB3     =  640000000;
boundsAD.MAX_TX.HB1     =  160000000;
boundsAD.MAX_TX.HB2     =  320000000;
boundsAD.MAX_TX.HB3     =  320000000;

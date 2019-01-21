classdef ad936x < handle
    % ad936x is reponsible for managing ad936x transceiver settings and
    % their validation
    %
    %   ad936x methods:
    %
    %   AutoSetRates   - Configure halfbands and other clocks automatically
    %
    properties
        %ChipType Chip Type
        %   Target AD936X variant for filter. Possible options are:
        %   {'AD9361','AD9363','AD9364','Pluto'}
        ChipType = 'AD9361';
        %DataRate Data Rate
        %   Complex data rate in Samples per second at output of FIR on
        %   receive path and input to FIR on transmit path
        DataRate = 7680000;
        %TxFIR Tx FIR
        %   Interpolation setting of FIR in transmit path. Possible options
        %   are [1,2,4]
        TxFIR = 2;
        %   Decimation setting of FIR in receive path. Possible options
        %   are [1,2,4]
        RxFIR = 2;
        %   Interpolation setting of HB1 in transmit path. Possible options
        %   are [1,2]
        TxHB1 = 2;
        %   Decimation setting of HB1 in receive path. Possible options
        %   are [1,2]
        RxHB1 = 2;
        %   Interpolation setting of HB2 in transmit path. Possible options
        %   are [1,2]
        TxHB2 = 2;
        %   Decimation setting of HB2 in receive path. Possible options
        %   are [1,2]
        RxHB2 = 2;
        %   Interpolation setting of HB3 in transmit path. Possible options
        %   are [1,2,3]. Decimation cannot be set to 3 when DACDivider is 2
        TxHB3 = 2;
        %   Decimation setting of HB3 in receive path. Possible options
        %   are [1,2,3]. Decimation cannot be set to 3 when DACDivider is 2
        RxHB3 = 2;
        %DACDivider DAC Divider
        %   Clock divider for DAC, which divides the ADC clock by this
        %   value to set the DAC rate
        DACDivider = 1;
        %PLLDivider PLL Divider
        %   Clock divider for ADC, which divides the PLLRate by this value
        %   to generate the ADC rate
        PLLDivider = 8;
        %TxAnalogCuttoffFrequency Tx Analog Cuttoff Frequency
        %   3dB cuttoff of front-end analog filter for transmit path in Hz
        TxAnalogCuttoffFrequency = 3.389e6;
        %RxAnalogCuttoffFrequency Rx Analog Cuttoff Frequency
        %   3dB cuttoff of front-end analog filter for receive path in Hz
        RxAnalogCuttoffFrequency = 8.069e6;
    end
    
    properties (Dependent)
        %TxAnalogCuttoffFrequencyActual Tx Analog Cuttoff Frequency Actual (Read-only)
        %   Actual achievable 3dB analog cutoff filter for transmit path
        TxAnalogCuttoffFrequencyActual
        %RxAnalogCuttoffFrequencyActual Rx Analog Cuttoff Frequency Actual (Read-only)
        %   Actual achievable 3dB analog cutoff filter for receive path
        RxAnalogCuttoffFrequencyActual
        %ValidConfiguration Valid Configuration (Read-only)
        %   Logical that determines if the current configuration based on
        %   the halfband and clock settings is functional
        ValidConfiguration
        %AvailableFIRTaps Available FIR Taps (Read-only)
        %   Maximum number of available FIR taps to be programmed based on
        %   clock configuration
        AvailableFIRTaps
    end
    
    properties (Dependent, Hidden)
        
        RxRateHB1
        RxRateHB2
        RxRateHB3
        RxRateFIR
        
        TxRateHB1
        TxRateHB2
        TxRateHB3
        TxRateFIR
        
        PLLRate
        ADCRate
        DACRate
        
        CALDividerTx
        CALDividerRx
    end
    
    properties (Hidden)
        % Filter objects
        TIARx
        ChannelFilterRx
        
        FilterStagesRX
        FilterStagesTX
        
        HB1RxResponse
        HB2RxResponse
        HB3RxResponse
        HB3Dec3RxResponse
        
        HB1TxResponse
        HB2TxResponse
        HB3TxResponse
        HB3Int3TxResponse
        
        
        AnalogStagesRX
        AnalogStagesTX
        
        AnalogFiltNum
        AnalogFiltDen
        TIANum
        TIADen
        
        TxBBLPFNum
        TxBBLPFDen
        TxSecLPFNum
        TxSecLPFDen
        
    end
    
    properties (Constant, Hidden)
        % Min/Max tested path rates
        MaxBBPLLRate  = 1430000000;
        MinBBPLLRate  = 715000000;
        MaxADCRate    = 640000000;
        MinADCRate    = 715000000/(2^6); % 11.2 MHz
        MaxDACRate    = 640000000/2; % (MaxADCRate/2)
        MinDACRate    = 25000000; % (MaxADCRate/2)
        MaxDataRate   = 61440000;
        MaxDataRateAD9363   = 20e6;
        MinDataRate   = 715000000/(48*(2^6));
        MaxFIR        = 61440000*2;
        MaxRxHB1      = 245760000;
        MaxRxHB2      = 320000000;
        MaxRxHB3      = 640000000;
        MaxTxHB1      = 160000000;
        MaxTxHB2      = 320000000;
        MaxTxHB3      = 320000000;
        % Define the digital filters with fixed coefficients
        AllPassCoeff = 1;
        RxHB1Coeff = 2^(-11)*[-8 0 42 0 -147 0 619 1013 619 0 -147 0 42 0 -8];
        RxHB2Coeff = 2^(-8)*[-9 0 73 128 73 0 -9];
        RxHB3Coeff = 2^(-4)*[1 4 6 4 1];
        RxHB3CoeffDec3 = 2^(-14)*[55 83 0 -393 -580 0 1914 4041 5120 4041 1914 0 -580 -393 0 83 55];
        TxHB1Coeff = 2^(-14)*[-53 0 313 0 -1155 0 4989 8192 4989 0 -1155 0 313 0 -53];
        TxHB2Coeff = 2^(-8)*[-9 0 73 128 73 0 -9];
        TxHB3Coeff = 2^(-2)*[1 2 1];
        TxHB3CoeffInt3 = (1/3)*2^(-13)*[36 -19 0 -156 -12 0 479 223 0 -1215 -993 0 3569 6277 8192 6277 3569 0 -993 -1215 0 223 479 0 -12 -156 0 -19 36];
        % Halfband numerical behavior descriptions
        HB1ConfigRx = struct(...
            'FullPrecisionOverride',false,...
            'OutputDataType','Custom',...
            'CustomOutputDataType',numerictype([],16,14),...
            'CoefficientsDataType','Custom',...
            'CustomCoefficientsDataType',numerictype([],16),...
            'ProductDataType','Custom',...
            'AccumulatorDataType','Custom',...
            'CustomProductDataType',numerictype([],31,30),...
            'CustomAccumulatorDataType',numerictype([],33,30));
        HB1ConfigTx = struct(...
            'FullPrecisionOverride',false,...
            'OutputDataType','Custom',...
            'CustomOutputDataType',numerictype([],16,14),...
            'CoefficientsDataType','Custom',...
            'CustomCoefficientsDataType',numerictype([],16),...
            'ProductDataType','Custom',...
            'AccumulatorDataType','Custom',...
            'CustomProductDataType',numerictype([],31,29),...
            'CustomAccumulatorDataType',numerictype([],31,29));
        HB2ConfigRx = struct(...
            'FullPrecisionOverride',false,...
            'OutputDataType','Custom',...
            'CustomOutputDataType',numerictype([],16,14),...
            'CoefficientsDataType','Custom',...
            'CustomCoefficientsDataType',numerictype([],16),...
            'ProductDataType','Custom',...
            'AccumulatorDataType','Custom',...
            'CustomProductDataType',numerictype([],31,29),...
            'CustomAccumulatorDataType',numerictype([],32,29));
        HB2ConfigTx = struct(...
            'FullPrecisionOverride',false,...
            'OutputDataType','Custom',...
            'CustomOutputDataType',numerictype([],16,14),...
            'CoefficientsDataType','Custom',...
            'CustomCoefficientsDataType',numerictype([],16),...
            'ProductDataType','Custom',...
            'AccumulatorDataType','Custom',...
            'CustomProductDataType',numerictype([],31,29),...
            'CustomAccumulatorDataType',numerictype([],31,29));
        HB3ConfigRx = struct(...
            'FullPrecisionOverride',false,...
            'OutputDataType','Custom',...
            'CustomOutputDataType',numerictype([],8,6),...
            'CoefficientsDataType','Custom',...
            'CustomCoefficientsDataType',numerictype([],16),...
            'ProductDataType','Custom',...
            'AccumulatorDataType','Custom',...
            'CustomProductDataType',numerictype([],19,18),...
            'CustomAccumulatorDataType',numerictype([],21,18));
        HB3ConfigTx = struct(...
            'FullPrecisionOverride',false,...
            'OutputDataType','Custom',...
            'CustomOutputDataType',numerictype([],8,6),...
            'CoefficientsDataType','Custom',...
            'CustomCoefficientsDataType',numerictype([],16),...
            'ProductDataType','Custom',...
            'AccumulatorDataType','Custom',...
            'CustomProductDataType',numerictype([],19,17),...
            'CustomAccumulatorDataType',numerictype([],19,17));
        HB3ConfigRxDec3 = struct(...
            'FullPrecisionOverride',false,...
            'OutputDataType','Custom',...
            'CustomOutputDataType',numerictype([],16,14),...
            'CoefficientsDataType','Custom',...
            'CustomCoefficientsDataType',numerictype([],16),...
            'ProductDataType','Custom',...
            'AccumulatorDataType','Custom',...
            'CustomProductDataType',numerictype([],19,18),...
            'CustomAccumulatorDataType',numerictype([],21,18));
        HB3ConfigTxInt3 = struct(...
            'FullPrecisionOverride',false,...
            'OutputDataType','Custom',...
            'CustomOutputDataType',numerictype([],16,14),...
            'CoefficientsDataType','Custom',...
            'CustomCoefficientsDataType',numerictype([],16),...
            'ProductDataType','Custom',...
            'AccumulatorDataType','Custom',...
            'CustomProductDataType',numerictype([],19,18),...
            'CustomAccumulatorDataType',numerictype([],20,18));
    end
    
    methods
        % Constructor
        function obj = ad936x(varargin)
        end
        function delete(~)
        end
        
        %% Validate user tunable parameters
        function set.DataRate(obj, val)
            if ~strcmp(obj.ChipType,'AD9363') %#ok<MCSUP>
                validateattributes(val, {'numeric'}, ...
                    {'scalar', 'real','integer', 'positive', 'nonnan', 'finite','>=',obj.MinDataRate,'<=',obj.MaxDataRate}, ...
                    '', 'DataRate');
            else
                validateattributes(val, {'numeric'}, ...
                    {'scalar', 'real','integer', 'positive', 'nonnan', 'finite','>=',obj.MinDataRate,'<=',obj.MaxDataRateAD9363}, ...
                    '', 'DataRate');
                
            end
            obj.DataRate = val;
        end
        function set.RxAnalogCuttoffFrequency(obj, val)
            if ~strcmp(obj.ChipType,'AD9363') %#ok<MCSUP>
                validateattributes(val, {'numeric'}, ...
                    {'scalar', 'real','positive', 'nonnan', 'finite','>=',200000,'<=',56000000}, ...
                    '', 'RxAnalogCuttoffFrequency');
            else
                validateattributes(val, {'numeric'}, ...
                    {'scalar', 'real','positive', 'nonnan', 'finite','>=',200000,'<=',20000000}, ...
                    '', 'RxAnalogCuttoffFrequency');
            end
            obj.RxAnalogCuttoffFrequency = val;
        end
        function set.TxAnalogCuttoffFrequency(obj, val)
            if ~strcmp(obj.ChipType,'AD9363') %#ok<MCSUP>
                validateattributes(val, {'numeric'}, ...
                    {'scalar', 'real', 'positive', 'nonnan', 'finite','>=',200000,'<=',40000000}, ...
                    '', 'TxAnalogCuttoffFrequency');
            else
                validateattributes(val, {'numeric'}, ...
                    {'scalar', 'real', 'positive', 'nonnan', 'finite','>=',200000,'<=',20000000}, ...
                    '', 'TxAnalogCuttoffFrequency');
                
            end
            obj.TxAnalogCuttoffFrequency = val;
        end
        function set.RxFIR(obj, val)
            obj.paramCheck(val,[1,2,4],'RxFIR');
            obj.RxFIR = val;
        end
        function set.RxHB1(obj, val)
            obj.paramCheck(val,[1,2],'RxHB1');
            obj.RxHB1 = val;
        end
        function set.RxHB2(obj, val)
            obj.paramCheck(val,[1,2],'RxHB2');
            obj.RxHB2 = val;
        end
        function set.RxHB3(obj, val)
            obj.paramCheck(val,[1,2,3],'RxHB3');
            obj.RxHB3 = val;
        end
        function set.TxFIR(obj, val)
            obj.paramCheck(val,[1,2,4],'TxFIR');
            obj.TxFIR = val;
        end
        function set.TxHB1(obj, val)
            obj.paramCheck(val,[1,2],'TxHB1');
            obj.TxHB1 = val;
        end
        function set.TxHB2(obj, val)
            obj.paramCheck(val,[1,2],'TxHB2');
            obj.TxHB2 = val;
        end
        function set.TxHB3(obj, val)
            obj.paramCheck(val,[1,2,3],'TxHB3');
            obj.TxHB3 = val;
        end
        function set.DACDivider(obj, val)
            obj.paramCheck(val,[1,2],'DACDivider');
            obj.DACDivider = val;
        end
        function set.PLLDivider(obj, val)
            obj.paramCheck(val,[1,2,4,8,16,32,64],'PLLDivider');
            obj.PLLDivider = val;
        end
        
        %% Determine dependent parameters
        function value = get.RxRateFIR(obj)
            value = obj.DataRate * obj.RxFIR;
        end
        function value = get.RxRateHB1(obj)
            value = obj.RxRateFIR * obj.RxHB1;
        end
        function value = get.RxRateHB2(obj)
            value = obj.RxRateHB1 * obj.RxHB2;
        end
        function value = get.RxRateHB3(obj)
            value = obj.RxRateHB2 * obj.RxHB3;
        end
        function value = get.TxRateFIR(obj)
            value = obj.DataRate * obj.TxFIR;
        end
        function value = get.TxRateHB1(obj)
            value = obj.TxRateFIR * obj.TxHB1;
        end
        function value = get.TxRateHB2(obj)
            value = obj.TxRateHB1 * obj.TxHB2;
        end
        function value = get.TxRateHB3(obj)
            value = obj.TxRateHB2 * obj.TxHB3;
        end
        function value = get.ADCRate(obj)
            value = obj.RxRateHB3;
        end
        function value = get.PLLRate(obj)
            value = obj.ADCRate*obj.PLLDivider;
        end
        function value = get.DACRate(obj)
            value = obj.ADCRate/obj.DACDivider;
        end
        
        function caldiv = get.CALDividerRx(obj)
            div = ceil((obj.PLLRate/obj.RxAnalogCuttoffFrequency)*(log(2)/(2*pi)));
            caldiv = min(max(div,1),511);
            [~, caldiv] = calculate_rfbw(obj.PLLRate, caldiv, 'Rx', true);
        end
        function caldiv = get.CALDividerTx(obj)
            div = ceil((obj.PLLRate/obj.TxAnalogCuttoffFrequency)*(log(2)/(2*pi)));
            caldiv = min(max(div,1),511);
            [~, caldiv] = calculate_rfbw(obj.PLLRate, caldiv, 'Tx', true);
        end
        
        function value = get.RxAnalogCuttoffFrequencyActual(obj)
            value = calculate_rfbw(obj.PLLRate, obj.CALDividerRx, 'Rx', true);
        end
        function value = get.TxAnalogCuttoffFrequencyActual(obj)
            value = calculate_rfbw(obj.PLLRate, obj.CALDividerTx, 'Tx', true);
        end
        
        
        function value = get.AvailableFIRTaps(obj)
            % RX
            if obj.RxHB3 == 3 || obj.RxHB3 == 1
                Nrx = min(16*floor(obj.ADCRate/(obj.DataRate)),128);
            else
                Nrx = min(16*floor(obj.ADCRate/(2*obj.DataRate)),128);
            end
            % TX
            switch obj.TxFIR
                case 1
                    Nmax = 64;
                case 2
                    Nmax = 128;
                case 4
                    Nmax = 128;
                otherwise
                    error('Unknown FIR decimation/interpolation');
            end
            Ntx = min(16*floor(obj.DACRate*obj.DACDivider/(2*obj.DataRate)),Nmax);
            
            % Since Tx side is more limited set taps based on what is
            % available
            if Ntx < Nrx
                value = Ntx;
            else
                value = Nrx;
            end
        end
        
        function value = get.ValidConfiguration(obj)
            [~,c1] = validatePathRates(obj);
            c2 = obj.ADCRate == obj.DACRate*obj.DACDivider;
            c3 = obj.TxRateHB3 == obj.DACRate;
            value = c1 && c2 && c3;
        end
        
        
        function e = getConfigurationError(obj)
            [valid,v] = validatePathRates(obj);
            e = struct;
            e.msg = '';
            e.pass = true;
            if v
                return;
            else
                fn = fieldnames(valid);
                for k = 1:length(fn)
                    f = valid.(fn{k});
                    if  ~f.pass
                        e = f;
                        return
                    end
                end
            end
            if obj.ADCRate ~= obj.DACRate*obj.DACDivider
                e = struct;
                e.msg = 'ADCRate != DACRate*DACDivider';
                e.pass = false;
            end
            if obj.TxRateHB3 ~= obj.DACRate
                e = struct;
                e.msg = 'TxHB3Rate != DACRate';
                e.pass = false;
            end
        end
        
        function AutoSetRates(obj)
            %AutoSetRates Configures halfband filters, PLL clock, and
            %converter rates in an optimal configuration based on the
            %current 'DataRate' property
            
            currentMaxADCRate = 0;
            savedDACDivider = 0;
            rate_gov = false;
            savedFIR = 0;
            savedPLL = 0;
            savedFilterConfig = 0;
            
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
                    obj.RxFIR = FilterConfig(k,4);
                    obj.RxHB1 = FilterConfig(k,3);
                    obj.RxHB2 = FilterConfig(k,2);
                    obj.RxHB3 = FilterConfig(k,1);
                    obj.TxFIR = FilterConfig(k,4);
                    obj.TxHB1 = FilterConfig(k,3);
                    obj.TxHB2 = FilterConfig(k,2);
                    obj.TxHB3 = FilterConfig(k,1);
                    
                    % HB3 cannot be 3 if rate_gov enabled
                    if rate_gov && (obj.RxHB3==3)
                        continue;
                    end
                    
                    % Check valid rates (ignore some section pre-dacdiv)
                    valid = obj.validatePathRates();
                    v = valid.RxRateHB3.pass && valid.RxRateHB2.pass && valid.RxRateHB1.pass && valid.RxRateFIR.pass;
                    v = v && valid.TxRateHB2.pass && valid.TxRateHB1.pass && valid.TxRateFIR.pass;
                    if v
                        % Determine PLL divider
                        pll = obj.determine_pll_div();
                        if pll>0
                            obj.PLLDivider = pll;
                            % Determine DAC divider setting and check ADC/DAC settings
                            dac_div = obj.check_dac_adc_config(pll,k);
                            if dac_div>0
                                obj.DACDivider = dac_div;
                                if obj.ADCRate>currentMaxADCRate
                                    % Final check
                                    [~,v] = obj.validatePathRates();
                                    if ~v
                                        continue
                                    end
                                    currentMaxADCRate = obj.ADCRate;
                                    savedFilterConfig = k;
                                    savedDACDivider = obj.DACDivider;
                                    savedFIR = FIR;
                                    savedPLL = pll;
                                end
                            end
                        end
                    end
                end
            end
            % Set HBs based on best config found
            obj.PLLDivider = savedPLL;
            obj.RxFIR = savedFIR;
            obj.RxHB1 = FilterConfig(savedFilterConfig,3);
            obj.RxHB2 = FilterConfig(savedFilterConfig,2);
            obj.RxHB3 = FilterConfig(savedFilterConfig,1);
            obj.TxFIR = savedFIR;
            obj.TxHB1 = FilterConfig(savedFilterConfig,3);
            obj.TxHB2 = FilterConfig(savedFilterConfig,2);
            obj.TxHB3 = FilterConfig(savedFilterConfig,1);
            % Apply dac div
            if (savedFilterConfig<3) && (savedDACDivider>1)
                obj.TxHB1 = 1;
            elseif (savedFilterConfig<5) && (savedDACDivider>1)
                obj.TxHB2 = 1;
            end
        end
    end
    
    methods (Hidden)
        
        % Path Rate Validations
        function [valid,v] = validatePathRates(obj)
            valid = struct;
            valid.ADCRate = obj.rangeCheck(obj.ADCRate,'ADC Rate',obj.MinADCRate,obj.MaxADCRate);
            % Rates into stage
            valid.RxRateHB3 = obj.rangeCheck(obj.RxRateHB3,'Rx HB3',0,obj.MaxRxHB3);
            valid.RxRateHB2 = obj.rangeCheck(obj.RxRateHB2,'Rx HB2',0,obj.MaxRxHB2);
            valid.RxRateHB1 = obj.rangeCheck(obj.RxRateHB1,'Rx HB1',0,obj.MaxRxHB1);
            valid.RxRateFIR = obj.rangeCheck(obj.RxRateFIR,'Rx FIR',0,obj.MaxFIR);
            valid.DACRate = obj.rangeCheck(obj.DACRate,'DAC Rate',0,obj.MaxDACRate);
            % Rates out of stage
            valid.TxRateHB3 = obj.rangeCheck(obj.TxRateHB3,'Tx HB3',0,obj.MaxTxHB3);
            valid.TxRateHB2 = obj.rangeCheck(obj.TxRateHB2,'Tx HB2',0,obj.MaxTxHB2);
            valid.TxRateHB1 = obj.rangeCheck(obj.TxRateHB1,'Tx HB1',0,obj.MaxTxHB1);
            valid.TxRateFIR = obj.rangeCheck(obj.TxRateFIR,'Tx FIR',0,obj.MaxFIR);
            valid.PLLRate = obj.rangeCheck(obj.PLLRate,'PLL Rate',obj.MinBBPLLRate,obj.MaxBBPLLRate);
            v = valid.ADCRate.pass && valid.RxRateHB3.pass && valid.RxRateHB2.pass && valid.RxRateHB1.pass && valid.RxRateFIR.pass;
            v = v && valid.DACRate.pass && valid.TxRateHB3.pass && valid.TxRateHB2.pass && valid.TxRateHB1.pass && valid.TxRateFIR.pass;
            v = v && valid.PLLRate.pass;
        end

    end
    
    methods (Access = protected, Hidden)
        
        function pll = determine_pll_div(obj)
            % Determine necessary PLL multiplier
            PLL_mult = 64; %MAX_BBPLL_DIV;
            rate1 = obj.ADCRate; % CHECK THIS RATE IS CORRECT
            
            while (PLL_mult>1)
                rate0 = rate1*PLL_mult;
                v = (rate0 >= obj.MinBBPLLRate) && (rate0 <= obj.MaxBBPLLRate);
                if v
                    pll = PLL_mult;
                    return
                end
                PLL_mult = PLL_mult/2;
            end
            pll = -1;
        end
        function r = check_dac_adc_config(obj,PLL_mult,dec_table_index)
            
            with_dd = obj.PLLRate/PLL_mult/2;
            without_dd = obj.PLLRate/PLL_mult/1;
            
            a = obj.rangeCheck(with_dd,'DAC Rate',obj.MinDACRate,obj.MaxDACRate);
            b = obj.rangeCheck(without_dd,'ADC Rate',obj.MinADCRate,obj.MaxADCRate);
            c = obj.rangeCheck(without_dd,'DAC Rate',obj.MinDACRate,obj.MaxDACRate);
            
            if (c.pass && b.pass)
                r = 1; %Run without dac div
            elseif (a.pass && b.pass && (dec_table_index<6))
                r = 2; % Run with dac div
            else
                r = -1; % All rates invalid
            end
        end
        
    end
    
    methods (Static, Hidden)
        function PV = getPVPair(s)
            F = fieldnames(s);
            C = struct2cell(s);
            PV = [reshape(F, 1, []); reshape(C, 1, [])];
        end
        
        % Generic parameter check
        function paramCheck(val,possible,name)
            if ~sum(val==possible)
                error('%s must be one of %s',name,num2str(possible));
            end
        end
        
        % coerces the normalized cutoff frequency passed between 0.0 and 1.0
        % for digital Butterworth filter designs
        function Wn = coerce_cutoff(freq)
            Wn = freq;
            if Wn < 0.0
                Wn = 0.0 + eps;
            elseif Wn > 1.0
                Wn = 1.0 - eps;
            end
        end
        
        function e = rangeCheck(val,name,min,max)
            % scalar', 'real','integer', 'nonnan'
            if (fix(val)~=val) || ~isreal(val) || isnan(val) || ~isscalar(val)
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
        
        
    end
end



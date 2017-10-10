classdef ad936xFilterDesigner < handle
    %ad936xFilterDesigner Base Class
    %
    
    properties
        % Chip specific configurations
        HB1 = 2;
        HB2 = 2;
        HB3 = 2;
        FIR = 2;
        FPGA_Filter = false; % Pluto FPGA filter enabled
        Rdata = 7680000; % Data Rate at chip output
        PLL_mult = 8;
        converter_rate = 122880000;
        PLL_rate = 983040000;
        caldiv = 25;
        RFbw = 6196961;
        FIRdBmin = 0;
        % Design criteria
        Fpass = 2560000;
        Fstop = 3200000;
        Fcenter = 0;
        Apass = 0.5000;
        Astop = 80;
        phEQ = -1;
        wnom = 4480000;
        Type = 'Lowpass';
        RxTx = 'Rx';
        % Results from designer
        GroupDelay = -1;
        GroupDelayVariance = -1;
        ApassGain = -1;
        AstopGain = -1;
        MaxInput = -1;
        NumTaps = -1;
        Taps = [];
        
    end
    
    properties (Constant, Access = private)
        AD9361_MAX_RATE = 61440000;
        AD9361_MIN_RATE = 520833;
        AD9361_MAX_RX_HB1 = 245760000;
        AD9361_MAX_RX_HB2 = 320000000;
        AD9361_MAX_RX_HB3 = 640000000;
        AD9361_MAX_TX_HB1 = 160000000;
        AD9361_MAX_TX_HB2 = 320000000;
        AD9361_MAX_TX_HB3 = 320000000;
        AD9361_MAX_FIR = 61440000 * 2
        AD9361_MAX_BBPLL_FREQ = 1430000000;
        AD9361_MIN_BBPLL_FREQ = 715000000;
        AD9361_MIN_ADC_CLK = 25000000;
        AD9361_MAX_ADC_CLK = 640000000;
        AD9361_MIN_DAC_CLK = 25000000;
        AD9361_MAX_DAC_CLK = 640000000 / 2
    end
    
    
    methods
        % Constructor
        function obj = ad936xFilterDesigner(varargin)
%                     setProperties(obj,nargin,varargin{:});
        end
        %% Check properties
        
        function set.Rdata(obj, value)
            validateattributes( value, { 'double','single' }, ...
                {'>=', 0.1e6, '<=', 61.44e6,'real', 'positive', 'scalar', 'finite',...
                'nonnan', 'nonempty'}, ...
                '', 'RData');
            obj.Rdata = value;
        end
        function set.HB1(obj, value)
            validateattributes( value, { 'double','single' }, ...
                {'>=', 1, '<=', 2,'real', 'positive', 'scalar', 'finite',...
                'nonnan', 'nonempty', 'integer'}, ...
                '','HB1');
            obj.HB1 = double(value);
        end
        function set.HB2(obj, value)
            validateattributes( value, { 'double','single' }, ...
                {'>=', 1, '<=', 2,'real', 'positive', 'scalar', 'finite',...
                'nonnan', 'nonempty', 'integer'}, ...
                '','HB2');
            obj.HB2 = double(value);
        end
        function set.HB3(obj, value)
            validateattributes( value, { 'double','single' }, ...
                {'>=', 1, '<=', 3,'real', 'positive', 'scalar', 'finite',...
                'nonnan', 'nonempty', 'integer'}, ...
                '','HB3');
            obj.HB3 = double(value);
        end
        function set.FIR(obj, value)
            validateattributes( value, { 'double','single' }, ...
                {'real', 'positive', 'scalar', 'finite',...
                'nonnan', 'nonempty', 'integer'}, ...
                '','FIR');
            vals = 2.^(0:2);
            if ~sum(vals==value)
                error(['Expected FIR to be any of: ',num2str(vals)]);
            end
            isa(value,vals);
            obj.FIR = double(value);
        end
        function set.PLL_mult(obj, value)
            validateattributes( value, { 'double','single' }, ...
                {'real', 'positive', 'scalar', 'finite',...
                'nonnan', 'nonempty', 'integer'}, ...
                '','PLL_mult');
            vals = 2.^(0:6);
            if ~sum(vals==value)
                error(['Expected PLL_mult to be any of: ',num2str(vals)]);
            end
            isa(value,vals);
            obj.PLL_mult = double(value);
        end
        
    end
    
    methods (Access = private)
        
        %% Check path rates
        function ValidatePathRates(obj)
            
        end
        
        %% Call designer
        function DesignFilter(obj)
        
        end
    end
end

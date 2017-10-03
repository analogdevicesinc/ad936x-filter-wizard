classdef ad936xFilterDesigner < handle
    %ad936xFilterDesigner Base Class
    %
    
    properties
%        TX = struct('HB1',1,'HB2',1);
%        RX
        HB1 = 1;
        HB2 = 1;
        HB3 = 1;
        FIR = 1;
        FPGA_Filter = false;
        RData = 30.72e6;
    end
    
    
    
    methods
%         % Constructor
%         function obj = ad936xFilterDesigner(varargin)
%             setProperties(obj,nargin,varargin{:});
%         end
        %% Check properties
        
        % Data rate
        function set.RData(obj, value)
            validateattributes( value, { 'double','single' }, ...
                {'>', 0.1, '<=', 61.44,'real', 'positive', 'scalar', 'finite',...
                'nonnan', 'nonempty'}, ...
                '', 'RData');
            % Scale to MHZ
            obj.RData = value*1e6;
        end
        
%         %% Utility function
%         function autoselect_rates(obj)
%             % sanity check the PLL rate and DAC divider values and alter them if necessary
%             if isfield(handles, 'tx') && isfield(handles, 'rx')
%                 if (obj.rx.PLL_rate ~= obj.tx.PLL_rate)
%                     % If dec3 is used by Rx or Tx, both must use it in order for the
%                     % PLL rates to match.
%                     if obj.tx.HB3 == 3 || obj.rx.HB3 == 3
%                         obj.rx.HB3 = 3;
%                         obj.tx.HB3 = 3;
%                     end
%                     
%                     ADC_rate = obj.rx.Rdata * obj.rx.FIR * ...
%                         obj.rx.HB1 * obj.rx.HB2 * obj.rx.HB3;
%                     DAC_rate = obj.tx.Rdata * obj.tx.FIR * ...
%                         obj.tx.HB1 * obj.tx.HB2 * obj.tx.HB3;
%                     DAC_div = ADC_rate / DAC_rate;
%                     if (obj.tx.DAC_div ~= DAC_div)
%                         if (DAC_div == 1 || DAC_div == 2)
%                             obj.tx.DAC_div = DAC_div;
%                             obj.tx.PLL_mult = obj.rx.PLL_mult;
%                             filter_type = get(handles.filter_type, 'Value');
%                             set(handles.filter_type, 'Value', 0);
%                             obj.tx.caldiv = default_caldiv(handles);
%                             set(handles.filter_type, 'Value', filter_type);
%                         end
%                     end
%                     
%                     obj.rx.PLL_mult = fastest_FIR([64 32 16 8 4 2 1], handles.bounds.MAX_BBPLL_FREQ, handles.bounds.MIN_BBPLL_FREQ, ...
%                         obj.rx.Rdata * obj.rx.FIR * obj.rx.HB1 * obj.rx.HB2 * obj.rx.HB3 * obj.rx.DAC_div);
%                     obj.tx.PLL_mult = obj.rx.PLL_mult;
%                     
%                     if obj.rx.PLL_mult > 64
%                         X = ['Date rate = ', num2str(tohwTx.TXSAMP), ' Hz. Tx BBPLL is too high for Rx to match.'];
%                         disp(X);
%                     end
%                     
%                     obj.rx.PLL_rate = obj.rx.Rdata * obj.rx.FIR * obj.rx.HB1 * ...
%                         obj.rx.HB2 * obj.rx.HB3 * obj.rx.PLL_mult;
%                 else
%                     ADC_rate = obj.rx.Rdata * obj.rx.FIR * ...
%                         obj.rx.HB1 * obj.rx.HB2 * obj.rx.HB3;
%                     DAC_rate = obj.tx.Rdata * obj.tx.FIR * ...
%                         obj.tx.HB1 * obj.tx.HB2 * obj.tx.HB3;
%                     DAC_div = ADC_rate / DAC_rate;
%                     if (obj.tx.DAC_div ~= DAC_div)
%                         if (DAC_div == 1 || DAC_div == 2)
%                             obj.tx.DAC_div = DAC_div;
%                             obj.tx.PLL_mult = obj.rx.PLL_mult;
%                             filter_type = get(handles.filter_type, 'Value');
%                             set(handles.filter_type, 'Value', 0);
%                             obj.tx.caldiv = default_caldiv(handles);
%                             set(handles.filter_type, 'Value', filter_type);
%                         end
%                     end
%                 end
%             end
%             
%         end
        
    end
end

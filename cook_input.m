%%
% Copyright 2014(c) Analog Devices, Inc.
%
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without modification,
%  are permitted provided that the following conditions are met:
%      - Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      - Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the
%        distribution.
%      - Neither the name of Analog Devices, Inc. nor the names of its
%        contributors may be used to endorse or promote products derived
%        from this software without specific prior written permission.
%      - The use of this software may or may not infringe the patent rights
%        of one or more patent holders.  This license does not release you
%        from the requirement that you obtain separate licenses from these
%        patent holders to use this software.
%      - Use of the software either in source or binary form or filter designs
%        resulting from the use of this software, must be connected to, run
%        on or loaded to an Analog Devices Inc. component.
%
%  THIS SOFTWARE IS PROVIDED BY ANALOG DEVICES "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
%  INCLUDING, BUT NOT LIMITED TO, NON-INFRINGEMENT, MERCHANTABILITY AND FITNESS FOR A
%  PARTICULAR PURPOSE ARE DISCLAIMED.
%
%  IN NO EVENT SHALL ANALOG DEVICES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, INTELLECTUAL PROPERTY
%  RIGHTS, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
%  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
%  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
%  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%%
% structure in/out can be the same. A partial structure can be passed in,
% and a full structure is passed out.
%
% General description:
%   struct.Type     = The type of filter required. one of:
%                       - 'Lowpass' (default)
%                       - 'Bandpass'
%                       - 'Equalize'
%                       - 'Root Raised Cosine'
%   struct.RxTx     = Is this 'Rx' or 'Tx'? Default is 'Rx'.
%
% Clock Settings:
%   struct.Rdata    = The output datarate.
%   struct.FIR      = FIR interpolation/decimation rates [1 2 4]
%   struct.HB1      = HB1 interpolation/decimation rates [1 2]
%   struct.HB2      = HB2 interpolation/decimation rates [1 2]
%   struct.HB3      = HB3 interpolation/decimation rates [1 2 3]
%   struct.DAC_div  = the ADC/DAC ratio, for Rx channels, this is
%                    always '1', for Tx, it is either '1' or '2'
%   struct.PLL_mult = the converter to PLL ratio
%   struct.PLL_rate = the PLL rate in Hz
%
% Description of the FIR
%   struct.Fpass    = Passband Frequency in Hz
%   struct.Fstop    = Stopband Frequency in Hz
%   struct.Fcenter  = Center Frequency in Hz (only used for Bandpass),
%                     otherwise 0
%   struct.Apass    = Passband ripple in dB (peak to peak)
%   struct.Astop    = Cascaded (FIR + HB + Analog) stop band attenuation (in dB)
%   struct.FIRdBmin = Minimum stop band attentuation of the FIR (in dB)
%                     un-cascaded. 0 if not used.
%
% Description of the Analog Filter settings
%   struct.caldiv   = The actual discrete register value that describes the
%                     rolloff for the analog filters
%   struct.Fcutoff  = the -3dB point of the Analog Filters expressed in
%                     baseband frequency (Hz)
%   struct.wnom     = the RF bandwidth of the Analog Filters
%   struct.phEQ     = the target for phase equalization in nanoseconds
%                     (-1 for none).

function cooked = cook_input(input)

% import various value bounds9361
bounds9361 = rate_bounds;

if ~isstruct(input)
    input = struct;
end

if ~isfield(input, 'Type')
    input.Type = 'Lowpass';
end

if (isfield(input, 'RxTx') && ~(strcmp(input.RxTx, 'Rx') || strcmp(input.RxTx, 'Tx')))
    error('RxTx parameter must be Rx or Tx');
end

if ~isfield(input, 'RxTx')
    input.RxTx = 'Rx';
end

% Make sure all the clock settings are there.
%   struct.Rdata    = The output datarate.
%   struct.FIR      = FIR interpolation/decimation rates [1 2 4]
%   struct.HB1      = HB1 interpolation/decimation rates [1 2]
%   struct.HB2      = HB2 interpolation/decimation rates [1 2]
%   struct.HB3      = HB3 interpolation/decimation rates [1 2 3]
%   struct.DAC_div  = the ADC/DAC ratio, for Rx channels, this is
%                    always '1', for Tx, it is either '1' or '2'
%   struct.PLL_mult = the converter to PLL ratio

if ~isfield(input, 'Rdata')
    if isfield(input, 'PLL_rate')
        input.Rdata = input.PLL_rate;
        while input.Rdata > bounds9361.MAX_DATA_RATE / 2
            input.Rdata = input.Rdata / 2;
        end
    else
        % Assume LTE5
        input.Rdata = 7680000;
    end
else
    if ~isfloat(input.Rdata)
        error('Rdata must be a floating point value!');
    end
end

if input.Rdata > bounds9361.MAX_DATA_RATE
    input.Rdata = bounds9361.MAX_DATA_RATE;
end
if input.Rdata < bounds9361.MIN_DATA_RATE
    input.Rdata = bounds9361.MIN_DATA_RATE;
end

input = autoselect_rates(input, bounds9361, false);

pll_out_of_bounds = ((input.PLL_rate > bounds9361.MAX_BBPLL_FREQ) || (input.PLL_rate < bounds9361.MIN_BBPLL_FREQ));
converter_out_of_bounds = 0;
if strcmp(input.RxTx, 'Rx')
    converter_out_of_bounds = (input.converter_rate > bounds9361.MAX_ADC_CLK || input.converter_rate < bounds9361.MIN_ADC_CLK);
elseif strcmp(input.RxTx, 'Tx')
    converter_out_of_bounds = (input.converter_rate > bounds9361.MAX_DAC_CLK || input.converter_rate < bounds9361.MIN_DAC_CLK);
end

% If PLL or converter rate bounds9361 aren't met, enable 3x dec/int for HB3.
if (pll_out_of_bounds || converter_out_of_bounds)
    input = autoselect_rates(input, bounds9361, true);
end

if strcmp(input.Type, 'Lowpass')
    if ~isfield(input, 'Fpass')
        % Asssume that Fpass is 1/3 datarate, which is about right for LTE5
        % works out to 2560000. Actual number is 2250000
        input.Fpass = input.Rdata / 3;
    end
    
    if ~isfield(input, 'Fstop')
        % Asssume that Fstop is 1.25 Fpass, again close to LTE5
        input.Fstop = input.Fpass * 1.25;
    end
    
    if ~isfield(input, 'Fcenter')
        input.Fcenter = 0;
    end
elseif strcmp(input.Type, 'Bandpass')
    error('Bandpass is not done yet');
end

% struct.Apass = Passband ripple (Apass) in dB (peak to peak)
% struct.Astop   = Cascaded (FIR + HB + Analog) stop band attenuation (in dB)
if ~isfield(input, 'Apass')
    input.Apass = .5;
end

if ~isfield(input, 'Astop')
    input.Astop = 80;
end

% Assume no phase equalization
if ~isfield(input, 'phEQ')
    input.phEQ = -1;
end

% Assume RF bandwidth (nominal frequency)
if ~isfield(input, 'wnom')
    if strcmp(input.RxTx, 'Rx')
        input.wnom = 1.4 * input.Fstop; % Rx
    else
        input.wnom = 1.6 * input.Fstop; % Tx
    end
end

if ~isfield(input, 'caldiv')
    div = ceil((input.PLL_rate/input.wnom)*(log(2)/(2*pi)));
    input.caldiv = min(max(div,1),511);
end

if ~isfield(input, 'RFbw')
    [input.RFbw, input.caldiv] = calculate_rfbw(input.PLL_rate, input.caldiv, input.RxTx, false);
end

% Assume no dBmin
if ~isfield(input, 'FIRdBmin')
    input.FIRdBmin = 0;
end

cooked = input;

function input = autoselect_rates(input, bounds9361, dec_int3)
if strcmp(input.RxTx, 'Rx')
    max_HB = bounds9361.MAX_RX;
else
    max_HB = bounds9361.MAX_TX;
end

if ~isfield(input, 'DAC_div')
    if strcmp(input.RxTx, 'Rx')
        input.DAC_div = 1;
    else
        input.DAC_div = 2;
    end
end

if dec_int3 || (~isfield(input, 'FIR') && ~isfield(input, 'HB1') && ~isfield(input, 'HB2') && ~isfield(input, 'HB3') && ~isfield(input, 'PLL_mult'))
    % Everything is blank, run as fast as possible
    if dec_int3
        input.HB3 = fastest_FIR([3 2 1], max_HB.HB3, 0, input.Rdata);
    else
        input.HB3 = fastest_FIR([2 1], max_HB.HB3, 0, input.Rdata);
    end
    input.HB2 = fastest_FIR([2 1], max_HB.HB2, 0, input.Rdata * input.HB3);
    input.HB1 = fastest_FIR([2 1], max_HB.HB1, 0, input.Rdata * input.HB3 * input.HB2);
    input.FIR = fastest_FIR([4 2 1], bounds9361.MAX_FIR, 0, input.Rdata * input.HB3 * input.HB2 * input.HB1);
    input.PLL_mult = fastest_FIR([64 32 16 8 4 2 1], bounds9361.MAX_BBPLL_FREQ, bounds9361.MIN_BBPLL_FREQ, input.Rdata * input.FIR * input.HB1 * input.HB2 * input.HB3 * input.DAC_div);
end

input.converter_rate = input.Rdata * input.FIR * input.HB1 * input.HB2 * input.HB3;
input.PLL_rate = input.converter_rate * input.DAC_div * input.PLL_mult;

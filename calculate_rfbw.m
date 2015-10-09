%  Copyright 2015(c) Analog Devices, Inc.
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

% calculate a channel's complex bandwidth related to the calibration divider value
function [rfbw, caldiv] = calculate_rfbw(pll_rate, caldiv, RxTx, hw)
if strcmp(RxTx, 'Rx')
    channel_factor = 1.4;
    % (1.4 * 2 * pi)/log(2) rounded to the same precision the driver uses
    rounded_factor = 12.6906;
else
    channel_factor = 1.6;
    % (1.6 * 2 * pi)/log(2) rounded to the same precision the driver uses
    rounded_factor = 14.5036;
end

if hw
    % avoid divide by zero on boundary case
    if caldiv == 1
        caldiv = 1 + eps;
    end
    % used to reproduce the divider value (caldiv) we expect on the driver
    rfbw = uint32(fix(((pll_rate - 1)/(caldiv - 1))*(2/rounded_factor)));
else
    % full precision RF bandwidth
    rfbw = round((pll_rate/caldiv)*(2/(channel_factor*(2*pi)/log(2))));
end

% min/max possible values for the RF bandwidth (2x baseband bandwidth) from the
% reference manual (values are in Hz since RFbw is in Hz)
if strcmp(RxTx, 'Rx')
    % Rx: 0.4 MHz <= RF bandwidth <= 56 MHz
    min_rfbw = 400000;
    max_rfbw = 56000000;
else
    % Tx: 1.25 MHz <= RF bandwidth <= 40 MHz
    min_rfbw = 1250000;
    max_rfbw = 40000000;
end

% If the RF bandwidth is outside the range of acceptable values we modify
% the divider value until it falls into an acceptable range.
while (rfbw < min_rfbw) || (rfbw > max_rfbw)
    if (rfbw < min_rfbw)
        caldiv = caldiv - 1;
    else
        caldiv = caldiv + 1;
    end

    if (caldiv < 1) || (caldiv > 511)
        msgbox(sprintf('Calibration divider out of bounds (1 - 511): %i', caldiv), 'Error', 'error');
        return;
    end

    rfbw = calculate_rfbw(pll_rate, caldiv, RxTx, hw);
end

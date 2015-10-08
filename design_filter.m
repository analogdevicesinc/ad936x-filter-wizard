%  Copyright 2014-2015(c) Analog Devices, Inc.
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

% Inputs (structure containing the following fields)
% ============================================
% Rdata      = input/output sample data rate (in Hz)
% FIR        = FIR interpolation/decimation factor
% PLL_mult   = PLL multiplication
% Fpass      = passband frequency (in Hz)
% Fstop      = stopband frequency (in Hz)
% Apass      = max ripple allowed in passband (in dB)
% Astop      = min attenuation in stopband (in dB)
% FIRdBmin   = min rejection that FIR is required to have (in dB)
% phEQ       = phase equalization on (not -1)/off (-1)
% int_FIR    = use AD9361 FIR on (1)/off (0)
% wnom       = analog cutoff frequency (in Hz)
% clkPLL     = PLL frequency (in HZ)
%
% Outputs (structure containing the following fields)
% ===============================================
% firtaps          = fixed point FIR coefficients
% filter           = system object for visualization (does not include analog filters)
% Apass_actual     = actual passband ripple
% Astop_actual     = actual stopband attentuation
% delay            = actual delay used in phase equalization

function output = designfilter(input)

input = cook_input(input);
converter_rate = input.Rdata * input.FIR * input.HB1 * input.HB2 * input.HB3;

% use the internal FIR if unspecified
if ~isfield(input, 'int_FIR')
    input.int_FIR = 1;
end

if strcmp(input.RxTx, 'Rx')
    wTIA = input.wnom*(2.5/1.4);

    % Define the analog filters (for design purpose)
    [b1,a1] = butter(1,2*pi*wTIA,'s');  % 1st order
    [b2,a2] = butter(3,2*pi*input.wnom,'s');    % 3rd order

    % Digital representation of the analog filters (It is an approximation for group delay calculation only)
    [z1,p1,k1] = butter(3,coerce_cutoff(input.wnom/(converter_rate/2)),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    Hd1=dsp.BiquadFilter('SOSMatrix',sos1,'ScaleValues',g1);
    [z2,p2,k2] = butter(1,coerce_cutoff(wTIA/(converter_rate/2)),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    Hd2=dsp.BiquadFilter('SOSMatrix',sos2,'ScaleValues',g2);
    Hanalog = cascade(Hd2,Hd1);

    % Define the digital filters with fixed coefficients
    hb1_coeff = 2^(-11)*[-8 0 42 0 -147 0 619 1013 619 0 -147 0 42 0 -8];
    hb2_coeff = 2^(-8)*[-9 0 73 128 73 0 -9];
    hb3_coeff = 2^(-4)*[1 4 6 4 1];
    dec_int3_coeff = 2^(-14)*[55 83 0 -393 -580 0 1914 4041 5120 4041 1914 0 -580 -393 0 83 55];

    dec_int_func = @dsp.FIRDecimator;
else
    wreal = input.wnom*(5.0/1.6);

    % Define the analog filters (for design purpose)
    [b1,a1] = butter(3,2*pi*input.wnom,'s');     % 3rd order
    [b2,a2] = butter(1,2*pi*wreal,'s');  % 1st order

    % Digital representation of the analog filters (It is an approximation for group delay calculation only)
    [z1,p1,k1] = butter(3,coerce_cutoff(input.wnom/(converter_rate/2)),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    Hd1=dsp.BiquadFilter('SOSMatrix',sos1,'ScaleValues',g1);
    [z2,p2,k2] = butter(1,coerce_cutoff(wreal/(converter_rate/2)),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    Hd2=dsp.BiquadFilter('SOSMatrix',sos2,'ScaleValues',g2);
    Hanalog = cascade(Hd1,Hd2);

    % Define the digital filters with fixed coefficients
    hb1_coeff = 2^(-14)*[-53 0 313 0 -1155 0 4989 8192 4989 0 -1155 0 313 0 -53];
    hb2_coeff = 2^(-8)*[-9 0 73 128 73 0 -9];
    hb3_coeff = 2^(-2)*[1 2 1];
    dec_int3_coeff = (1/3)*2^(-13)*[36 -19 0 -156 -12 0 479 223 0 -1215 -993 0 3569 6277 8192 6277 3569 0 -993 -1215 0 223 479 0 -12 -156 0 -19 36];

    dec_int_func = @dsp.FIRInterpolator;
end

Hm1 = dec_int_func(2, hb1_coeff);
Hm2 = dec_int_func(2, hb2_coeff);
Hm3 = dec_int_func(2, hb3_coeff);
Hm4 = dec_int_func(3, dec_int3_coeff);

hb1 = input.HB1;
hb2 = input.HB2;
if input.HB3 == 2
    hb3 = 2;
    dec_int3 = 1;
elseif input.HB3 == 3
    hb3 = 1;
    dec_int3 = 3;
else
    hb3 = 1;
    dec_int3 = 1;
end

% convert the enables into a string
enables = strrep(num2str([hb1 hb2 hb3 dec_int3]), ' ', '');
if strcmp(input.RxTx, 'Rx')
    switch enables
        case '1111' % only FIR
            filter = 1;
        case '2111' % Hb1
            filter = Hm1;
        case '1211' % Hb2
            filter = Hm2;
        case '1121' % Hb3
            filter = Hm3;
        case '2211' % Hb2,Hb1
            filter = cascade(Hm2,Hm1);
        case '2121' % Hb3,Hb1
            filter = cascade(Hm3,Hm1);
        case '2221' % Hb3,Hb2,Hb1
            filter = cascade(Hm3,Hm2,Hm1);
        case '1113' % Dec3
            filter = Hm4;
        case '2113' % Dec3,Hb1
            filter = cascade(Hm4,Hm1);
        case '2213' % Dec3,Hb2,Hb1
            filter = cascade(Hm4,Hm2,Hm1);
        case '1221' % Hb3,Hb2
            filter = cascade(Hm3,Hm2);
        case '1213' % Dec3,Hb2
            filter = cascade(Hm4,Hm2);
        otherwise
            error('ddcresponse:IllegalOption', 'At least one of the stages must be there.')
    end
else
    switch enables
        case '1111' % only TFIR
            filter = 1;
        case '2111' % Hb1
            filter = Hm1;
        case '1211' % Hb2
            filter = Hm2;
        case '1213' % Hb2,Int3
            filter = cascade(Hm2,Hm4);
        case '1121' % Hb3
            filter = Hm3;
        case '1221' % Hb2,Hb3
            filter = cascade(Hm2,Hm3);
        case '2211' % Hb1,Hb2
            filter = cascade(Hm1,Hm2);
        case '2121' % Hb1,Hb3
            filter = cascade(Hm1,Hm3);
        case '2221' % Hb1,Hb2,Hb3
            filter = cascade(Hm1,Hm2,Hm3);
        case '1113' % Int3
            filter = Hm4;
        case '2113' % Hb1,Int3
            filter = cascade(Hm1,Hm4);
        case '2213' % Hb1,Hb2,Int3
            filter = cascade(Hm1,Hm2,Hm4);
        otherwise
            error('ddcresponse:IllegalOption', 'At least one of the stages must be there.')
    end
end

Hmiddle = clone(filter);
if strcmp(input.RxTx, 'Rx')
    if strcmp(enables,'1111') ||  strcmp(enables,'2111') || strcmp(enables,'1211') || strcmp(enables,'1121') || strcmp(enables,'1113')
        Hmiddle = cascade(Hd1,Hmiddle);
    else
        addStage(Hmiddle,Hd1,1);
    end
    addStage(Hmiddle,Hd2,1);
else
    if strcmp(enables,'1111') ||  strcmp(enables,'2111') || strcmp(enables,'1211') || strcmp(enables,'1121') || strcmp(enables,'1113')
        Hmiddle = cascade(Hmiddle,Hd1);
    else
        addStage(Hmiddle,Hd1);
    end
    addStage(Hmiddle,Hd2);
end

% Find out the best fit delay on passband
Nw = 2048;
w = zeros(1,Nw);
phi = zeros(1,Nw);
invariance = zeros(1,Nw);

w(1) = -input.Fpass;
for i = 2:(Nw)
    w(i) = w(1)-2*w(1)*i/(Nw);
end

if strcmp(input.RxTx, 'Rx')
    response = analogresp('Rx',w,converter_rate,b1,a1,b2,a2).*freqz(filter,w,converter_rate);
else
    response = freqz(filter,w,converter_rate).*analogresp('Tx',w,converter_rate,b1,a1,b2,a2);
end
for i = 1:(Nw)
    invariance(i) = real(response(i))^2+imag(response(i))^2;
end

phi(1)=atan2(imag(response(1)),real(response(1)));
for i = 2:(Nw)
    phi(i) = phi(i-1)+alias_b(atan2(imag(response(i)),real(response(i)))-phi(i-1),2*pi);
end

sigma = sum(invariance);
sigmax = sum(w.*invariance);
sigmay = sum(phi.*invariance);
sigmaxx = sum(w.*w.*invariance);
sigmaxy = sum(w.*phi.*invariance);
delta = sigma*sigmaxx-sigmax^2;
b = (sigma*sigmaxy-sigmax*sigmay)/delta;
if input.phEQ == 0 || input.phEQ == -1
    delay = -b/(2*pi);
else
    delay = input.phEQ*(1e-9);
end

% Design the FIR
G = 16384;
clkFIR = input.Rdata*input.FIR;
Gpass = floor(G*input.Fpass/clkFIR);
Gstop=ceil(G*input.Fstop/clkFIR);
Gpass = min(Gpass,Gstop-1);
fg = zeros(1,Gpass);
omega = zeros(1,Gpass);

% passband
for i = 1:(Gpass+1)
    fg(i) = (i-1)/G;
    omega(i) = fg(i)*clkFIR;
end
if strcmp(input.RxTx, 'Rx')
    rg1 = analogresp('Rx',omega,converter_rate,b1,a1,b2,a2).*freqz(filter,omega,converter_rate);
else
    rg1 = freqz(filter,omega,converter_rate).*analogresp('Tx',omega,converter_rate,b1,a1,b2,a2);
end
phase = unwrap(angle(rg1));
gd1 = GroupDelay(omega,phase); % group delay on passband for Analog + Converter + HB
omega1 = omega;                % frequency grid on pass band
rg2 = exp(-1i*2*pi*omega*delay);
rg = rg2./rg1;
w = abs(rg1)/(dBinv(input.Apass/2)-1);

g = Gpass+1;
% stop band
for m = Gstop:(G/2)
    g = g+1;
    fg(g) = m/G;
    omega(g) = fg(g)*clkFIR;
    rg(g) = 0;
end
if strcmp(input.RxTx, 'Rx')
    wg1 = abs(analogresp('Rx',omega(Gpass+2:end),converter_rate,b1,a1,b2,a2).*freqz(filter,omega(Gpass+2:end),converter_rate));
    wg2 = (wg1)/(dBinv(-input.Astop));
else
    wg1 = abs(freqz(filter,omega(Gpass+2:end),converter_rate).*analogresp('Tx',omega(Gpass+2:end),converter_rate,b1,a1,b2,a2));
    wg2 = (sqrt(input.FIR)*wg1)/(dBinv(-input.Astop));
end
wg3 = dBinv(input.FIRdBmin);
wg = max(wg2,wg3);
grid = fg;
if input.phEQ == -1
    resp = abs(rg);
else resp = rg;
end
weight = [w wg];
weight = weight/max(weight);

% design FIR filter
cr = real(resp);
B = 2;
F1 = grid(1:Gpass+1)*2;
F2 = grid(Gpass+2:end)*2;
A1 = cr(1:Gpass+1);
A2 = cr(Gpass+2:end);
W1 = weight(1:Gpass+1);
W2 = weight(Gpass+2:end);

% Determine the number of taps for FIR
if strcmp(input.RxTx, 'Rx')
    if hb3 == 1
        N = min(16*floor(converter_rate/(input.Rdata)),128);
    else
        N = min(16*floor(converter_rate/(2*input.Rdata)),128);
    end
else
    switch input.FIR
        case 1
            Nmax = 64;
        case 2
            Nmax = 128;
        case 4
            Nmax = 128;
    end
    N = min(16*floor(converter_rate*input.DAC_div/(2*input.Rdata)),Nmax);
end

tap_store = zeros(N/16,N);
Apass_actual_vector = zeros(N/16,1);
Astop_actual_vector = zeros(N/16,1);
i = 1;

while (1)
    if input.int_FIR
        d = fdesign.arbmag('N,B,F,A',N-1,B,F1,A1,F2,A2);
    else
        d = fdesign.arbmag('B,F,A,R');
        d.NBands = 2;
        d.B1Frequencies = F1;
        d.B1Amplitudes = A1;
        d.B1Ripple = db2mag(-input.Astop);
        d.B2Frequencies = F2;
        d.B2Amplitudes = A2;
        d.B2Ripple = db2mag(-input.Astop);
    end
    Hd = design(d,'equiripple','B1Weights',W1,'B2Weights',W2,'SystemObject',false);
    ccoef = Hd.Numerator;
    M = length(ccoef);

    if input.phEQ ~= -1
        sg = 0.5-grid(end:-1:1);
        sr = imag(resp(end:-1:1));
        sw = weight(end:-1:1);
        F3 = sg(1:G/2-Gstop+1)*2;
        F4 = sg(G/2-Gstop+2:end)*2;
        A3 = sr(1:G/2-Gstop+1);
        A4 = sr(G/2-Gstop+2:end);
        W3 = sw(1:G/2-Gstop+1);
        W4 = sw(G/2-Gstop+2:end);
        if input.int_FIR
            d2 = fdesign.arbmag('N,B,F,A',N-1,B,F3,A3,F4,A4);
        else
            d2 = fdesign.arbmag('N,B,F,A',M-1,B,F3,A3,F4,A4);
        end
        Hd2 = design(d2,'equiripple','B1Weights',W3,'B2Weights',W4,'SystemObject',false);
        scoef = Hd2.Numerator;
        for k = 1:length(scoef)
            scoef(k) = -scoef(k)*(-1)^(k-1);
        end
    else
        scoef = 0;
    end
    tap_store(i,1:M)=ccoef+scoef;

    Hmd = dec_int_func(input.FIR,tap_store(i,1:M));
    if ~isempty(ver('fixedpoint'))
        Hmd.Numerator = double(fi(Hmd.Numerator,true,16));
    end
    if strcmp(input.RxTx, 'Rx')
        if strcmp(enables,'1111') || strcmp(enables,'2111') || strcmp(enables,'1211') || strcmp(enables,'1121') || strcmp(enables,'1113')
            filter = cascade(filter,Hmd);
        else
            addStage(filter,Hmd);
        end
        rg_pass = abs(analogresp('Rx',omega(1:Gpass+1),converter_rate,b1,a1,b2,a2).*freqz(filter,omega(1:Gpass+1),converter_rate));
        rg_stop = abs(analogresp('Rx',omega(Gpass+2:end),converter_rate,b1,a1,b2,a2).*freqz(filter,omega(Gpass+2:end),converter_rate));
    else
        if strcmp(enables,'1111') || strcmp(enables,'2111') || strcmp(enables,'1211') || strcmp(enables,'1121') || strcmp(enables,'1113')
            filter = cascade(Hmd,filter);
        else
            addStage(filter, Hmd, 1);
        end
        rg_pass = abs(freqz(filter,omega(1:Gpass+1),converter_rate).*analogresp('Tx',omega(1:Gpass+1),converter_rate,b1,a1,b2,a2));
        rg_stop = abs(freqz(filter,omega(Gpass+2:end),converter_rate).*analogresp('Tx',omega(Gpass+2:end),converter_rate,b1,a1,b2,a2));
    end

    % quantitative values about actual passband and stopband
    Apass_actual_vector(i) = mag2db(max(rg_pass))-mag2db(min(rg_pass));
    Astop_actual_vector(i) = -mag2db(max(rg_stop));

    if input.int_FIR == 0
        h = tap_store(1,1:M);
        Apass_actual = Apass_actual_vector(1);
        Astop_actual = Astop_actual_vector(1);
        if strcmp(input.RxTx, 'Rx')
            removeStage(filter);
        else
            removeStage(filter, 1);
        end
        break
    elseif Apass_actual_vector(1) > input.Apass || Astop_actual_vector(1) < input.Astop
        h = tap_store(1,1:N);
        Apass_actual = Apass_actual_vector(1);
        Astop_actual = Astop_actual_vector(1);
        if strcmp(input.RxTx, 'Rx')
            removeStage(filter);
        else
            removeStage(filter, 1);
        end
        break
    elseif Apass_actual_vector(i) > input.Apass || Astop_actual_vector(i) < input.Astop
        h = tap_store(i-1,1:N+16);
        Apass_actual = Apass_actual_vector(i-1);
        Astop_actual = Astop_actual_vector(i-1);
        if strcmp(input.RxTx, 'Rx')
            removeStage(filter);
        else
            removeStage(filter, 1);
        end
        break
    else
        N = N-16;
        i = i+1;
        if strcmp(input.RxTx, 'Rx')
            removeStage(filter);
        else
            removeStage(filter, 1);
        end
    end
end

if input.RxTx == 'Tx'
    if input.int_FIR == 1 && input.FIR == 2
        R = rem(length(h),32);
        if R ~= 0
            h = [zeros(1,8),h,zeros(1,8)];
        end
    elseif input.int_FIR == 1 && input.FIR == 4
        R = rem(length(h),64);
        if R ~= 0
            newlength = ceil(length(h)/64)*64;
            addlength = (newlength-length(h))/2;
            h = [zeros(1,addlength),h,zeros(1,addlength)];
        end
    end
end

Hmd = dec_int_func(input.FIR,h);
if ~isempty(ver('fixedpoint'))
    Hmd.Numerator = double(fi(Hmd.Numerator,true,16));
end
if strcmp(input.RxTx, 'Rx')
    addStage(filter, Hmd);
else
    addStage(filter, Hmd, 1);
end
gd2 = grpdelay(Hmd,omega1,clkFIR).*(1/clkFIR);
if input.phEQ == -1
    groupdelay = gd1 + gd2;
else
    groupdelay = gd1 + gd2';
end
grpdelayvar = max(groupdelay)-min(groupdelay);

aTFIR = 1 + ceil(log2(max(Hmd.Numerator)));
switch aTFIR
    case 2
        gain = 6;
    case 1
        gain = 0;
    case 0
        gain = -6;
    otherwise
        gain = -12;
end

if strcmp(input.RxTx, 'Rx')
    if aTFIR > 2
        gain = 6;
    end
else
    if input.FIR == 2
        gain = gain+6;
    elseif input.FIR == 4
        gain = gain+12;
    end
    if gain > 0
        gain = 0;
    elseif gain < -6
        gain = -6;
    end
end

bTFIR = 16 - aTFIR;
firtaps = Hmd.Numerator.*(2^bTFIR);

if length(firtaps) < 128
    firtaps = [firtaps,zeros(1,128-length(firtaps))];
end

output = input;

% webinar.Hm1_rx = Hm1;
% webinar.Hm2_rx = Hm2;
% webinar.Hm3_rx = Hm3;
% webinar.Hm4_rx = Hm4;
% webinar.Hmd_rx = Hmd;
% webinar.enable_rx = enables;
%
% tohw.RF = input.data_rate * input.FIR_interp;
% tohw.R1 = tohw.RF * input.HB1;
% tohw.R2 = tohw.R1 * input.HB2;
% tohw.ADC = input.converter_rate;
% tohw.BBPLL = input.clkPLL;

output.firtaps = firtaps;
output.taps_length = length(h);
output.filter = filter;
output.Hanalog = Hanalog;
output.Apass_actual = Apass_actual;
output.Astop_actual = Astop_actual;
output.delay = delay;
output.grpdelayvar = grpdelayvar;
output.gain = gain;
output.Hd1 = Hd1;
output.Hd2 = Hd2;
output.Hmd = Hmd;
output.Hmiddle = Hmiddle;
output.a1 = a1;
output.b1 = b1;
output.a2 = a2;
output.b2 = b2;

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
%
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
%
% Outputs (structure containing the following fields)
% ===============================================
% firtaps          = fixed point FIR coefficients
% filter           = system object for visualization (does not include analog filters)
% Apass_actual     = actual passband ripple
% Astop_actual     = actual stopband attentuation
% delay            = actual delay used in phase equalization

function output = internal_design_filter(input)

% support a simple data rate input otherwise it must be a structure
if isfloat(input)
    input = struct('Rdata', input);
end

input = cook_input(input);

% use the internal FIR if unspecified
if ~isfield(input, 'int_FIR')
    input.int_FIR = 1;
end

% nominal frequency can't be zero
if ~input.wnom
    input.wnom = double(calculate_rfbw(input.PLL_rate, input.caldiv, input.RxTx, true));
end

if strcmp(input.RxTx, 'Rx')
    wTIA = input.wnom*(2.5/1.4);
    
    % Define the analog filters (for design purpose)
    [b1,a1] = butter(1,2*pi*wTIA,'s');  % 1st order
    [b2,a2] = butter(3,2*pi*input.wnom,'s');    % 3rd order
    
    % Digital representation of the analog filters (It is an approximation for group delay calculation only)
    [z1,p1,k1] = butter(3,coerce_cutoff(input.wnom/(input.converter_rate/2)),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    Hd1=dsp.BiquadFilter('SOSMatrix',sos1,'ScaleValues',g1);
    [z2,p2,k2] = butter(1,coerce_cutoff(wTIA/(input.converter_rate/2)),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    Hd2=dsp.BiquadFilter('SOSMatrix',sos2,'ScaleValues',g2);
    Hanalog = cascade(Hd2,Hd1);
    
    % Define the Pluto DEC8 filter
    ast = 80;
    n = 128;
    f = fdesign.decimator(8, 'Nyquist', 8, 'N,Ast', n, ast);
    hf = design(f,'SystemObject',true);
    
    % Define the digital filters with fixed coefficients
    allpass_coeff = 1;
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
    [z1,p1,k1] = butter(3,coerce_cutoff(input.wnom/(input.converter_rate/2)),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    Hd1=dsp.BiquadFilter('SOSMatrix',sos1,'ScaleValues',g1);
    [z2,p2,k2] = butter(1,coerce_cutoff(wreal/(input.converter_rate/2)),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    Hd2=dsp.BiquadFilter('SOSMatrix',sos2,'ScaleValues',g2);
    Hanalog = cascade(Hd1,Hd2);
    
    % Define the Pluto INT8 filter
    ast = 80;
    n = 128;
    f = fdesign.interpolator(8,'Nyquist', 8,'N,Ast', n, ast);
    hf = design(f,'kaiserwin','SystemObject',true);
    hf.Numerator = hf.Numerator./8;
    
    % Define the digital filters with fixed coefficients
    allpass_coeff = 1;
    hb1_coeff = 2^(-14)*[-53 0 313 0 -1155 0 4989 8192 4989 0 -1155 0 313 0 -53];
    hb2_coeff = 2^(-8)*[-9 0 73 128 73 0 -9];
    hb3_coeff = 2^(-2)*[1 2 1];
    dec_int3_coeff = (1/3)*2^(-13)*[36 -19 0 -156 -12 0 479 223 0 -1215 -993 0 3569 6277 8192 6277 3569 0 -993 -1215 0 223 479 0 -12 -156 0 -19 36];
    
    dec_int_func = @dsp.FIRInterpolator;
end

Hallpass = dec_int_func(1, allpass_coeff);
Hm1 = dec_int_func(2, hb1_coeff);
Hm1.FullPrecisionOverride = false;
Hm1.OutputDataType='Custom';
Hm1.CustomOutputDataType=numerictype([],16,14);
Hm1.CoefficientsDataType='Custom';
Hm1.CustomCoefficientsDataType=numerictype([],16);
Hm1.ProductDataType='Custom';
Hm1.AccumulatorDataType = 'Custom';
if strcmp(input.RxTx, 'Rx')
    Hm1.CustomProductDataType=numerictype([],31,30);
    Hm1.CustomAccumulatorDataType=numerictype([],33,30);
else
    Hm1.CustomProductDataType=numerictype([],31,29);
    Hm1.CustomAccumulatorDataType=numerictype([],31,29);
end

Hm1c34 = dec_int_func(2, hb1_coeff);
Hm1c34.FullPrecisionOverride = false;
Hm1c34.OutputDataType='Custom';
Hm1c34.CustomOutputDataType=numerictype([],4,2);
Hm1c34.CoefficientsDataType='Custom';
Hm1c34.CustomCoefficientsDataType=numerictype([],16);
Hm1c34.ProductDataType='Custom';
Hm1c34.AccumulatorDataType = 'Custom';
if strcmp(input.RxTx, 'Rx')
    Hm1c34.CustomProductDataType=numerictype([],31,30);
    Hm1c34.CustomAccumulatorDataType=numerictype([],33,30);
else
    Hm1c34.CustomProductDataType=numerictype([],31,29);
    Hm1c34.CustomAccumulatorDataType=numerictype([],31,29);
end

Hm2 = dec_int_func(2, hb2_coeff);
Hm2.FullPrecisionOverride = false;
Hm2.OutputDataType='Custom';
Hm2.CustomOutputDataType=numerictype([],16,14);
Hm2.CoefficientsDataType='Custom';
Hm2.CustomCoefficientsDataType=numerictype([],16);
Hm2.ProductDataType='Custom';
Hm2.CustomProductDataType=numerictype([],31,29);
Hm2.AccumulatorDataType = 'Custom';
if strcmp(input.RxTx, 'Rx')
    Hm2.CustomAccumulatorDataType=numerictype([],32,29);
else
    Hm2.CustomAccumulatorDataType=numerictype([],31,29);
end

Hm2c34 = dec_int_func(2, hb2_coeff);
Hm2c34.FullPrecisionOverride = false;
Hm2c34.OutputDataType='Custom';
Hm2c34.CustomOutputDataType=numerictype([],4,2);
Hm2c34.CoefficientsDataType='Custom';
Hm2c34.CustomCoefficientsDataType=numerictype([],16);
Hm2c34.ProductDataType='Custom';
Hm2c34.CustomProductDataType=numerictype([],31,29);
Hm2c34.AccumulatorDataType = 'Custom';
if strcmp(input.RxTx, 'Rx')
    Hm2c34.CustomAccumulatorDataType=numerictype([],32,29);
else
    Hm2c34.CustomAccumulatorDataType=numerictype([],31,29);
end

Hm3 = dec_int_func(2, hb3_coeff);
Hm3.FullPrecisionOverride = false;
Hm3.OutputDataType='Custom';
Hm3.CustomOutputDataType=numerictype([],8,6);
Hm3.CoefficientsDataType='Custom';
Hm3.CustomCoefficientsDataType=numerictype([],16);
Hm3.ProductDataType='Custom';
Hm3.AccumulatorDataType = 'Custom';
if strcmp(input.RxTx, 'Rx')
    Hm3.CustomProductDataType=numerictype([],19,18);
    Hm3.CustomAccumulatorDataType=numerictype([],21,18);
else
    Hm3.CustomProductDataType=numerictype([],19,17);
    Hm3.CustomAccumulatorDataType=numerictype([],19,17);
end

Hm4 = dec_int_func(3, dec_int3_coeff);
Hm4.FullPrecisionOverride = false;
Hm4.OutputDataType='Custom';
Hm4.CustomOutputDataType=numerictype([],16,14);
Hm4.CoefficientsDataType='Custom';
Hm4.CustomCoefficientsDataType=numerictype([],16);
Hm4.ProductDataType='Custom';
Hm4.CustomProductDataType=numerictype([],19,18);
Hm4.AccumulatorDataType = 'Custom';
if strcmp(input.RxTx, 'Rx')
    Hm4.CustomAccumulatorDataType=numerictype([],21,18);
else
    Hm4.CustomAccumulatorDataType=numerictype([],20,18);
end

hf.FullPrecisionOverride = false;
hf.OutputDataType='Custom';
hf.CustomOutputDataType=numerictype([],16,15);
hf.CoefficientsDataType='Custom';
hf.CustomCoefficientsDataType=numerictype([],16,15);
hf.ProductDataType='Custom';
hf.CustomProductDataType=numerictype([],16,15);
hf.AccumulatorDataType = 'Custom';
hf.CustomAccumulatorDataType=numerictype([],16,15);

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
switch enables
    case '1111' % only FIR
        filter_stages = {Hallpass};
    case '2111' % Hb1
        filter_stages = {Hm1};
    case '1211' % Hb2
        filter_stages = {Hm2};
    case '1121' % Hb3
        filter_stages = {Hm3};
    case '2211' % Hb2,Hb1
        filter_stages = {Hm2,Hm1};
    case '2121' % Hb3,Hb1
        filter_stages = {Hm3,Hm1c34};
    case '1221' % Hb3,Hb2
        filter_stages = {Hm3,Hm2c34};
    case '2221' % Hb3,Hb2,Hb1
        filter_stages = {Hm3,Hm2c34,Hm1};
    case '1113' % Dec/Int3
        filter_stages = {Hm4};
    case '2113' % Dec/Int3,Hb1
        filter_stages = {Hm4,Hm1c34};
    case '1213' % Dec/Int3,Hb2
        filter_stages = {Hm4,Hm2c34};
    case '2213' % Dec/Int3,Hb2,Hb1
        filter_stages = {Hm4,Hm2c34,Hm1};
    otherwise
        error('ddcresponse:IllegalOption', 'At least one of the stages must be there.')
end

% filter stages are reversed for Tx path
if strcmp(input.RxTx, 'Tx')
    filter_stages = fliplr(filter_stages);
end
dfilter = cascade(filter_stages{:});

Hmiddle = clone(dfilter);
if strcmp(input.RxTx, 'Rx')
    if strcmp(enables,'1111') ||  strcmp(enables,'2111') || strcmp(enables,'1211') || strcmp(enables,'1121') || strcmp(enables,'1113')
        addStage(Hmiddle,Hd1,1);
    else
        addStage(Hmiddle,Hd1,1);
    end
    addStage(Hmiddle,Hd2,1);
else
    if strcmp(enables,'1111') ||  strcmp(enables,'2111') || strcmp(enables,'1211') || strcmp(enables,'1121') || strcmp(enables,'1113')
        addStage(Hmiddle,Hd1);
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
    response = analogresp('Rx',w,input.converter_rate,b1,a1,b2,a2).*freqz(dfilter,w,input.converter_rate);
else
    response = freqz(dfilter,w,input.converter_rate).*analogresp('Tx',w,input.converter_rate,b1,a1,b2,a2);
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
    rg1 = analogresp('Rx',omega,input.converter_rate,b1,a1,b2,a2).*freqz(dfilter,omega,input.converter_rate);
else
    rg1 = freqz(dfilter,omega,input.converter_rate).*analogresp('Tx',omega,input.converter_rate,b1,a1,b2,a2);
end
phase = unwrap(angle(rg1));
gd1 = group_delay(omega,phase); % group delay on passband for Analog + Converter + HB
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
    wg1 = abs(analogresp('Rx',omega(Gpass+2:end),input.converter_rate,b1,a1,b2,a2).*freqz(dfilter,omega(Gpass+2:end),input.converter_rate));
    wg2 = (wg1)/(dBinv(-input.Astop));
else
    wg1 = abs(freqz(dfilter,omega(Gpass+2:end),input.converter_rate).*analogresp('Tx',omega(Gpass+2:end),input.converter_rate,b1,a1,b2,a2));
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
        N = min(16*floor(input.converter_rate/(input.Rdata)),128);
    else
        N = min(16*floor(input.converter_rate/(2*input.Rdata)),128);
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
    N = min(16*floor(input.converter_rate/(input.Rdata)),Nmax);
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
        Hdeq = design(d2,'equiripple','B1Weights',W3,'B2Weights',W4,'SystemObject',false);
        scoef = Hdeq.Numerator;
        for k = 1:length(scoef)
            scoef(k) = -scoef(k)*(-1)^(k-1);
        end
    else
        scoef = 0;
    end
    tap_store(i,1:M)=ccoef+scoef;
    
    Hmd = dec_int_func(input.FIR,tap_store(i,1:M));
    if ~isempty(ver('fixedpoint')) % Make sure fixed-point toolbox is installed
        if license('test','fixed_point_toolbox') % Try to checkout a license
            Hmd.Numerator = double(fi(Hmd.Numerator,true,16));
        end
    end
    if strcmp(input.RxTx, 'Rx')
        if strcmp(enables,'1111') || strcmp(enables,'2111') || strcmp(enables,'1211') || strcmp(enables,'1121') || strcmp(enables,'1113')
            addStage(dfilter,Hmd);
        else
            addStage(dfilter,Hmd);
        end
        rg_pass = abs(analogresp('Rx',omega(1:Gpass+1),input.converter_rate,b1,a1,b2,a2).*freqz(dfilter,omega(1:Gpass+1),input.converter_rate));
        rg_stop = abs(analogresp('Rx',omega(Gpass+2:end),input.converter_rate,b1,a1,b2,a2).*freqz(dfilter,omega(Gpass+2:end),input.converter_rate));
    else
        if strcmp(enables,'1111') || strcmp(enables,'2111') || strcmp(enables,'1211') || strcmp(enables,'1121') || strcmp(enables,'1113')
            addStage(dfilter, Hmd, 1);
        else
            addStage(dfilter, Hmd, 1);
        end
        rg_pass = abs(freqz(dfilter,omega(1:Gpass+1),input.converter_rate).*analogresp('Tx',omega(1:Gpass+1),input.converter_rate,b1,a1,b2,a2));
        rg_stop = abs(freqz(dfilter,omega(Gpass+2:end),input.converter_rate).*analogresp('Tx',omega(Gpass+2:end),input.converter_rate,b1,a1,b2,a2));
    end
    
    % quantitative values about actual passband and stopband
    Apass_actual_vector(i) = mag2db(max(rg_pass))-mag2db(min(rg_pass));
    Astop_actual_vector(i) = -mag2db(max(rg_stop));
    
    if input.int_FIR == 0
        h = tap_store(1,1:M);
        Apass_actual = Apass_actual_vector(1);
        Astop_actual = Astop_actual_vector(1);
        if strcmp(input.RxTx, 'Rx')
            removeStage(dfilter);
        else
            removeStage(dfilter, 1);
        end
        break
    elseif Apass_actual_vector(1) > input.Apass || Astop_actual_vector(1) < input.Astop
        h = tap_store(1,1:N);
        Apass_actual = Apass_actual_vector(1);
        Astop_actual = Astop_actual_vector(1);
        if strcmp(input.RxTx, 'Rx')
            removeStage(dfilter);
        else
            removeStage(dfilter, 1);
        end
        break
    elseif Apass_actual_vector(i) > input.Apass || Astop_actual_vector(i) < input.Astop
        h = tap_store(i-1,1:N+16);
        Apass_actual = Apass_actual_vector(i-1);
        Astop_actual = Astop_actual_vector(i-1);
        if strcmp(input.RxTx, 'Rx')
            removeStage(dfilter);
        else
            removeStage(dfilter, 1);
        end
        break
    else
        N = N-16;
        i = i+1;
        if strcmp(input.RxTx, 'Rx')
            removeStage(dfilter);
        else
            removeStage(dfilter, 1);
        end
    end
end

if strcmp(input.RxTx, 'Tx')
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
if ~isempty(ver('fixedpoint')) && license('test','fixed_point_toolbox') && license('checkout','fixed_point_toolbox')
    Hmd.Numerator = double(fi(Hmd.Numerator,true,16));
end
if strcmp(input.RxTx, 'Rx')
    addStage(dfilter, Hmd);
else
    addStage(dfilter, Hmd, 1);
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

bTFIR = min([16 - aTFIR,16]);
firtaps = Hmd.Numerator.*(2^bTFIR);
if ~isequal(double(int16(firtaps)),double(int32(firtaps)))
    firtaps = Hmd.Numerator.*(2^(bTFIR-1));
end

if length(firtaps) < 128
    firtaps = [firtaps,zeros(1,128-length(firtaps))];
end

output = input;

% externally accessible fields
output.firtaps = firtaps;
output.nfirtaps = length(h);
output.filter = dfilter;
output.gain = gain;
output.Hm1 = Hm1;
output.Hm2 = Hm2;
output.Hm3 = Hm3;
output.Hm4 = Hm4;
output.Hmd = Hmd;
output.enables = enables;
if isfield(input,'FPGAfilter')
    output.FPGAfilter = input.FPGAfilter;
else
    output.FPGAfilter = false;
end

% internal fields used by the GUI
output.Hanalog = Hanalog;
output.Apass_actual = Apass_actual;
output.Astop_actual = Astop_actual;
output.delay = delay;
output.grpdelayvar = grpdelayvar;
output.Hd1 = Hd1;
output.Hd2 = Hd2;
output.Hmiddle = Hmiddle;
output.a1 = a1;
output.b1 = b1;
output.a2 = a2;
output.b2 = b2;

function output = alias_b(f,fs)
output = f-fs*floor(f/fs+0.5);

% coerces the normalized cutoff frequency passed between 0.0 and 1.0
% for digital Butterworth filter designs
function Wn = coerce_cutoff(freq)
Wn = freq;
if Wn < 0.0
    Wn = 0.0 + eps;
elseif Wn > 1.0
    Wn = 1.0 - eps;
end

function dBoutput = dBinv(dBinput)
dBmin = -150;
if dBinput>dBmin
    dBoutput = 10^(dBinput/20);
else
    dBoutput = 0;
end

function t = group_delay(freq,phase)
% calculates the group delay from frequency data (in Hz) and phase data (in radians)

k = length(phase);

% unwrap phase data
phase = (180/pi)*unwrap(phase);

% calculate group delay
for n = 2:k-1
    t(n) = (-1/720) * (((phase(n) - phase(n - 1)) / (freq(n) - freq(n - 1)))+ ((phase(n + 1) - phase(n)) / (freq(n + 1) - freq(n))));
end
t(1) = (-1/360) * (((phase(2) - phase(1))/(freq(2) - freq(1))));
t(k) = (-1/360) * (((phase(k) - phase(k - 1))/(freq(k) - freq(k - 1))));
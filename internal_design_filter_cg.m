% internal_design_filter_cg
%
% This implementation of the ADI ad936x-filter-wizard supports
% code generation from MATLAB Coder
%
% Edits by: Travis F. Collins <travisfcollins@gmail.com>
%
% When calling this function utilize the function "process_input" to make
% sure all the necessary fields exist
%
% Todo:
% - Set fractionalLength based on FSR of inputs (ML has no doc on this)

%
% Inputs
% ============================================
% Rdata      = Input/output sample data rate (in Hz)
% Fpass      = Passband frequency (in Hz)
% Fstop      = Stopband frequency (in Hz)
% caldiv     = The actual discrete register value that describes the
%              rolloff for the analog filters
% FIR        = FIR interpolation/decimation factor
% HB1        = HB1 interpolation/decimation rates
% PLL_mult   = PLL multiplication
% Apass      = Max ripple allowed in passband (in dB)
% Astop      = Min attenuation in stopband (in dB)
% phEQ       = Phase equalization on (not -1)/off (-1)
% HB2        = HB2 interpolation/decimation rates
% HB3        = HB3 interpolation/decimation rates
% Type       = The type of filter required. one of:
%              'Lowpass' (default,and only support)
%              'Bandpass' (Not implemented)
%              'Equalize' (Not implemented)
%              'Root Raised Cosine' (Not implemented)
% RxTx       = Is this 'Rx' or 'Tx'?.
% RFbw
% DAC_div    = The ADC/DAC ratio, for Rx channels, this is
%              always '1', for Tx, it is either '1' or '2'
% converter_rate = Rate of converter
% PLL_rate   = the PLL rate in Hz
% Fcenter    = Center Frequency in Hz (only used for Bandpass),
%              otherwise 0
% wnom       = analog cutoff frequency (in Hz)
% FIRdBmin   = min rejection that FIR is required to have (in dB)
% int_FIR    = use AD9361 FIR on (1)/off (0)
%
% Outputs
% ===============================================
% firtaps          = fixed point FIR coefficients

%#codegen

function [outputTaps] = internal_design_filter_cg(...
    Rdata,...
    Fpass,...
    Fstop,...
    caldiv,...
    FIR,...
    HB1,...
    PLL_mult,...
    Apass,...
    Astop,...
    phEQ,...
    HB2,...
    HB3,...
    Type,...
    RxTx,...
    RFbw,...
    DAC_div,...
    converter_rate,...
    PLL_rate,...
    Fcenter,...
    wnom,...
    FIRdBmin,...
    int_FIR)
%%%%%%%%%%%%%%%%%%%% Perform Input Checks

% Scalar checks
doubleTypes = {Rdata,Fpass,Fstop,caldiv,FIR,HB1,PLL_mult,Apass,Astop,...
    phEQ,HB2,HB3,RFbw,DAC_div,converter_rate,PLL_rate,Fcenter,wnom,...
    FIRdBmin,int_FIR};
for t = 1:length(doubleTypes)
    assert(isa(doubleTypes{t},'double') && ...
        isreal(doubleTypes{t}) && ...
        all(size(doubleTypes{t}) == [1,1]));
end

% String checks
assert(isa(Type,'char') && isreal(Type) && all(size(Type) == [1,7])); % This is unused
assert(isa(RxTx,'char') && isreal(RxTx) && all(size(RxTx) == [1,2]));


%%%%%%%%%%%%%%%%%%%% Build processing struct

input = struct;
input.Rdata = Rdata;
input.Fpass = Fpass;
input.Fstop = Fstop;
input.caldiv = caldiv;
input.FIR = FIR;
input.HB1 = HB1;
input.PLL_mult = PLL_mult;
input.Apass = Apass;
input.Astop = Astop;
input.phEQ = phEQ;
input.HB2 = HB2;
input.HB3 = HB3;
input.Type = Type; % Not used
input.RxTx = RxTx;
input.RFbw = RFbw;
input.DAC_div = DAC_div;
input.converter_rate = converter_rate;
input.PLL_rate = PLL_rate;
input.Fcenter = Fcenter;
input.wnom = wnom;
input.FIRdBmin = FIRdBmin;
input.int_FIR = int_FIR;

%% Design analog filters
if strcmp(input.RxTx, 'Rx')
    wTIA = input.wnom*(2.5/1.4);

    % Define the analog filters (for design purpose)
    [b1,a1] = butter_cg(1,2*pi*wTIA,'s');  % 1st order
    [b2,a2] = butter_cg(3,2*pi*input.wnom,'s');    % 3rd order

    % Define the digital filters with fixed coefficients
    allpass_coeff = 1;
    hb1_coeff = 2^(-11)*[-8 0 42 0 -147 0 619 1013 619 0 -147 0 42 0 -8];
    hb2_coeff = 2^(-8)*[-9 0 73 128 73 0 -9];
    hb3_coeff = 2^(-4)*[1 4 6 4 1];
    dec_int3_coeff = 2^(-14)*[55 83 0 -393 -580 0 1914 4041 5120 4041 1914 0 -580 -393 0 83 55];

else
    wreal = input.wnom*(5.0/1.6);

    % Define the analog filters (for design purpose)
    [b1,a1] = butter_cg(3,2*pi*input.wnom,'s');     % 3rd order
    [b2,a2] = butter_cg(1,2*pi*wreal,'s');  % 1st order

    % Define the digital filters with fixed coefficients
    allpass_coeff = 1;
    hb1_coeff = 2^(-14)*[-53 0 313 0 -1155 0 4989 8192 4989 0 -1155 0 313 0 -53];
    hb2_coeff = 2^(-8)*[-9 0 73 128 73 0 -9];
    hb3_coeff = 2^(-2)*[1 2 1];
    dec_int3_coeff = (1/3)*2^(-13)*[36 -19 0 -156 -12 0 479 223 0 -1215 -993 0 3569 6277 8192 6277 3569 0 -993 -1215 0 223 479 0 -12 -156 0 -19 36];

end

%% Configure staging of filters
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
enables = [to_char(hb1),to_char(hb2),to_char(hb3),to_char(dec_int3)];

% Find out the best fit delay on passband
Nw = 2048;
w = zeros(1,Nw);
phi = zeros(1,Nw);

w(1) = -input.Fpass;
for i = 2:(Nw)
    w(i) = w(1)-2*w(1)*i/(Nw);
end

%% Generate target responses used in filter design phase

% Generate responses then convolve
combinedResponse = generateCascadedResponseRx(enables,w,input.converter_rate,...
    allpass_coeff,...
    hb1_coeff,...
    hb2_coeff,...
    hb3_coeff,...
    dec_int3_coeff,[]);

% Determine overall response with analog filters inline
assert( strcmp(input.RxTx, 'Rx') || strcmp(input.RxTx, 'Tx'), 'RxTx must be set to Rx or Tx');
response = combinedResponse.*analogresp(input.RxTx,w,input.converter_rate,b1,a1,b2,a2);

invariance = real(response).^2+imag(response).^2;
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
fg = zeros(1,Gpass+1);
omega = zeros(1,Gpass+1);

% passband
for i = 1:(Gpass+1)
    fg(i) = (i-1)/G;
    omega(i) = fg(i)*clkFIR;
end
% Generate responses then convolve
combinedResponse = generateCascadedResponseRx(enables,omega,input.converter_rate,...
    allpass_coeff,...
    hb1_coeff,...
    hb2_coeff,...
    hb3_coeff,...
    dec_int3_coeff,[]);

% Determine overall response with analog filters inline
assert( strcmp(input.RxTx, 'Rx') || strcmp(input.RxTx, 'Tx'), 'RxTx must be set to Rx or Tx');
rg1 = combinedResponse.*analogresp(input.RxTx,omega,input.converter_rate,b1,a1,b2,a2);

rg2 = exp(-1i*2*pi*omega*delay);
rg = rg2./rg1;
w = abs(rg1)/(dBinv(input.Apass/2)-1);

g = Gpass+1;

% Expand memory correctly
fg2 = zeros(1,length(Gstop:(G/2))+length(fg));
fg2(1:length(fg)) = fg;
omega2 = zeros(1,length(Gstop:(G/2))+length(omega));
omega2(1:length(omega)) = omega;
rgN = complex(zeros(1,length(Gstop:(G/2))+length(rg)));
rgN(1:length(rg)) = rg;

% stop band
for m = Gstop:(G/2)
    g = g+1;
    fg2(g) = m/G;
    omega2(g) = fg2(g)*clkFIR;
    rgN(g) = 0;
end

% Generate responses then convolve
combinedResponse = generateCascadedResponseRx(enables,omega2(Gpass+2:end),input.converter_rate,...
    allpass_coeff,...
    hb1_coeff,...
    hb2_coeff,...
    hb3_coeff,...
    dec_int3_coeff,[]);

assert( strcmp(input.RxTx, 'Rx') || strcmp(input.RxTx, 'Tx'), 'RxTx must be set to Rx or Tx');
wg1 = abs(combinedResponse.*analogresp(input.RxTx,omega2(Gpass+2:end),input.converter_rate,b1,a1,b2,a2));
if strcmp(input.RxTx, 'Rx')
    wg2 = (wg1)/(dBinv(-input.Astop));
else
    wg2 = (sqrt(input.FIR)*wg1)/(dBinv(-input.Astop));
end
wg3 = dBinv(input.FIRdBmin);
wg = max(wg2,wg3);
grid = fg2;
if input.phEQ == -1
    resp = abs(rgN);
else
    resp = rgN;
end
weight = [w wg];
weight = weight/max(weight);

%% Set up design for FIR filter
cr = real(resp);
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
        otherwise
            error('Wrong FIR Type');
    end
    N = min(16*floor(input.converter_rate*input.DAC_div/(2*input.Rdata)),Nmax);
end

tap_store = zeros(N/16,N);
Apass_actual_vector = zeros(N/16,1);
Astop_actual_vector = zeros(N/16,1);
i = 1;

%% Design filter
while (1)

    if input.int_FIR
        ccoef = firpm_cg(N-1, [F1(1),F1(end),F2(1),F2(end)], [A1,A2], [F1,F2], [W1,W2]);
    else
        % Check different designs until we reach required ripple condition
        R = db2mag(-input.Astop); % Peak Ripple
        ccoef = 0; % Predef type
        for k = 3:128
            [ccoef,valid,err] = firpm_cg(k, [F1(1),F1(end),F2(1),F2(end)], [A1,A2], [F1,F2], [W1,W2]);
            % Check if design meets specs
            if (err<R(1) && valid)
                break
            end
        end
    end
    M = length(ccoef);
    % Enable phase equalization and apply update to taps
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
            MN = N-1;
        else
            MN = M-1;
        end
        scoef = firpm_cg(MN, [F3(1),F3(end),F4(1),F4(end)], [A3,A4], [F3,F4], [W3,W4]);

        for k = 1:length(scoef)
            scoef(k) = -scoef(k)*(-1)^(k-1);
        end
    else
        scoef = zeros(size(ccoef));
    end
    tap_store(i,1:M)=ccoef+scoef; % scoef ==0 when no EQ

    tap_store(i,1:M) = determineBestFractionLength(tap_store,i,M);

    rg_pass = 0; %#ok<NASGU>
    rg_stop = 0; %#ok<NASGU>
    if strcmp(input.RxTx, 'Rx')
        combinedResponsePass = generateCascadedResponseRx(enables,omega2(1:Gpass+1),input.converter_rate,...
            allpass_coeff,...
            hb1_coeff,...
            hb2_coeff,...
            hb3_coeff,...
            dec_int3_coeff,tap_store(i,1:M));
        combinedResponseStop = generateCascadedResponseRx(enables,omega2(Gpass+2:end),input.converter_rate,...
            allpass_coeff,...
            hb1_coeff,...
            hb2_coeff,...
            hb3_coeff,...
            dec_int3_coeff,tap_store(i,1:M));

        rg_pass = abs(analogresp('Rx',omega2(1:Gpass+1),input.converter_rate,b1,a1,b2,a2).*combinedResponsePass);
        rg_stop = abs(analogresp('Rx',omega2(Gpass+2:end),input.converter_rate,b1,a1,b2,a2).*combinedResponseStop);
    else % TX
        combinedResponsePass = generateCascadedResponseRx(enables,omega2(1:Gpass+1),input.converter_rate,...
            allpass_coeff,...
            hb1_coeff,...
            hb2_coeff,...
            hb3_coeff,...
            dec_int3_coeff,tap_store(i,1:M));
        combinedResponseStop = generateCascadedResponseRx(enables,omega2(Gpass+2:end),input.converter_rate,...
            allpass_coeff,...
            hb1_coeff,...
            hb2_coeff,...
            hb3_coeff,...
            dec_int3_coeff,tap_store(i,1:M));

        rg_pass = abs(combinedResponsePass.*analogresp('Tx',omega2(1:Gpass+1),input.converter_rate,b1,a1,b2,a2));
        rg_stop = abs(combinedResponseStop.*analogresp('Tx',omega2(Gpass+2:end),input.converter_rate,b1,a1,b2,a2));

    end

    % quantitative values about actual passband and stopband
    Apass_actual_vector(i) = mag2db(max(rg_pass))-mag2db(min(rg_pass));
    Astop_actual_vector(i) = -mag2db(max(rg_stop));

    if input.int_FIR == 0
        h = tap_store(1,1:M);
        %Apass_actual = Apass_actual_vector(1);
        %Astop_actual = Astop_actual_vector(1);
        break
    elseif Apass_actual_vector(1) > input.Apass || Astop_actual_vector(1) < input.Astop
        h = tap_store(1,1:N);
        %Apass_actual = Apass_actual_vector(1);
        %Astop_actual = Astop_actual_vector(1);
        break
    elseif Apass_actual_vector(i) > input.Apass || Astop_actual_vector(i) < input.Astop
        h = tap_store(i-1,1:N+16);
        %Apass_actual = Apass_actual_vector(i-1);
        %Astop_actual = Astop_actual_vector(i-1);
        break
    else
        N = N-16;
        i = i+1;
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

% There will always be 128 taps output
numTaps = length(h);
firTapsPreScale = zeros(1,128);
firTapsPreScale(1:numTaps) = h;

%% Calculate group delay
% Hmd = dec_int_func(input.FIR,h(1:128));
%
% if ~isempty(ver('fixedpoint')) && license('test','fixed_point_toolbox') && license('checkout','fixed_point_toolbox')
%     Hmd.Numerator = double(fi(Hmd.Numerator,true,16));
% end
% if strcmp(input.RxTx, 'Rx')
%     addStage(dfilter, Hmd);
% else
%     addStage(dfilter, Hmd, 1);
% end

%gd2c = grpdelay(Hmd,omega1,clkFIR).*(1/clkFIR);
% gd2 = grpdelay_cg(firTapsPreScale,1,omega1,clkFIR).'.*(1/clkFIR);
%
% if input.phEQ == -1
%     groupdelay = gd1 + gd2;
% else
%     groupdelay = gd1 + gd2;
% end
% grpdelayvar = max(groupdelay)-min(groupdelay);

%% Determine Gains
aTFIR = 1 + ceil(log2(max(firTapsPreScale)));
% switch aTFIR
%     case 2
%         gain = 6;
%     case 1
%         gain = 0;
%     case 0
%         gain = -6;
%     otherwise
%         gain = -12;
% end
% 
% if strcmp(input.RxTx, 'Rx')
%     if aTFIR > 2
%         gain = 6;
%     end
% else
%     if input.FIR == 2
%         gain = gain+6;
%     elseif input.FIR == 4
%         gain = gain+12;
%     end
%     if gain > 0
%         gain = 0;
%     elseif gain < -6
%         gain = -6;
%     end
% end

%% Scale taps
bTFIR = 16 - aTFIR;
firtaps = int16(firTapsPreScale.*(2^bTFIR));

%output = input;

% %% Non-codegen outputs
% % externally accessible fields
% output.firtaps = firtaps;
% output.nfirtaps = length(h);
% %output.filter = dfilter;
% output.gain = gain;
% %output.Hm1 = Hm1;
% %output.Hm2 = Hm2;
% %output.Hm3 = Hm3;
% %output.Hm4 = Hm4;
% %output.Hmd = Hmd;
% output.enables = enables;
%
% % internal fields used by the GUI
% %output.Hanalog = Hanalog;
% output.Apass_actual = Apass_actual;
% output.Astop_actual = Astop_actual;
% %output.delay = delay;
% %output.grpdelayvar = grpdelayvar;
% %output.Hd1 = Hd1;
% %output.Hd2 = Hd2;
% %output.Hmiddle = Hmiddle;
% output.a1 = a1;
% output.b1 = b1;
% output.a2 = a2;
% output.b2 = b2;

%% For codegen only output taps
outputTaps = firtaps;


function output = alias_b(f,fs)
output = f-fs*floor(f/fs+0.5);

% coerces the normalized cutoff frequency passed between 0.0 and 1.0
% for digital Butterworth filter designs
function Wn = coerce_cutoff(freq) %#ok<DEFNU>
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

function t = group_delay(freq,phase) %#ok<DEFNU>
% calculates the group delay from frequency data (in Hz) and phase data (in radians)

k = length(phase);

% unwrap phase data
phase = (180/pi)*unwrap(phase);

t = zeros(1,k);

% calculate group delay
for n = 2:k-1
    t(n) = (-1/720) * (((phase(n) - phase(n - 1)) / (freq(n) - freq(n - 1)))+ ((phase(n + 1) - phase(n)) / (freq(n + 1) - freq(n))));
end
t(1) = (-1/360) * (((phase(2) - phase(1))/(freq(2) - freq(1))));
t(k) = (-1/360) * (((phase(k) - phase(k - 1))/(freq(k) - freq(k - 1))));


function d = to_char(c)
d = char(48+int8(c));

function taps = determineBestFractionLength(tap_store,i,M)
    % Codegen workaround for fixed fi call requirements
    signed = true; wordlength = 16;
    org = tap_store(i,1:M);
    e = zeros(16,1);
    k=1;
    r  = zeros(16,M);
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,1));
    e(1) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,2));
    e(2) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,3));
    e(3) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,4));
    e(4) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,5));
    e(5) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,6));
    e(6) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,7));
    e(7) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,8));
    e(8) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,9));
    e(9) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,10));
    e(10) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,11));
    e(11) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,12));
    e(12) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,13));
    e(13) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,14));
    e(14) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,15));
    e(15) = sum(abs(r(k,1:M)-org));k = k+1;
    r(k,1:M) = double(fi(tap_store(i,1:M),signed,wordlength,16));
    e(16) = sum(abs(r(k,1:M)-org));k = k+1;
    [~,fractionLength] = min(e);
    taps = r(fractionLength,1:M);



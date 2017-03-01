function [num, den, z, p] = butter_cg(n, Wn, varargin)
%BUTTER_CG Butterworth digital and analog filter design.  Codegen support
%
% This function is based on 'butter' by The MathWorks Inc.
%#codegen

% Set types we will only use
btype = 1; analog = 1;

if n>500
    error(message('signal:butter:InvalidRange'))
end
% Cast to enforce precision rules
Wn = double(Wn);

% step 1: get analog, pre-warped frequencies
if ~analog,
    fs = 2;
    u = 2*fs*tan(pi*Wn/fs);
else
    u = Wn;
end

Bw=[];
% step 2: convert to low-pass prototype estimate
if btype == 1	% lowpass
    Wn = u;
elseif btype == 2	% bandpass
    Bw = u(2) - u(1);
    Wn = sqrt(u(1)*u(2));	% center frequency
elseif btype == 3	% highpass
    Wn = u;
elseif btype == 4	% bandstop
    Bw = u(2) - u(1);
    Wn = sqrt(u(1)*u(2));	% center frequency
end

% Cast to enforce precision rules
validateattributes(n,{'numeric'},{'scalar','integer','positive'},'butter','N');
n = double(n);

% step 3: Get N-th order Butterworth analog lowpass prototype
[z,p,k] = buttap_cg(n);

% Transform to state-space
[a,b,c,d] = zp2ss_cg(z,p,k);

% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
if btype == 1		% Lowpass
    [a,b,c,d] = lp2lp_cg(a,b,c,d,Wn);

elseif btype == 2	% Bandpass
    [a,b,c,d] = lp2bp(a,b,c,d,Wn,Bw);

elseif btype == 3	% Highpass
    [a,b,c,d] = lp2hp(a,b,c,d,Wn);

elseif btype == 4	% Bandstop
    [a,b,c,d] = lp2bs(a,b,c,d,Wn,Bw);
end

% step 5: Use Bilinear transformation to find discrete equivalent:
if ~analog,
    [a,b,c,d] = bilinear(a,b,c,d,fs);
end

if nargout == 4
    num = a;
    den = b;
    z = c;
    p = d;
else	% nargout <= 3
    % Transform to zero-pole-gain and polynomial forms:
    if nargout == 3
        [z,p,k] = ss2zp(a,b,c,d,1); %#ok
        z = buttzeros(btype,n,Wn,analog);
        num = z;
        den = p;
        z = k;
    else % nargout <= 2
        den = poly(a);
        num = buttnum(btype,n,Wn,Bw,analog,den);
        % num = poly(a-b*c)+(d-1)*den;

    end
end

%---------------------------------
function b = buttnum(btype,n,Wn,Bw,analog,den)
% This internal function returns more exact numerator vectors
% for the num/den case.
% Wn input is two element band edge vector
if analog
    switch btype
        case 1  % lowpass
            b = [zeros(1,n) n^(-n)];
            b = real( b*polyval(den,-1i*0)/polyval(b,-1i*0) );
        case 2  % bandpass
            b = [zeros(1,n) Bw^n zeros(1,n)];
            b = real( b*polyval(den,-1i*Wn)/polyval(b,-1i*Wn) );
        case 3  % highpass
            b = [1 zeros(1,n)];
            b = real( b*den(1)/b(1) );
        case 4  % bandstop
            r = 1i*Wn*((-1).^(0:2*n-1)');
            b = poly(r);
            b = real( b*polyval(den,-1i*0)/polyval(b,-1i*0) );
    end
else
    Wn = 2*atan2(Wn,4);
    switch btype
        case 1  % lowpass
            r = -ones(n,1);
            w = 0;
        case 2  % bandpass
            r = [ones(n,1); -ones(n,1)];
            w = Wn;
        case 3  % highpass
            r = ones(n,1);
            w = pi;
        case 4  % bandstop
            r = exp(1i*Wn*( (-1).^(0:2*n-1)' ));
            w = 0;
    end
    b = poly(r);
    % now normalize so |H(w)| == 1:
    kern = exp(-1i*w*(0:length(b)-1));
    b = real(b*(kern*den(:))/(kern*b(:)));
end

function z = buttzeros(btype,n,Wn,analog)
% This internal function returns more exact zeros.
% Wn input is two element band edge vector
if analog
    % for lowpass and bandpass, don't include zeros at +Inf or -Inf
    switch btype
        case 1  % lowpass
            z = zeros(0,1);
        case 2  % bandpass
            z = zeros(n,1);
        case 3  % highpass
            z = zeros(n,1);
        case 4  % bandstop
            z = 1i*Wn*((-1).^(0:2*n-1)');
    end
else
    Wn = 2*atan2(Wn,4);
    switch btype
        case 1  % lowpass
            z = -ones(n,1);
        case 2  % bandpass
            z = [ones(n,1); -ones(n,1)];
        case 3  % highpass
            z = ones(n,1);
        case 4  % bandstop
            z = exp(1i*Wn*( (-1).^(0:2*n-1)' ));
    end
end

function [hh,w,s,options] = freqz_cg(b,a,w,Fs)
%FREQZ_CG Frequency response of digital filter with codegen support
%
% This function is based on 'freqz' by The MathWorks Inc.
%#codegen

isTF = true;

% Cast to enforce precision rules
options.nfft = double(w);
options.Fs = double(Fs);
options.w = double(w);
% Remaining are default or for advanced use
options.fvflag = double(1);
options.range = 'onesided';
options.centerdc = 0;
options.configlevel = 'omitted';


if isTF
    if length(a) == 1
        [h,w,options] = firfreqz(b/a,options);
    else
        [h,w,options] = iirfreqz(b,a,options);
    end
else
    validateattributes(b,{'double','single'},{},'freqz');
    [h,w,options] = sosfreqz(b,options);
end

% Generate the default structure to pass to freqzplot
s = struct;
s.plot = 'both';
s.fvflag = options.fvflag;
s.yunits = 'db';
s.xunits = 'rad/sample';
s.Fs     = options.Fs; % If rad/sample, Fs is empty
if ~isempty(options.Fs)
    s.xunits = 'Hz';
end

hh = h;
if isa(h,'single')
    % Cast to enforce precision rules. Cast when output is requested.
    % Otherwise, plot using double precision frequency vector.
    w = single(w);
end


%--------------------------------------------------------------------------
function [h,w,options] = firfreqz(b,options)

% Make b a row
b = b(:).';
n  = length(b);

w      = options.w;
Fs     = options.Fs;
nfft   = options.nfft;
fvflag = options.fvflag;

% Actual Frequency Response Computation
%if fvflag,
%   Frequency vector specified.  Use Horner's method of polynomial
%   evaluation at the frequency points and divide the numerator
%   by the denominator.
%
%   Note: we use positive i here because of the relationship
%            polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1))
%               ( assuming w = 2*pi*(0:length(a)-1)/length(a) )
%
if ~isempty(Fs) % Fs was specified, freq. vector is in Hz
    digw = 2.*pi.*w./Fs; % Convert from Hz to rad/sample for computational purposes
else
    digw = w;
end

s = exp(1i*digw); % Digital frequency must be used for this calculation
h = polyval(b,s)./exp(1i*digw*(n-1));

%--------------------------------------------------------------------------
function [h,w,options] = iirfreqz(b,a,options)
% Make b and a rows
b = b(:).';
a = a(:).';

nb = length(b);
na = length(a);
a  = [a zeros(1,nb-na)];  % Make a and b of the same length
b  = [b zeros(1,na-nb)];
n  = length(a); % This will be the new length of both num and den

w      = options.w;
Fs     = options.Fs;
nfft   = options.nfft;
fvflag = options.fvflag;

% Actual Frequency Response Computation
if fvflag,
    %   Frequency vector specified.  Use Horner's method of polynomial
    %   evaluation at the frequency points and divide the numerator
    %   by the denominator.
    %
    %   Note: we use positive i here because of the relationship
    %            polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1))
    %               ( assuming w = 2*pi*(0:length(a)-1)/length(a) )
    %
    if ~isempty(Fs), % Fs was specified, freq. vector is in Hz
        digw = 2.*pi.*w./Fs; % Convert from Hz to rad/sample for computational purposes
    else
        digw = w;
    end
    
    s = exp(1i*digw); % Digital frequency must be used for this calculation
    h = polyval(b,s) ./ polyval(a,s);
else
    % freqvector not specified, use nfft and RANGE in calculation
    s = find(strncmpi(options.range, {'twosided','onesided'}, length(options.range)));
    
    if s*nfft < n,
        % Data is larger than FFT points, wrap modulo s*nfft
        b = datawrap(b,s.*nfft);
        a = datawrap(a,s.*nfft);
    end
    
    % dividenowarn temporarily shuts off warnings to avoid "Divide by zero"
    h = dividenowarn(fft(b,s.*nfft),fft(a,s.*nfft)).';
    % When RANGE = 'half', we computed a 2*nfft point FFT, now we take half the result
    h = h(1:nfft);
    h = h(:); % Make it a column only when nfft is given (backwards comp.)
    w = freqz_freqvec(nfft, Fs, s);
    w = w(:); % Make it a column only when nfft is given (backwards comp.)
end
%--------------------------------------------------------------------------
function [h,w,options] = sosfreqz(sos,options)

[h, w] = iirfreqz(sos(1,1:3), sos(1,4:6), options);
for indx = 2:size(sos, 1)
    h = h.*iirfreqz(sos(indx,1:3), sos(indx,4:6), options);
end


% [EOF] freqz_cg.m

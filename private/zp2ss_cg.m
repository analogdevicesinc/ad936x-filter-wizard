function [a,b,c,d] = zp2ss(z,p,k)
%ZP2SS  Zero-pole to state-space conversion.
%   [A,B,C,D] = ZP2SS(Z,P,K)  calculates a state-space representation:
%       .
%       x = Ax + Bu
%       y = Cx + Du
%
%   for a system given a set of pole locations in column vector P,
%   a matrix Z with the zero locations in as many columns as there are
%   outputs, and the gains for each numerator transfer function in
%   vector K.  The A,B,C,D matrices are returned in block diagonal
%   form.  
%
%   The poles and zeros must correspond to a proper system. If the poles
%   or zeros are complex, they must appear in complex conjugate pairs,
%   i.e., the corresponding transfer function must be real.
%     
%   See also SS2ZP, ZP2TF, TF2ZP, TF2SS, SS2TF.

%   Thanks to G.F. Franklin
%   Copyright 1984-2008 The MathWorks, Inc.
[z,p,k,isSIMO] = parse_input(z,p,k);
     
if isSIMO
    % If it's multi-output, we can't use the nice algorithm
    % that follows, so use the numerically unreliable method
    % of going through polynomial form, and then return.
    [num,den] = zp2tf(z,p,k); % Suppress compile-time diagnostics
    [a,b,c,d] = tf2ss(num,den);
    return
end

% Strip infinities and throw away.
p = p(isfinite(p));
z = z(isfinite(z));

% Group into complex pairs
np = length(p);
nz = length(z);
try
    % z and p should have real elements and exact complex conjugate pair.
    z = cplxpair(z,0);
    p = cplxpair(p,0);
catch
    % If fail, revert to use the old default tolerance.
    % The use of tolerance in checking for real entries and conjugate pairs
    % may result in misinterpretation for edge cases. Please review the
    % process of how z and p are generated.
    z = cplxpair(z,1e6*nz*norm(z)*eps + eps);
    p = cplxpair(p,1e6*np*norm(p)*eps + eps);
end

% Initialize state-space matrices for running series
a=[]; b=zeros(0,1); c=ones(1,0); d=1;

% If odd number of poles AND zeros, convert the pole and zero
% at the end into state-space.
%   H(s) = (s-z1)/(s-p1) = (s + num(2)) / (s + den(2))
if rem(np,2) && rem(nz,2)
    a = p(np);
    b = 1;
    c = p(np) - z(nz);
    d = 1;
    np = np - 1;
    nz = nz - 1;
end

% If odd number of poles only, convert the pole at the
% end into state-space.
%  H(s) = 1/(s-p1) = 1/(s + den(2)) 
if rem(np,2)
    a = p(np);
    b = 1;
    c = 1;
    d = 0;
    np = np - 1;
end 

% If odd number of zeros only, convert the zero at the
% end, along with a pole-pair into state-space.
%   H(s) = (s+num(2))/(s^2+den(2)s+den(3)) 
if rem(nz,2)
    num = real(poly(z(nz)));
    den = real(poly(p(np-1:np)));
    wn = sqrt(prod(abs(p(np-1:np))));
    if wn == 0, wn = 1; end
    t = diag([1 1/wn]); % Balancing transformation
    a = t\[-den(2) -den(3); 1 0]*t;
    b = t\[1; 0];
    c = [1 num(2)]*t;
    d = 0;
    nz = nz - 1;
    np = np - 2;
end

% Now we have an even number of poles and zeros, although not 
% necessarily the same number - there may be more poles.
%   H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3))
% Loop through rest of pairs, connecting in series to build the model.
i = 1;
while i < nz
    index = i:i+1;
    num = real(poly(z(index)));
    den = real(poly(p(index)));
    wn = sqrt(prod(abs(p(index))));
    if wn == 0, wn = 1; end
    t = diag([1 1/wn]); % Balancing transformation
    a1 = t\[-den(2) -den(3); 1 0]*t;
    b1 = t\[1; 0];
    c1 = [num(2)-den(2) num(3)-den(3)]*t;
    d1 = 1;
    % [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1); 
    % Next lines perform series connection 
    ma1 = size(a,1);
    na2 = size(a1,2);
    a = [a zeros(ma1,na2); b1*c a1];
    b = [b; b1*d];
    c = [d1*c c1];
    d = d1*d;

    i = i + 2;
end

% Take care of any left over unmatched pole pairs.
%   H(s) = 1/(s^2+den(2)s+den(3))
while i < np
    den = real(poly(p(i:i+1)));
    wn = sqrt(prod(abs(p(i:i+1))));
    if wn == 0, wn = 1; end
    t = diag([1 1/wn]); % Balancing transformation
    a1 = t\[-den(2) -den(3); 1 0]*t;
    b1 = t\[1; 0];
    c1 = [0 1]*t;
    d1 = 0;
    % [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1);
    % Next lines perform series connection 
    ma1 = size(a,1);
    na2 = size(a1,2);
    a = [a zeros(ma1,na2); b1*c a1];
    b = [b; b1*d];
    c = [d1*c c1];
    d = d1*d;

    i = i + 2;
end

% Apply gain k:
c = c*k;
d = d*k;

%----------------------------------------------------------------------------
function [z,p,k,isSIMO] = parse_input(z,p,k)
%PARSE_INPUT   Make sure input args are valid.

% Initially assume it is a SISO system
isSIMO = 0;

% Check that p is a vector
if ~any(size(p)<2),
   ctrlMsgUtils.error('Controllib:general:pNotVector')
end
% Columnize p
p = p(:);

% Check that k is a vector
if ~any(size(k)<2),
   ctrlMsgUtils.error('Controllib:general:kNotVector')
end
% Columnize k
k = k(:);

% Check size of z
if any(size(z)<2),
   % z is a vector or an empty, columnize it
   z = z(:);
else
   % z is a matrix
   isSIMO = 1;
end

% Check for properness
if size(z,1) > length(p),
   % improper
   ctrlMsgUtils.error('Controllib:general:improperSystem')
end

% Check for the appropriate length of k
if length(k) ~= size(z,2) && (~isempty(z))
   ctrlMsgUtils.error('Controllib:general:zkLengthMismatch')
end

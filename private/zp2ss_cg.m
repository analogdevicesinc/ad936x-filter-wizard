function [a,b,c,d] = zp2ss_cg(z,p,k)
%ZP2SS  Zero-pole to state-space conversion. Codegen support
%[z,p,k,isSIMO] = parse_input(z,p,k);
isSIMO = 0;

% if isSIMO
%     % If it's multi-output, we can't use the nice algorithm
%     % that follows, so use the numerically unreliable method
%     % of going through polynomial form, and then return.
%     [num,den] = zp2tf(z,p,k); % Suppress compile-time diagnostics
%     [a,b,c,d] = tf2ss(num,den);
%     return
% end

% Strip infinities and throw away.
pF = p(isfinite(p));
zF = z(isfinite(z));

% Group into complex pairs
np = length(pF);
nz = length(zF);
% try
%     % z and p should have real elements and exact complex conjugate pair.
%     z = cplxpair(zF,0);
%     p = cplxpair(pF,0);
% catch
%     % If fail, revert to use the old default tolerance.
%     % The use of tolerance in checking for real entries and conjugate pairs
%     % may result in misinterpretation for edge cases. Please review the
%     % process of how z and p are generated.
%     z = cplxpair(zF,1e6*nz*norm(zF)*eps + eps);
%     p = cplxpair(pF,1e6*np*norm(pF)*eps + eps);
% end

% Initialize state-space matrices for running series
a=[]; b=zeros(0,1); c=ones(1,0); d=1;

% If odd number of poles AND zeros, convert the pole and zero
% at the end into state-space.
%   H(s) = (s-z1)/(s-p1) = (s + num(2)) / (s + den(2))
if rem(np,2) && rem(nz,2)
    a = pF(np);
    b = 1;
    c = pF(np) - zF(nz);
    d = 1;
    np = np - 1;
    nz = nz - 1;
end

% If odd number of poles only, convert the pole at the
% end into state-space.
%  H(s) = 1/(s-p1) = 1/(s + den(2)) 
if rem(np,2)
    a = pF(np);
    b = 1;
    c = 1;
    d = 0;
    np = np - 1;
end 

% If odd number of zeros only, convert the zero at the
% end, along with a pole-pair into state-space.
%   H(s) = (s+num(2))/(s^2+den(2)s+den(3)) 
if rem(nz,2)
    num = real(poly(zF(nz)));
    den = real(poly(pF(np-1:np)));
    wn = sqrt(prod(abs(pF(np-1:np))));
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
    num = real(poly(zF(index)));
    den = real(poly(pF(index)));
    wn = sqrt(prod(abs(pF(index))));
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
    den = real(poly(pF(i:i+1)));
    wn = sqrt(prod(abs(pF(i:i+1))));
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

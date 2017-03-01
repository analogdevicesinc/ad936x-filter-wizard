function [z,pV,k] = buttap_cg(n)
%BUTTAP Butterworth analog lowpass filter prototype. Codegen support
%
% This function is based on 'butterap' by The MathWorks Inc.
%#codegen

% Cast to enforce precision rules
nD = double(n);
% Poles are on the unit circle in the left-half plane.
z = [];
p = exp(1i*(pi*(1:2:nD-1)/(2*nD) + pi/2));
pR = [p; conj(p)];
pC = pR(:);
if rem(n,2)==1   % n is odd
    pV = [pC; -1];
else
    pV = pC;
end
k = real(prod(-pV));

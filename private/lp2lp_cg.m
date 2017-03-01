function [at,bt,ct,dt] = lp2lp_cg(a,b,c,d,wo)
%LP2LP Lowpass to lowpass analog filter transformation. Codegen support
%
% This function is based on 'lp2lp' by The MathWorks Inc.
%#codegen

% Transform lowpass to lowpass
at = wo*a;
bt = wo*b;
ct = c;
dt = d;

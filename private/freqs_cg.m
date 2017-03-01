function [h,ww] = freqs_cg(b,a,w)
%FREQS Laplace-transform (s-domain) frequency response with codegen support
%
% This function is based on 'freqs' by The MathWorks Inc.
%#codegen

narginchk(2,3);
nargoutchk(0,2);

if isempty(b)
  b = 0;
end
if isempty(a)
  a = 0;
end


validateattributes(b,{'numeric'},{'vector'},'freqs','B');
validateattributes(a,{'numeric'},{'vector'},'freqs','A');

[b,a] = removeTrailingZero(b,a);

if nargin == 2
  w = 200;
end

w = double(w);

s = 1i*w;
hh = polyval(b,s) ./ polyval(a,s);

if nargout == 0,
  % make sure W is never decreasing
  if w(1) > w(length(w))
    error(message('signal:freqs:IncreasingW','W','W'));
  end
  newplot;
  mag = abs(hh);   phase = angle(hh)*180/pi;
  subplot(211),loglog(w,mag),set(gca,'xgrid','on','ygrid','on')
  set(gca,'xlim',[w(1) w(length(w))])
  xlabel(getString(message('signal:freqs:Frequencyrads')))
  ylabel(getString(message('signal:freqs:Magnitude')))
  ax = gca;
  subplot(212), semilogx(w,phase),set(gca,'xgrid','on','ygrid','on')
  set(gca,'xlim',[w(1) w(length(w))])
  xlabel(getString(message('signal:freqs:Frequencyrads')))
  ylabel(getString(message('signal:freqs:Phasedegrees')))
  axes(ax)
elseif nargout == 1,
  h = hh;
elseif nargout == 2,
  h = hh;
  % Cast to enforce precision rules
  if isa(h,'single')
    ww = single(w);
  else
    ww = w;
  end
end
% end freqs


function [bR,aR] = removeTrailingZero(b,a)

b_len = numel(b);
a_len = numel(a);
b_lastnzidx = find(b,1,'last');
a_lastnzidx = find(a,1,'last');
trz_len = min(b_len-b_lastnzidx,a_len-a_lastnzidx);
if trz_len(1) > 0
    bR = b(1:b_len-trz_len(1));
    aR = a(1:a_len-trz_len(1));
else
    bR = b;
    aR = a;
end


% end removeTrailingZero

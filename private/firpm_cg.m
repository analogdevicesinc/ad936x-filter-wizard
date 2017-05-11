function [h,valid,err] = firpm_cg(order, ff, amplitudes, frequencies, weights, varargin)
%FIRPM Parks-McClellan optimal equiripple FIR filter design.
%
% This function is based on 'firpm' by The MathWorks Inc.
%#codegen

[grid,des,wt] = genWeights(order, ff, frequencies, weights, amplitudes);
% Workaround
nfilt = order + 1;
%ftype = 2;
%sign_val = 1;
hilbert = 0; % Always bandpass designs
neg = 0;

% cast to enforce precision rules
wt = double(wt);
des = double(des);
grid = double(grid);

% Call actual design algorithm
[h,err,valid] = remezm(nfilt,ff/2,grid/2,des,wt,neg);

err = abs(err);
h = h(:).';  % make it a row
h = [h sign(.5-(neg))*h(length(h)-rem(nfilt,2):-1:1)];
h = h(length(h):-1:1);
if neg && ~hilbert
    h = -h;  %make sure differentiator has correct sign
end


%%
function [grid,des,wt] = genWeights(filterOrder, bands, frequencies, weights, amplitudes)

    lgrid = 16;
    grid = firpmgrid_cg(filterOrder+1,lgrid,bands,0,0);

    orgFreqIndx = frequencies;

    positionsOfNewFreqIndx = zeros(size(grid));
    for ind = 1:length(positionsOfNewFreqIndx)

        [~,indx] = min( abs(orgFreqIndx-grid(ind)) );

        positionsOfNewFreqIndx(ind) = indx;
    end

    wt = weights(positionsOfNewFreqIndx);
    des = amplitudes(positionsOfNewFreqIndx);


%%
function [h,dev,valid]  = remezm(nfilt,edge,grid,des,wt,neg)
% remezm function
% Inputs
%     nfilt - filter length
%     edge - vector of band edges (between 0 and .5)
%     grid - frequency grid (between 0 and .5)
%     des - desired function on frequency grid
%     wt - weight function on frequency grid
%     neg == 1 ==> antisymmetric imp resp,
%         == 0 ==> symmetric imp resp
% Outputs
%     h - coefficients of basis functions
%     dev - computed error
%     iext - indices of extremal frequencies

valid = true;

nbands = length(edge)/2;
jb = 2*nbands;
nodd = rem(nfilt,2);      % nodd == 1 ==> filter length is odd
% nodd == 0 ==> filter length is even
nfcns = fix(nfilt/2);
if nodd == 1 && neg == 0
    nfcns = nfcns + 1;
end

ngrid = length(grid);

if neg <= 0
    if nodd ~= 1
        des = des./cos(pi*grid);
        wt = wt.*cos(pi*grid);
    end
elseif nodd ~= 1
    des = des./sin(pi*grid);
    wt = wt.*sin(pi*grid);
else
    des = des./sin(2*pi*grid);
    wt = wt.*sin(2*pi*grid);
end
temp = (ngrid-1)/nfcns;
j=1:nfcns;
iext = [fix([temp*(j-1)+1 ngrid])';0];
nz = nfcns + 1;

% Remez exchange loop
comp = -1;dtemp = -1;y1 = -1;luck = -1;nut1=-1;err=-1;y=-1;dev=-1;
itrmax = 250;
devl = -1;
nzz = nz + 1; x = zeros(1,nz);
niter = 0;
jchnge = 1;
jet = fix((nfcns - 1)/15) + 1;
ad = zeros(1,nz);


% index manager(s)
inextLen = nzz;

while jchnge > 0
    iext(nzz) = ngrid + 1;
    niter = niter + 1;
    if niter > itrmax
        break;
    end
    l = iext(1:nz)';
    x = cos(2*pi*grid(l));
    for nn = 1:nz
        ad(nn) = remezdd(nn, nz, jet, x);
    end
    add = ones(size(ad));
    add(2:2:nz) = -add(2:2:nz);
    dnum = ad*des(l)';
    dden = add*(ad./wt(l))';
    dev = dnum/dden;
    nu = 1;
    if dev > 0
        nu = -1;
    end
    dev = -nu*dev;
    y = des(l) + nu*dev*add./wt(l);
    if dev <= devl
        %warning(message('signal:firpm:DidNotConverge',niter))
        fprintf('%s\n','DidNotConverge');
        h = zeros(nfilt,1);
        dev = -1;
        %iext
        valid = false;
        return;
        %break;
    end
    devl = dev;
    jchnge = 0;
    k1 = iext(1);
    knz = iext(nz);
    klow = 0;
    nut = -nu;
    j = 1;
    flag34 = 1;
    while j < nzz
        kup = iext(j+1);
        l = iext(j) + 1;
        nut = -nut;
        if j == 2
            y1 = comp;
        end
        comp = dev;
        flag = 1;
        if l < kup
            % gee
            c = ad./(cos(2*pi*grid(l))-x);
            err = (c*y'/sum(c) - des(l))*wt(l);
            dtemp = nut*err - comp;
            if dtemp > 0
                comp = nut*err;
                l = l + 1;
                while l < kup
                    % gee
                    c = ad./(cos(2*pi*grid(l))-x);
                    err = (c*y'/sum(c) - des(l))*wt(l);
                    dtemp = nut*err - comp;
                    if dtemp > 0
                        comp = nut*err;
                        l = l + 1;
                    else
                        break;
                    end
                end
                iext(j) = l - 1;
                j = j + 1;
                klow = l - 1;
                jchnge = jchnge + 1;
                flag = 0;
            end
        end
        if flag
            l = l - 2;
            while l > klow
                % gee
                c = ad./(cos(2*pi*grid(l))-x);
                err = (c*y'/sum(c) - des(l))*wt(l);
                dtemp = nut*err - comp;
                if dtemp > 0 || jchnge > 0
                    break;
                end
                l = l - 1;
            end
            if l <= klow
                l = iext(j) + 1;
                if jchnge > 0
                    iext(j) = l - 1;
                    j = j + 1;
                    klow = l - 1;
                    jchnge = jchnge + 1;
                else
                    l = l + 1;
                    while l < kup
                        % gee
                        c = ad./(cos(2*pi*grid(l))-x);
                        err = (c*y'/sum(c) - des(l))*wt(l);
                        dtemp = nut*err - comp;
                        if dtemp > 0
                            break;
                        end
                        l = l + 1;
                    end
                    if l < kup && dtemp > 0
                        comp = nut*err;
                        l = l + 1;
                        while l < kup
                            % gee
                            c = ad./(cos(2*pi*grid(l))-x);
                            err = (c*y'/sum(c) - des(l))*wt(l);
                            dtemp = nut*err - comp;
                            if dtemp > 0
                                comp = nut*err;
                                l = l + 1;
                            else
                                break;
                            end
                        end
                        iext(j) = l - 1;
                        j = j + 1;
                        klow = l - 1;
                        jchnge = jchnge + 1;
                    else
                        klow = iext(j);
                        j = j + 1;
                    end
                end
            elseif dtemp > 0
                comp = nut*err;
                l = l - 1;
                while l > klow
                    % gee
                    c = ad./(cos(2*pi*grid(l))-x);
                    err = (c*y'/sum(c) - des(l))*wt(l);
                    dtemp = nut*err - comp;
                    if dtemp > 0
                        comp = nut*err;
                        l = l - 1;
                    else
                        break;
                    end
                end
                klow = iext(j);
                iext(j) = l + 1;
                j = j + 1;
                jchnge = jchnge + 1;
            else
                klow = iext(j);
                j = j + 1;
            end
        end
    end
    while j == nzz
        ynz = comp;
        k1 = min([k1 iext(1)]);
        knz = max([knz iext(nz)]);
        nut1 = nut;
        nut = -nu;
        l = 0;
        kup = k1;
        comp = ynz*1.00001;
        luck = 1;
        flag = 1;
        l = l + 1;
        while l < kup
            % gee
            c = ad./(cos(2*pi*grid(l))-x);
            err = (c*y'/sum(c) - des(l))*wt(l);
            dtemp = err*nut - comp;
            if dtemp > 0
                comp = nut*err;
                j = nzz;
                l = l + 1;
                while l < kup
                    % gee
                    c = ad./(cos(2*pi*grid(l))-x);
                    err = (c*y'/sum(c) - des(l))*wt(l);
                    dtemp = nut*err - comp;
                    if dtemp > 0
                        comp = nut*err;
                        l = l + 1;
                    else
                        break;
                    end
                end
                iext(j) = l - 1;
                j = j + 1;
                jchnge = jchnge + 1;
                flag = 0;
                break;
            end
            l = l + 1;
        end
        if flag
            luck = 6;
            l = ngrid + 1;
            klow = knz;
            nut = -nut1;
            comp = y1*1.00001;
            l = l - 1;
            while l > klow
                % gee
                c = ad./(cos(2*pi*grid(l))-x);
                err = (c*y'/sum(c) - des(l))*wt(l);
                dtemp = err*nut - comp;
                if dtemp > 0
                    j = nzz;
                    comp = nut*err;
                    luck = luck + 10;
                    l = l - 1;
                    while l > klow
                        % gee
                        c = ad./(cos(2*pi*grid(l))-x);
                        err = (c*y'/sum(c) - des(l))*wt(l);
                        dtemp = nut*err - comp;
                        if dtemp > 0
                            comp = nut*err;
                            l = l - 1;
                        else
                            break;
                        end
                    end
                    iext(j) = l + 1;
                    j = j + 1;
                    jchnge = jchnge + 1;
                    flag = 0;
                    break;
                end
                l = l - 1;
            end
            if flag
                flag34 = 0;
                if luck ~= 6
                    inextLen = length([k1 iext(2:nz-nfcns)' iext(nz-nfcns:nz-1)']'); % Update index
                    iext(1:inextLen) = [k1 iext(2:nz-nfcns)' iext(nz-nfcns:nz-1)']';
                    jchnge = jchnge + 1;
                end
                break;
            end
        end
    end
    if flag34 && j > nzz,
        if luck > 9
            inextLen = length([iext(2:nfcns+1)' iext(nfcns+1:nz-1)' iext(nzz) iext(nzz)]'); % Update index
            iext(1:inextLen) = [iext(2:nfcns+1)' iext(nfcns+1:nz-1)' iext(nzz) iext(nzz)]';
            jchnge = jchnge + 1;
        else
            y1 = max([y1 comp]);
            k1 = iext(nzz);
            l = ngrid + 1;
            klow = knz;
            nut = -nut1;
            comp = y1*1.00001;
            l = l - 1;
            while l > klow
                % gee
                c = ad./(cos(2*pi*grid(l))-x);
                err = (c*y'/sum(c) - des(l))*wt(l);
                dtemp = err*nut - comp;
                if dtemp > 0
                    j = nzz;
                    comp = nut*err;
                    luck = luck + 10;
                    l = l - 1;
                    while l > klow
                        % gee
                        c = ad./(cos(2*pi*grid(l))-x);
                        err = (c*y'/sum(c) - des(l))*wt(l);
                        dtemp = nut*err - comp;
                        if dtemp > 0
                            comp = nut*err;
                            l = l - 1;
                        else
                            break;
                        end
                    end
                    iext(j) = l + 1;
                    jchnge = jchnge + 1;
                    inextLen = length([iext(2:nfcns+1)' iext(nfcns+1:nz-1)' iext(nzz) iext(nzz)]'); % Update index
                    iext(1:inextLen) = [iext(2:nfcns+1)' iext(nfcns+1:nz-1)' iext(nzz) iext(nzz)]';
                    break;
                end
                l = l - 1;
            end
            if luck ~= 6
                inextLen = length([k1 iext(2:nz-nfcns)' iext(nz-nfcns:nz-1)']'); % Update index
                iext(1:inextLen) = [k1 iext(2:nz-nfcns)' iext(nz-nfcns:nz-1)']';
                jchnge = jchnge + 1;
            end
        end
    end
end

% Inverse Fourier transformation
bb = -1;aa = -1;% initialize memory
nm1 = nfcns - 1;
fsh = 1.0e-6;
%x(nzz) = -2;
x2 = [x,-2];
cn = 2*nfcns - 1;
delf = 1/cn;
l = 1;
kkk = 0;
if (edge(1) == 0 && edge(jb) == .5) || nfcns <= 3
    kkk = 1;
end
if kkk ~= 1
    dtemp = cos(2*pi*grid(1));
    dnum = cos(2*pi*grid(ngrid));
    aa = 2/(dtemp - dnum);
    bb = -(dtemp+dnum)/(dtemp - dnum);
end
a = zeros(1,nfcns);
for j = 1:nfcns
    ft = (j-1)*delf;
    xt = cos(2*pi*ft);
    if kkk ~= 1
        xt = (xt-bb)/aa;
        ft = acos(xt)/(2*pi);
    end
    xe = x2(l);
    while (xt <= xe) && (xe-xt >= fsh)
        l = l + 1;
        xe = x2(l);
    end
    if abs(xt-xe) < fsh
        a(j) = y(l);
    else
        grid(1) = ft;
        % gee
        c = ad./(cos(2*pi*ft)-x2(1:nz));
        a(j) = c*y'/sum(c);
    end
    l = max([1, l-1]);
end

dden = 2*pi/cn;

alpha = zeros(1,nfcns);
for j = 1:nfcns
    dnum = (j-1)*dden;
    if nm1 < 1
        alpha(j) = a(1);
    else
        alpha(j) = a(1) + 2*cos(dnum*(1:nm1))*a(2:nfcns)';
    end
end
alpha = [alpha(1) 2*alpha(2:nfcns)]'/cn;
p = zeros(1,nfcns);
q = zeros(1,nm1);
if kkk ~= 1
    p(1) = 2*alpha(nfcns)*bb + alpha(nm1);
    p(2) = 2*aa*alpha(nfcns);
    q(1) = alpha(nfcns-2) - alpha(nfcns);
    for j = 2:nm1
        if j == nm1
            aa = aa/2;
            bb = bb/2;
        end
        p(j+1) = 0;
        sel = 1:j;
        a(sel) = p(sel);
        p(sel) = 2*bb*a(sel);
        p(2) = p(2) + 2*a(1)*aa;
        jm1 = j - 1;
        sel = 1:jm1;
        p(sel) = p(sel) + q(sel) + aa*a(sel+1);
        jp1 = j + 1;
        sel = 3:jp1;
        p(sel) = p(sel) + aa*a(sel-1);
        if j ~= nm1
            sel = 1:j;
            q(sel) = -a(sel);
            q(1) = q(1) + alpha(nfcns - 1 - j);
        end
    end
    alpha(1:nfcns) = p(1:nfcns);
end

% alpha must be at lease >=3
if nfcns <= 3
    %alpha(nfcns + 1) = 0;
    %alpha(nfcns + 2) = 0;
    alphaFull = [alpha; 0; 0]';
else
    alphaFull = alpha';
end

%alpha=alpha';
% now that's done!

if neg <= 0
    if nodd ~= 0
        h = [.5*alphaFull(nz-1:-1:nz-nm1) alphaFull(1)];
    else
        h = .25*[alphaFull(nfcns) alphaFull(nz-2:-1:nz-nm1)+alphaFull(nfcns:-1:nfcns-nm1+2) ...
            2*alphaFull(1)+alphaFull(2)];
    end
elseif nodd ~= 0
    h = .25*[alphaFull(nfcns) alphaFull(nm1)];
    h = [h .25*[alphaFull(nz-3:-1:nz-nm1)-alphaFull(nfcns:-1:nfcns-nm1+3) ...
        2*alphaFull(1)-alphaFull(3) ] 0];
else
    h = .25*[alphaFull(nfcns) alphaFull(nz-2:-1:nz-nm1)-alphaFull(nfcns:-1:nfcns-nm1+2) ...
        2*alphaFull(1)-alphaFull(2)];
end

%%
function y = remezdd(k, n, m, x)
%REMEZDD Lagrange interpolation coefficients.

%   Author: T. Krauss 1993
%       Was Revision: 1.4, Date: 1994/01/25 17:59:44

y = 1;
q = x(k);
for l=1:m
    xx = 2*(q - x(l:m:n));
    y = y*prod(xx(xx ~= 0));
end
y=1/y;
% EOF

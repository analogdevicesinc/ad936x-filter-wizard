function combinedResponse = generateCascadedResponseRx(enables,w,Fs,...
    allpass_coeff,...
    hb1_coeff,...
    hb2_coeff,...
    hb3_coeff,...
    dec_int3_coeff,extraTaps)
%#codegen

wd = double(w); % Cast
stages = 0; %#ok<NASGU>

switch enables
    case '1111' % only FIR
        combinedResponse = freqz_cg(us(allpass_coeff,1),1,wd,Fs);
        stages = 1;
        
    case '2111' % Hb1
        combinedResponse = freqz_cg(us(hb1_coeff,1),1,wd,Fs);
        stages = 1;
        
    case '1211' % Hb2
        combinedResponse = freqz_cg(us(hb2_coeff,1),1,wd,Fs);
        stages = 1;
        
    case '1121' % Hb3
        combinedResponse = freqz_cg(us(hb3_coeff,1),1,wd,Fs);
        stages = 1;
        
    case '2211' % Hb2,Hb1
        d1 = freqz_cg(us(hb1_coeff,2),1,wd,Fs);
        d2 = freqz_cg(us(hb2_coeff,1),1,wd,Fs);
        combinedResponse = d1.*d2;
        stages = 2;
        
    case '2121' % Hb3,Hb1
        d1 = freqz_cg(us(hb1_coeff,2),1,wd,Fs);
        d2 = freqz_cg(us(hb3_coeff,1),1,wd,Fs);
        combinedResponse = d1.*d2;
        stages = 2;
        
    case '1221' % Hb3,Hb2
        d1 = freqz_cg(us(hb2_coeff,2),1,wd,Fs);
        d2 = freqz_cg(us(hb3_coeff,1),1,wd,Fs);
        combinedResponse = d1.*d2;
        stages = 2;
        
    case '2221' % Hb3,Hb2,Hb1
        d1 = freqz_cg(us(hb1_coeff,4),1,wd,Fs);
        d2 = freqz_cg(us(hb2_coeff,2),1,wd,Fs);
        d3 = freqz_cg(us(hb3_coeff,1),1,wd,Fs);
        combinedResponse = d1.*d2.*d3;
        stages = 3;
        
    case '1113' % Dec/Int3
        combinedResponse = freqz_cg(us(dec_int3_coeff,1),1,wd,Fs);
        stages = 1;
        
    case '2113' % Dec/Int3,Hb1
        d1 = freqz_cg(us(hb1_coeff,2),1,wd,Fs);
        d2 = freqz_cg(us(dec_int3_coeff,1),1,wd,Fs);
        combinedResponse = d1.*d2;
        stages = 2;
        
    case '1213' % Dec/Int3,Hb2
        d1 = freqz_cg(us(hb2_coeff,2),1,wd,Fs);
        d2 = freqz_cg(us(dec_int3_coeff,1),1,wd,Fs);
        combinedResponse = d1.*d2;
        stages = 2;
        
    case '2213' % Dec/Int3,Hb2,Hb1
        d1 = freqz_cg(us(hb1_coeff,4),1,wd,Fs);
        d2 = freqz_cg(us(hb2_coeff,2),1,wd,Fs);
        d3 = freqz_cg(us(dec_int3_coeff,1),1,wd,Fs);
        combinedResponse = d1.*d2.*d3;
        stages = 3;
        
    otherwise
        error('At least one set of stages must be selected.')
end

% Add filter extra filter to end of cascade
if ~isempty(extraTaps)
    dn = freqz_cg(us(extraTaps,2^stages),1,wd,Fs);
    combinedResponse = combinedResponse.*dn;
end


end

function u = us(o,n)
u = upsample(o,n);
u = u(1:end-n+1);
end


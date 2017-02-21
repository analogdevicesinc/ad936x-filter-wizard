function combinedResponse = generateCascadedResponseRx(enables,w,Fs,...
    allpass_coeff,...
    hb1_coeff,...
    hb2_coeff,...
    hb3_coeff,...
    dec_int3_coeff,extraTaps,frontFilt)
%#codegen

wd = double(w); % Cast
stages = 0; %#ok<NASGU>

if ~isempty(frontFilt) && ~isempty(extraTaps)
    error('Cannot add to both front and back');
end

% Add filter extra filter to front of cascade
% if ~isempty(frontFilt)
%     scaler = freqz_cg(frontFilt,1,wd,Fs);
%     upSamp = 2;
% else
%     scaler = 1;
%     upSamp = 1;
% end
scaler = 1;
upSamp = 1;

switch enables
    case '1111' % only FIR
        combinedResponse = scaler.*freqz_cg(us(allpass_coeff,upSamp),1,wd,Fs);
        stages = 1;
        
    case '2111' % Hb1
        combinedResponse = scaler.*freqz_cg(us(hb1_coeff,upSamp),1,wd,Fs);
        stages = 1;
        
    case '1211' % Hb2
        combinedResponse = scaler.*freqz_cg(us(hb2_coeff,upSamp),1,wd,Fs);
        stages = 1;
        
    case '1121' % Hb3
        combinedResponse = scaler.*freqz_cg(us(hb3_coeff,upSamp),1,wd,Fs);
        stages = 1;
        
    case '2211' % Hb2,Hb1
        d1 = freqz_cg(us(hb1_coeff,upSamp*2),1,wd,Fs);
        d2 = freqz_cg(us(hb2_coeff,upSamp*1),1,wd,Fs);
        combinedResponse = scaler.*d1.*d2;
        stages = 2;
        
    case '2121' % Hb3,Hb1
        d1 = freqz_cg(us(hb1_coeff,upSamp*2),1,wd,Fs);
        d2 = freqz_cg(us(hb2_coeff,upSamp*1),1,wd,Fs);
        combinedResponse = scaler.*d1.*d2;
        stages = 2;
        
    case '1221' % Hb3,Hb2
        d1 = freqz_cg(us(hb1_coeff,upSamp*2),1,wd,Fs);
        d2 = freqz_cg(us(hb3_coeff,upSamp*1),1,wd,Fs);
        combinedResponse = scaler.*d1.*d2;
        stages = 2;
        
    case '2221' % Hb3,Hb2,Hb1
        d1 = freqz_cg(us(hb1_coeff,upSamp*4),1,wd,Fs);
        d2 = freqz_cg(us(hb2_coeff,upSamp*2),1,wd,Fs);
        d3 = freqz_cg(us(hb3_coeff,upSamp*1),1,wd,Fs);
        combinedResponse = scaler.*d1.*d2.*d3;
        stages = 3;
        
    case '1113' % Dec/Int3
        combinedResponse = scaler.*freqz_cg(us(dec_int3_coeff,upSamp),1,wd,Fs);
        stages = 1;
        
    case '2113' % Dec/Int3,Hb1
        d1 = freqz_cg(us(hb1_coeff,upSamp*2),1,wd,Fs);
        d2 = freqz_cg(us(dec_int3_coeff,upSamp*1),1,wd,Fs);
        combinedResponse = scaler.*d1.*d2;
        stages = 2;
        
    case '1213' % Dec/Int3,Hb2
        d1 = freqz_cg(us(hb2_coeff,upSamp*2),1,wd,Fs);
        d2 = freqz_cg(us(dec_int3_coeff,upSamp*1),1,wd,Fs);
        combinedResponse = scaler.*d1.*d2;
        stages = 2;
        
    case '2213' % Dec/Int3,Hb2,Hb1
        d1 = freqz_cg(us(hb1_coeff,upSamp*4),1,wd,Fs);
        d2 = freqz_cg(us(hb2_coeff,upSamp*2),1,wd,Fs);
        d3 = freqz_cg(us(dec_int3_coeff,upSamp*1),1,wd,Fs);
        combinedResponse = scaler.*d1.*d2.*d3;
        stages = 3;
        
    otherwise
        error('ddcresponse:IllegalOption', 'At least one of the stages must be there.')
end

% Add filter extra filter to end of cascade
if ~isempty(extraTaps)
    dn = freqz_cg(us(extraTaps,2^stages),1,wd,Fs);
    combinedResponse = combinedResponse.*dn;
end


end

function u = us(o,n)
u = upsample(o,n);
end


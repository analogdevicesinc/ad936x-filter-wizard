function abc = analogresp(type,f,Fconverter,b1,a1,b2,a2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch type
    case 'Tx'
        abc = sinc(f/Fconverter).*freqs_cg(b1,a1,2*pi*f).*freqs_cg(b2,a2,2*pi*f);
    case 'Rx'
        abc = freqs_cg(b1,a1,2*pi*f).*freqs_cg(b2,a2,2*pi*f).*(sinc(f/Fconverter).^3);
    otherwise % Default to Rx
        abc = freqs_cg(b1,a1,2*pi*f).*freqs_cg(b2,a2,2*pi*f).*(sinc(f/Fconverter).^3);
end
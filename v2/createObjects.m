function [FilterStagesTX, FilterStagesRX] = createObjects(designer)

%% RX
% Digital representation of the analog filters (It is an approximation for group delay calculation only)
wTIA = designer.RxAnalogCuttoffFrequency*(2.5/1.4);
[z1,p1,k1] = butter(3,coerce_cutoff(designer.RxAnalogCuttoffFrequency/(designer.ADCRate/2)),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
Hd1=dsp.BiquadFilter('SOSMatrix',sos1,'ScaleValues',g1);
[z2,p2,k2] = butter(1,coerce_cutoff(wTIA/(designer.ADCRate/2)),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
Hd2=dsp.BiquadFilter('SOSMatrix',sos2,'ScaleValues',g2);
analogRX = cascade(Hd2,Hd1);

filterStagesRX = cell(1,3);
Configs = {designer.HB1ConfigRx,designer.HB2ConfigRx,designer.HB3ConfigRx,designer.HB3ConfigRxDec3};
HBSetup = [designer.RxHB1,designer.RxHB2,designer.RxHB3];
Coeffs = {designer.RxHB1Coeff,designer.RxHB2Coeff,designer.RxHB3Coeff,designer.RxHB3CoeffDec3};
for stage = 1:3
    if HBSetup(stage) == 2
        PV = designer.getPVPair(Configs{stage});
        filterStagesRX{stage} = dsp.FIRDecimator(2,PV{:},'Numerator',Coeffs{stage});
    elseif HBSetup(stage) == 3 % Only HB3 can be setup with Dec3
        PV = designer.getPVPair(Configs{4});
        filterStagesRX{stage} = dsp.FIRDecimator(3,PV{:},'Numerator',Coeffs{4});
    else
        filterStagesRX{stage} = dsp.FIRDecimator(1,designer.AllPassCoeff);
    end
end
%% TX
wreal = designer.TxAnalogCuttoffFrequency*(5.0/1.6);
% Digital representation of the analog filters (It is an approximation for group delay calculation only)
[z1,p1,k1] = butter(3,coerce_cutoff(designer.TxAnalogCuttoffFrequency/(designer.ADCRate/2)),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
Hd1=dsp.BiquadFilter('SOSMatrix',sos1,'ScaleValues',g1);
[z2,p2,k2] = butter(1,coerce_cutoff(wreal/(designer.ADCRate/2)),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
Hd2=dsp.BiquadFilter('SOSMatrix',sos2,'ScaleValues',g2);
analogTX = cascade(Hd1,Hd2);

filterStagesTX = cell(1,3);
Configs = {designer.HB1ConfigTx,designer.HB2ConfigTx,designer.HB3ConfigTx,designer.HB3ConfigTxInt3};
HBSetup = [designer.TxHB1,designer.TxHB2,designer.TxHB3];
Coeffs = {designer.TxHB1Coeff,designer.TxHB2Coeff,designer.TxHB3Coeff,designer.TxHB3CoeffInt3};
for stage = 1:3
    if HBSetup(stage) == 2
        PV = designer.getPVPair(Configs{stage});
        filterStagesTX{stage} = dsp.FIRInterpolator(2,PV{:},'Numerator',Coeffs{stage});
    elseif HBSetup(stage) == 3 % Only HB3 can be setup with Int3
        PV = designer.getPVPair(Configs{4});
        filterStagesTX{stage} = dsp.FIRInterpolator(3,PV{:},'Numerator',Coeffs{4});
    else
        filterStagesTX{stage} = dsp.FIRInterpolator(1,designer.AllPassCoeff);
    end
end
%% Save into cascades
filterStagesRX = fliplr(filterStagesRX);

s = [{analogRX},filterStagesRX(:)'];
FilterStagesRX = cascade(s{:});
s = [filterStagesTX(:)',{analogTX}];
FilterStagesTX = cascade(s{:});

end

% coerces the normalized cutoff frequency passed between 0.0 and 1.0
% for digital Butterworth filter designs
function Wn = coerce_cutoff(freq)
Wn = freq;
if Wn < 0.0
    Wn = 0.0 + eps;
elseif Wn > 1.0
    Wn = 1.0 - eps;
end
end

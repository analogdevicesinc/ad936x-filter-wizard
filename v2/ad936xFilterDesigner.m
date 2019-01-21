classdef ad936xFilterDesigner < ad936x
    %ad936xFilterDesigner Design specific filters for AD936X transceivers
    %   The designer creates a design constrain based on the configuration
    %   of the transceiver and additional objectives defined in this class.
    %
    %   FD = ad936xFilterDesigner will create an instance of the 
    %   ad936xFilterDesigner class. Filters are created using the
    %   'designFilter' method, which requires a single input argument to
    %   design either a transmit ('Tx') or receiver ('Rx') filter.
    %   Four output arguments are returned by the 'designFilter' method:
    %   [found, taps, numTaps, preScaledTaps] = fd.designFilter('Tx');
    %       found: Logical which denotes if the design object was achieved
    %       taps: int16 vector of length 128 which contains the designed
    %       filter
    %       numTaps: Integer denoting number of taps to be used of the 128
    %       preScaleTaps: double vector of designed taps before converted
    %       to int16 fullscale
    %
    %   The 'designFilter' will always return the shortest length filter
    %   which meets the design criteria. If design criteria is not met,
    %   the longest possible filter is returned which was designed that can
    %   be correctly loaded into the transceiver.
    %
    %   % Example usage:
    %
    %   fd = ad936xFilterDesigner;
    %   fd.DataRate = 1e6;
    %   fd.AutoSetRates();
    %   [found,taps] = fd.designFilter('Tx');

    properties
        %TxApass Tx Amplitude Passband
        %   Passband max ripple allowed in passband in dB for transmitter 
        %   path
        TxApass = 0.125;
        %TxAstop Tx Attenuation Stopband
        %   Stopband min attenuation in dB for transmitter path
        TxAstop = 85;
        %TxFpass Tx Start Frequency Transistion
        %	Start frequency in Hz for transistion band of designed FIR for
        %	transmitter path
        TxFpass = 2560000;
        %TxFstop Tx Stop Frequency Transistion
        %	Stop frequency in Hz for transistion band of designed FIR for 
        %   transmitter path
        TxFstop = 3200000;
        %TxFIRdBmin Tx FIR dB Minimum
        %   Minimum required rejection for FIR in dB for transmitter path
        TxFIRdBmin = 0;
        %TxUseFIR Tx Use FIR
        %   Use FIR in transmitter path. Otherwise designer will implement
        %   an unconstrained filter using up to 128 taps.
        TxUseFIR = true;
        %TxPhaseEQ Tx Phase Equalize
        %   Equalize phase response in transmitter path with generated FIR
        TxPhaseEQ = false;
        %TxTargetDelay Tx Target Delay
        %   Target phase delay of FIR in nanoseconds for transmitter path.
        %   This parameter is only used when TxPhaseEQ is true.
        TxTargetDelay = 0;
        %RxApass Rx Amplitude Passband
        %   Passband max ripple allowed in passband in dB for receive 
        %   path
        RxApass = 0.125;
        %RxAstop Rx Attenuation Stopband
        %   Stopband min attenuation in dB for receive path
        RxAstop = 85;
        %RxFpass Rx Start Frequency Transistion
        %	Start frequency in Hz for transistion band of designed FIR for
        %	receive path
        RxFpass = 2560000;
        %RxFstop Rx Stop Frequency Transistion
        %	Stop frequency in Hz for transistion band of designed FIR for 
        %   receive path
        RxFstop = 3200000;
        %RxFIRdBmin Rx FIR dB Minimum
        %   Minimum required rejection for FIR in dB for receive path
        RxFIRdBmin = 0;
        %RxUseFIR Rx Use FIR
        %   Use FIR in receive path. Otherwise designer will implement
        %   an unconstrained filter using up to 128 taps.
        RxUseFIR = true;
        %RxPhaseEQ Rx Phase Equalize
        %   Equalize phase response in receive path with generated FIR
        RxPhaseEQ = false;
        %RxTargetDelay Rx Target Delay
        %   Target phase delay of FIR in nanoseconds for receive path.
        %   This parameter is only used when TxPhaseEQ is true.
        RxTargetDelay = 0;    
    end
    
    properties (Hidden)
        F1
        F2
        W1
        W2
        A1
        A2
        
        grid
        weight
        omega2
        resp
        
        Gpass
        Gstop
    end
    
    properties (Constant, Hidden)
       G = 16384; 
    end
    
    methods
        function obj = FilterDesignerCurrent()
        end
        
        function [found, h, numTaps, firTapsPreScale,stats] = designFilter(obj,type)
            if obj.ValidConfiguration
                availableTaps = 16:16:obj.AvailableFIRTaps;
                delay = obj.createObjective(type);
                [found, h, numTaps, firTapsPreScale,stats] = obj.designer(type, availableTaps);
                % Collect stats from design
                stats.DelayVarianceNanoSeconds = delay;
                stats.UsedTaps = numTaps;
            else
                found = false; %#ok<*NASGU>
                h = int16([]);
                numTaps = 0;
                firTapsPreScale = [];
                stats = [];
                e = getConfigurationError(obj);
                error(['Chip configuration invalid ',e.msg]);
            end
        end
    end
    
    methods (Access = protected, Hidden)
        function [found, firtaps, numTaps, firTapsPreScale, stats] = designer(obj,type,availableTaps)
            
            found = false;
            
            if strcmpi(type,'RX')
                Astop = obj.RxAstop;
                Apass = obj.RxApass;
                FIR = obj.RxFIR;
                phEQ = obj.RxPhaseEQ;
                UseFIR = obj.RxUseFIR;
            else
                Astop = obj.TxAstop;
                Apass = obj.TxApass;
                FIR = obj.TxFIR;
                phEQ = obj.TxPhaseEQ;
                UseFIR = obj.TxUseFIR;
            end
            
            %% Design filter
            taps = zeros(1,128);
            for Nindx = 1:length(availableTaps)
                N = availableTaps(Nindx);
                if UseFIR
                    ccoef = firpm_cg(N-1, [obj.F1(1),obj.F1(end),obj.F2(1),obj.F2(end)], [obj.A1,obj.A2], [obj.F1,obj.F2], [obj.W1,obj.W2]);
                else
                    % Design arbitrary filter between lengths [3 128]
                    % Check different designs until we reach required ripple condition
                    R = db2mag(-Astop); % Peak Ripple
                    ccoef = 0; % Predef type
                    for k = 3:128
                        [ccoef,valid,err] = firpm_cg(k, [obj.F1(1),obj.F1(end),obj.F2(1),obj.F2(end)], [obj.A1,obj.A2], [obj.F1,obj.F2], [obj.W1,obj.W2]);
                        % Check if design meets specs
                        if (err<R(1) && valid)
                            break
                        end
                    end
                end
                M = length(ccoef);
                % Enable phase equalization and apply update to taps
                if phEQ ~= -1
                    sg = 0.5-obj.grid(end:-1:1);
                    sr = imag(obj.resp(end:-1:1));
                    sw = obj.weight(end:-1:1);
                    F3 = sg(1:obj.G/2-obj.Gstop+1)*2;
                    F4 = sg(obj.G/2-obj.Gstop+2:end)*2;
                    A3 = sr(1:obj.G/2-obj.Gstop+1);
                    A4 = sr(obj.G/2-obj.Gstop+2:end);
                    W3 = sw(1:obj.G/2-obj.Gstop+1);
                    W4 = sw(obj.G/2-obj.Gstop+2:end);
                    if UseFIR
                        MN = N-1;
                    else
                        MN = M-1;
                    end
                    scoef = firpm_cg(MN, [F3(1),F3(end),F4(1),F4(end)], [A3,A4], [F3,F4], [W3,W4]);
                    
                    for k = 1:length(scoef)
                        scoef(k) = -scoef(k)*(-1)^(k-1);
                    end
                else
                    scoef = zeros(size(ccoef));
                end
                taps = ccoef+scoef; % scoef ==0 when no EQ
                
                taps  = obj.determineBestFractionLength(taps, M);
                
                rg_pass = 0;
                rg_stop = 0;
                rg_pass = abs(generateCombinedResponse(obj, obj.omega2(1:obj.Gpass+1), type, true, taps));
                rg_stop = abs(generateCombinedResponse(obj, obj.omega2(obj.Gpass+2:end), type, true, taps));
                
                % quantitative values about actual passband and stopband
                Apass_actual = mag2db(max(rg_pass))-mag2db(min(rg_pass));
                Astop_actual = -mag2db(max(rg_stop));
                
                if UseFIR == 0
                    break
                elseif Apass_actual <= Apass && Astop_actual >= Astop
                    found = true;
                    break
                end
            end
            
            h = taps;
            
            if strcmpi(type, 'TX')
                if UseFIR && FIR == 2
                    R = rem(length(h),32);
                    if R ~= 0
                        h = [zeros(1,8),h,zeros(1,8)];
                    end
                elseif UseFIR && FIR == 4
                    R = rem(length(h),64);
                    if R ~= 0
                        newlength = ceil(length(h)/64)*64;
                        addlength = (newlength-length(h))/2;
                        h = [zeros(1,addlength),h,zeros(1,addlength)];
                    end
                end
            end
            
            % There will always be 128 taps output
            numTaps = length(h);
            firTapsPreScale = zeros(1,128);
            firTapsPreScale(1:numTaps) = h;
            firTapsPreScale = obj.determineBestFractionLength(firTapsPreScale,128);
            aTFIR = 1 + ceil(log2(max(firTapsPreScale)));
            bTFIR = 16 - aTFIR;
            firtaps = int16(firTapsPreScale.*(2^bTFIR));
            % Collect stats
            stats = struct;
            stats.Apass = Apass_actual;
            stats.Astop = Astop_actual;

        end
        
        
        function delay = createObjective(obj,type)
            
            if strcmpi(type,'RX')
                w1 = -obj.RxFpass;
                phEQ = obj.RxPhaseEQ;
                Fpass = obj.RxFpass;
                Fstop = obj.RxFstop;
                FIR = obj.RxFIR;
                Apass = obj.RxApass;
                Astop = obj.RxAstop;
                FIRdBmin = obj.RxFIRdBmin;
            else
                w1 = -obj.TxFpass;
                phEQ = obj.TxPhaseEQ;
                Fpass = obj.TxFpass;
                Fstop = obj.TxFstop;
                FIR = obj.TxFIR;
                Apass = obj.TxApass;
                Astop = obj.TxAstop;
                FIRdBmin = obj.TxFIRdBmin;
            end
            
            %% Determine response delay
            % Build design grids
            Nw = 2048;
            w = w1 - 2*w1.*(2:Nw)/Nw;
            w = [w1 w];
            % Design
            response = generateCombinedResponse(obj, w, type);
            % Measure delay
            invariance = real(response).^2+imag(response).^2;
            phi = unwrap(angle(response));
            sigma = sum(invariance);
            sigmax = sum(w.*invariance);
            sigmay = sum(phi.*invariance);
            sigmaxx = sum(w.*w.*invariance);
            sigmaxy = sum(w.*phi.*invariance);
            delta = sigma*sigmaxx-sigmax^2;
            b = (sigma*sigmaxy-sigmax*sigmay)/delta;
            if phEQ == 0 || phEQ == -1
                delay = -b/(2*pi);
            else
                delay = phEQ*(1e-9);
            end
            
            
            %% Build requirements for FIR
            clkFIR = obj.DataRate*FIR;
            Gpass = floor(obj.G*Fpass/clkFIR);
            Gstop=ceil(obj.G*Fstop/clkFIR);
            Gpass = min(Gpass,Gstop-1);
            fg = (0:Gpass)/obj.G;
            omega = fg*clkFIR;
            
            %% Passband weights
            % Generate response for passband
            response = generateCombinedResponse(obj, omega, type);
            rg2 = exp(-1i*2*pi*omega*delay);
            rg = rg2./response;
            w = abs(response)/(obj.dBinv(Apass/2)-1);
            
            %% Stopband weights
            g = Gpass+1;
            % Expand memory correctly
            fg2 = zeros(1,length(Gstop:(obj.G/2))+length(fg));
            fg2(1:length(fg)) = fg;
            omega2 = zeros(1,length(Gstop:(obj.G/2))+length(omega));
            omega2(1:length(omega)) = omega;
            rgN = complex(zeros(1,length(Gstop:(obj.G/2))+length(rg)));
            rgN(1:length(rg)) = rg;
            % stop band
            for m = Gstop:(obj.G/2)
                g = g+1;
                fg2(g) = m/obj.G;
                omega2(g) = fg2(g)*clkFIR;
                rgN(g) = 0;
            end
            % Generate response for stopband
            wg1 = abs(generateCombinedResponse(obj, omega2(Gpass+2:end), type));
            if strcmpi(type, 'RX')
                wg2 = (wg1)/(obj.dBinv(-Astop));
            else
                wg2 = (sqrt(FIR)*wg1)/(obj.dBinv(-Astop));
            end
            wg3 = obj.dBinv(FIRdBmin);
            wg = max(wg2,wg3);

            %% Combine passband and stopband weights
            weight = [w wg];
            weight = weight/max(weight);
            
            %% Set up design for FIR filter
            grid = fg2;
            if phEQ == -1
                resp = abs(rgN);
            else
                resp = rgN;
            end
            cr = real(resp); %#ok<*PROPLC>
            
            %% Save design criteria
            obj.F1 = grid(1:Gpass+1)*2;
            obj.F2 = grid(Gpass+2:end)*2;
            obj.A1 = cr(1:Gpass+1);
            obj.A2 = cr(Gpass+2:end);
            obj.W1 = weight(1:Gpass+1);
            obj.W2 = weight(Gpass+2:end);
            
            obj.Gpass = Gpass;
            obj.Gstop = Gstop;
            
            obj.grid = grid;
            obj.weight = weight;
            obj.resp = resp;
            obj.omega2 = omega2;
        end
        
        
        function response = generateCombinedResponse(obj, weights, type, addFIR, FIRTaps)
            
            highestSampleRate = obj.ADCRate;
            N = length(weights);
            
            if nargin==3
                addFIR = false;
            end
            
            %% Analog sections
            %%FIXME Do we use TxAnalogCuttoffFrequencyActual or TxAnalogCuttoffFrequency
            %%FIXME Do we use RxAnalogCuttoffFrequencyActual or RxAnalogCuttoffFrequency
            
            if strcmpi(type,'TX')
                % Define analog filter models
                wreal = obj.TxAnalogCuttoffFrequency*(5.0/1.6);
                
                % Define the analog filters (for design purpose)
                [b1,a1] = butter_cg(3,2*pi*obj.TxAnalogCuttoffFrequency,'s');% 3rd order
                [b2,a2] = butter_cg(1,2*pi*wreal,'s');  % 1st order
                TXAnalog = ...
                    sinc(weights/obj.DACRate).*...
                    freqs_cg(b1,a1,2*pi*weights).*...
                    freqs_cg(b2,a2,2*pi*weights);
            else
                wTIA = obj.RxAnalogCuttoffFrequency*(2.5/1.4);
                % Define the analog filters (for design purpose)
                [b1,a1] = butter_cg(1,2*pi*wTIA,'s');  % 1st order
                [b2,a2] = butter_cg(3,2*pi*obj.RxAnalogCuttoffFrequency,'s');% 3rd order
                RXAnalog = ...
                    (sinc(weights/obj.ADCRate).^3).* ...
                    freqs_cg(b1,a1,2*pi*weights).*...
                    freqs_cg(b2,a2,2*pi*weights);
            end
            
            %% Digital sections
            if strcmpi(type,'RX')
                %% RX
                upsampleFactor = 1;
                responseRX = ones(1,N);
                
                if obj.RxHB3 == 2
                    responseRX = freqz_cg(obj.us(obj.RxHB3Coeff,upsampleFactor), 1, weights, highestSampleRate);
                    upsampleFactor = upsampleFactor * 2;
                elseif obj.RxHB3 == 3
                    responseRX = freqz_cg(obj.us(obj.RxHB3CoeffDec3,upsampleFactor), 1, weights, highestSampleRate);
                    upsampleFactor = upsampleFactor * 3;
                end
                if obj.RxHB2 == 2
                    responseRX = responseRX .* freqz_cg(obj.us(obj.RxHB2Coeff,upsampleFactor), 1, weights, highestSampleRate);
                    upsampleFactor = upsampleFactor * 2;
                end
                if obj.RxHB1 == 2
                    responseRX = responseRX .* freqz_cg(obj.us(obj.RxHB1Coeff,upsampleFactor), 1, weights, highestSampleRate);
                    upsampleFactor = upsampleFactor * 2;
                end
                
                response = responseRX .* RXAnalog;
                                
            else
                %% TX
                upsampleFactor = 1;
                responseTX = ones(1,N);
                highestSampleRate = highestSampleRate/obj.DACDivider;
                
                if obj.TxHB3 == 2
                    responseTX = freqz_cg(obj.us(obj.TxHB3Coeff,upsampleFactor), 1, weights, highestSampleRate);
                    upsampleFactor = upsampleFactor * 2;
                elseif obj.TxHB3 == 3
                    responseTX = freqz_cg(obj.us(obj.TxHB3CoeffInt3,upsampleFactor), 1, weights, highestSampleRate);
                    upsampleFactor = upsampleFactor * 3;
                end
                if obj.TxHB2 == 2
                    responseTX = responseTX .* freqz_cg(obj.us(obj.TxHB2Coeff,upsampleFactor), 1, weights, highestSampleRate);
                    upsampleFactor = upsampleFactor * 2;
                end
                if obj.TxHB1 == 2
                    responseTX = responseTX .* freqz_cg(obj.us(obj.TxHB1Coeff,upsampleFactor), 1, weights, highestSampleRate);
                    upsampleFactor = upsampleFactor * 2;
                end
                
                response = responseTX .* TXAnalog;
            
            end
            
            if addFIR
                response = response .* freqz_cg(obj.us(FIRTaps,upsampleFactor),1,weights, highestSampleRate);
            end

            
        end
        
        
    end
    
    methods (Access = protected)
    end
    
    methods (Static)
        function u = us(o,n)
            u = upsample(o,n);
            u = u(1:end-n+1);
        end
        function dBoutput = dBinv(dBinput)
            dBmin = -150;
            if dBinput>dBmin
                dBoutput = 10^(dBinput/20);
            else
                dBoutput = 0;
            end
        end
        function tapsR = determineBestFractionLength(taps,M)
            % Codegen workaround for fixed fi call requirements
            signed = true; wordlength = 16;
            org = taps;
            e = zeros(16,1);
            k=1;
            r  = zeros(16,M);
            r(k,1:M) = double(fi(taps,signed,wordlength,1));
            e(1) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,2));
            e(2) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,3));
            e(3) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,4));
            e(4) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,5));
            e(5) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,6));
            e(6) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,7));
            e(7) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,8));
            e(8) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,9));
            e(9) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,10));
            e(10) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,11));
            e(11) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,12));
            e(12) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,13));
            e(13) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,14));
            e(14) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,15));
            e(15) = sum(abs(r(k,1:M)-org));k = k+1;
            r(k,1:M) = double(fi(taps,signed,wordlength,16));
            e(16) = sum(abs(r(k,1:M)-org));k = k+1;
            [~,fractionLength] = min(e);
            tapsR = r(fractionLength,1:M);
        end
    end
    
end


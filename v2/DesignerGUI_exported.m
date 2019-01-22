classdef DesignerGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ad936xFilterWizardUIFigure     matlab.ui.Figure
        TXPathPanel                    matlab.ui.container.Panel
        HB1DropDownLabel               matlab.ui.control.Label
        TxHB1                          matlab.ui.control.DropDown
        HB2DropDownLabel               matlab.ui.control.Label
        TxHB2                          matlab.ui.control.DropDown
        HB3DropDownLabel               matlab.ui.control.Label
        TxHB3                          matlab.ui.control.DropDown
        FIRDropDownLabel               matlab.ui.control.Label
        TxFIR                          matlab.ui.control.DropDown
        TxRateFIR                      matlab.ui.control.NumericEditField
        TxRateHB1                      matlab.ui.control.NumericEditField
        TxRateHB2                      matlab.ui.control.NumericEditField
        TxRateHB3                      matlab.ui.control.NumericEditField
        RXPathPanel                    matlab.ui.container.Panel
        HB1DropDown_2Label             matlab.ui.control.Label
        RxHB1                          matlab.ui.control.DropDown
        HB2DropDown_2Label             matlab.ui.control.Label
        RxHB2                          matlab.ui.control.DropDown
        HB3DropDown_2Label             matlab.ui.control.Label
        RxHB3                          matlab.ui.control.DropDown
        FIRDropDown_2Label             matlab.ui.control.Label
        RxFIR                          matlab.ui.control.DropDown
        RxRateFIR                      matlab.ui.control.NumericEditField
        RxRateHB1                      matlab.ui.control.NumericEditField
        RxRateHB2                      matlab.ui.control.NumericEditField
        RxRateHB3                      matlab.ui.control.NumericEditField
        ClocksPanel                    matlab.ui.container.Panel
        DACDividerDropDownLabel        matlab.ui.control.Label
        DACDivider                     matlab.ui.control.DropDown
        PLLDividerDropDownLabel        matlab.ui.control.Label
        PLLDivider                     matlab.ui.control.DropDown
        StatusLight                    matlab.ui.control.Lamp
        ErrorStatus                    matlab.ui.control.Label
        MHzLabel_13                    matlab.ui.control.Label
        DACRateEditField_2Label        matlab.ui.control.Label
        DACRate                        matlab.ui.control.NumericEditField
        MHzLabel_14                    matlab.ui.control.Label
        ADCRateEditFieldLabel          matlab.ui.control.Label
        ADCRate                        matlab.ui.control.NumericEditField
        MHzLabel_15                    matlab.ui.control.Label
        PLLRateEditFieldLabel          matlab.ui.control.Label
        PLLRate                        matlab.ui.control.NumericEditField
        AvailableFIRTapsLabel          matlab.ui.control.Label
        AvailableFIRTaps               matlab.ui.control.NumericEditField
        DeviceSettingsPanel            matlab.ui.container.Panel
        DataRateEditFieldLabel         matlab.ui.control.Label
        DataRate                       matlab.ui.control.NumericEditField
        DeviceDropDownLabel            matlab.ui.control.Label
        DeviceDropDown                 matlab.ui.control.DropDown
        ProfileDropDownLabel           matlab.ui.control.Label
        ProfileDropDown                matlab.ui.control.DropDown
        AutoConfigureHBFIRButton       matlab.ui.control.Button
        MHzLabel_21                    matlab.ui.control.Label
        TabGroup                       matlab.ui.container.TabGroup
        TXFilterTab                    matlab.ui.container.Tab
        UIAxesTX                       matlab.ui.control.UIAxes
        FilterSpecificationsPanel      matlab.ui.container.Panel
        FrequencyPanel                 matlab.ui.container.Panel
        FpassEditFieldLabel            matlab.ui.control.Label
        TxFpass                        matlab.ui.control.NumericEditField
        FstopEditFieldLabel            matlab.ui.control.Label
        TxFstop                        matlab.ui.control.NumericEditField
        AnalogFcutoffEditFieldLabel    matlab.ui.control.Label
        TxAnalogCuttoffFrequency       matlab.ui.control.NumericEditField
        MHzLabel                       matlab.ui.control.Label
        MHzLabel_3                     matlab.ui.control.Label
        MHzLabel_2                     matlab.ui.control.Label
        MagnitudePanel                 matlab.ui.container.Panel
        ApassEditFieldLabel            matlab.ui.control.Label
        TxApass                        matlab.ui.control.NumericEditField
        dBLabel                        matlab.ui.control.Label
        dBLabel_2                      matlab.ui.control.Label
        AstopEditFieldLabel            matlab.ui.control.Label
        TxAstop                        matlab.ui.control.NumericEditField
        AdvancedPanel                  matlab.ui.container.Panel
        TargetDelayEditFieldLabel      matlab.ui.control.Label
        TxphEQ                         matlab.ui.control.NumericEditField
        EnTxphEQ                       matlab.ui.control.CheckBox
        nsLabel                        matlab.ui.control.Label
        DesignResultsTX                matlab.ui.control.Table
        DesignTXFilterButton           matlab.ui.control.Button
        ExportDesignPanelTX            matlab.ui.container.Panel
        ToWorkspaceButtonTX            matlab.ui.control.Button
        ToftrFileButtonTX              matlab.ui.control.Button
        ToNoOSFileButtonTX             matlab.ui.control.Button
        ProgramBoardButtonTX           matlab.ui.control.Button
        IPAddressEditField_2Label      matlab.ui.control.Label
        IPAddressEditField_2           matlab.ui.control.EditField
        RXFilterTab                    matlab.ui.container.Tab
        UIAxesRX                       matlab.ui.control.UIAxes
        FilterSpecificationsPanel_2    matlab.ui.container.Panel
        FrequencyPanel_2               matlab.ui.container.Panel
        FpassEditField_2Label          matlab.ui.control.Label
        RxFpass                        matlab.ui.control.NumericEditField
        FstopEditField_2Label          matlab.ui.control.Label
        RxFstop                        matlab.ui.control.NumericEditField
        AnalogFcutoffEditField_2Label  matlab.ui.control.Label
        RxAnalogCuttoffFrequency       matlab.ui.control.NumericEditField
        MHzLabel_18                    matlab.ui.control.Label
        MHzLabel_19                    matlab.ui.control.Label
        MHzLabel_20                    matlab.ui.control.Label
        MagnitudePanel_2               matlab.ui.container.Panel
        ApassEditField_2Label          matlab.ui.control.Label
        RxApass                        matlab.ui.control.NumericEditField
        dBLabel_3                      matlab.ui.control.Label
        dBLabel_4                      matlab.ui.control.Label
        AstopEditField_2Label          matlab.ui.control.Label
        RxAstop                        matlab.ui.control.NumericEditField
        AdvancedPanel_2                matlab.ui.container.Panel
        TargetDelayEditField_2Label    matlab.ui.control.Label
        TargetDelayEditField_2         matlab.ui.control.NumericEditField
        PhaseEqualizeCheckBox_2        matlab.ui.control.CheckBox
        nsLabel_2                      matlab.ui.control.Label
        DesignResultsRX                matlab.ui.control.Table
        DesignRXFilterButton           matlab.ui.control.Button
        ExportDesignPanelRX            matlab.ui.container.Panel
        ToWorkspaceButtonRX            matlab.ui.control.Button
        ToftrFileButtonRX              matlab.ui.control.Button
        ToNoOSFileButtonRX             matlab.ui.control.Button
        ProgramBoardButtonRX           matlab.ui.control.Button
        IPAddressEditFieldLabel        matlab.ui.control.Label
        IPAddressEditField             matlab.ui.control.EditField
        MHzLabel_16                    matlab.ui.control.Label
        MHzLabel_17                    matlab.ui.control.Label
        Button                         matlab.ui.control.Button
    end


    properties (Access = private)
        FilterDesigner % Description
        Black = [0,0,0];
        Red = [1,0,0];
        Green = [0,1,0];
        Grey = [0.65,0.65,0.65];
        RxFIRTaps % Description
        TxFIRTaps % Description
        RxDesignComplete = false;
        TxDesignComplete = false;
    end

    methods (Access = private)
    
        function [pass,E] = updateFilterDesigner(app,parameterName)
            E = [];
            try
                needConv = {'DataRate','TxFpass','TxFstop','RxFpass','RxFstop','TxAnalogCuttoffFrequency'};
                if isa(app.(parameterName).Value,'char')
                    v = str2double(app.(parameterName).Value);
                else
                    v = app.(parameterName).Value;
                end
                if ismember(parameterName,needConv)
                    app.FilterDesigner.(parameterName) = v*1e6;
                else
                    app.FilterDesigner.(parameterName) = v;
                end
                pass = false;
            catch E
                pass = true;
            end
        end
        
        function updateGUIRates(app)
            attributesToUpdate = {...
                'TxRateHB1','TxRateHB2','TxRateHB3','TxRateFIR',...
                'RxRateHB1','RxRateHB2','RxRateHB3','RxRateFIR',...
                'PLLRate','ADCRate','DACRate','PLLRate'};
            
            checks = app.FilterDesigner.validatePathRates();
                
            for attr = attributesToUpdate
                app.([attr{:}]).Value = app.FilterDesigner.(attr{:})/1e6;
                % Color rate if invalid
                if getfield(getfield(checks,attr{:}),'pass')
                    app.([attr{:}]).FontColor = app.Black;
                else
                    app.([attr{:}]).FontColor = app.Red;
                end
            end
            app.updateConfigStatus();
        end
    
        function updateGUIHBs(app)
            attributesToUpdate = {...
                'TxHB1','TxHB2','TxHB3','TxFIR',...
                'RxHB1','RxHB2','RxHB3','RxFIR',...
                'PLLDivider','DACDivider'};
            
            for attr = attributesToUpdate
                app.([attr{:}]).Value = num2str(app.FilterDesigner.(attr{:}));
            end
            app.updateConfigStatus();
        end
    
        function updateGUIFilterConfigs(app)
            attributesToUpdate = {...
                'TxApass','TxAstop',...
                'RxApass','RxAstop'};
            
            for attr = attributesToUpdate
                app.([attr{:}]).Value = (app.FilterDesigner.(attr{:}));
            end
            
            attributesToUpdate = {...
                'TxFpass','TxFstop','TxAnalogCuttoffFrequency',...
                'RxFpass','RxFstop','RxAnalogCuttoffFrequency'};
            
            for attr = attributesToUpdate
                app.([attr{:}]).Value = (app.FilterDesigner.(attr{:}))/1e6;
            end
        end
    
        
        function updateConfigStatus(app)
            % Check if configuration is valid
            validFlag = app.FilterDesigner.ValidConfiguration;
            if validFlag
                app.StatusLight.Color = app.Green;
                app.ErrorStatus.Text = '';
                app.DesignRXFilterButton.Enable = true;
                app.DesignTXFilterButton.Enable = true;
            else
                app.StatusLight.Color = app.Red;
                e = app.FilterDesigner.getConfigurationError();
                app.ErrorStatus.Text = e.msg;
                app.DesignRXFilterButton.Enable = false;
                app.DesignTXFilterButton.Enable = false;
            end
            app.updateAvailableTaps();
        end
    
        
        function updateAvailableTaps(app)
            app.AvailableFIRTaps.Value = app.FilterDesigner.AvailableFIRTaps;
        end
        
        function UpdatePlot(app,FilterStages,Axes,ConverterRate,Apass,Astop,Fpass,Fstop)
            G = 8192;
            f = linspace(0,app.FilterDesigner.DataRate/2,G);
            r = mag2db(abs(freqz(FilterStages,f,ConverterRate)));
            f = f./1e6;
            plot(Axes, f(1:end-1), r(1:end-1) );           
            hold(Axes,'on');
            line(Axes, [Fpass Fpass]./1e6, [-(Apass/2) -100], 'Color', 'Red');
            line(Axes, [0 Fpass/1e6], [-(Apass/2) -(Apass/2)], 'Color', 'Red');
            line(Axes, [0 Fstop/1e6], [Apass/2 Apass/2], 'Color', 'Red');
            line(Axes, [Fstop Fstop]./1e6, [Apass/2 -Astop], 'Color', 'Red');
            line(Axes, [Fstop app.FilterDesigner.DataRate]./1e6, [-Astop -Astop], 'Color', 'Red');
            xlim(Axes, [0 app.FilterDesigner.DataRate/2]./1e6);
            ylim(Axes, [-inf, inf]);
            hold(Axes,'off');
            
            UpdateExportButtons(app);
        end
        
        function SaveToFTR(app)
            [filename,path] = uiputfile('*.ftr', 'Save coefficients as');
            if filename == 0
                return;
            else
                newpath = strcat(path,filename);
            end
            
            fid = fopen(newpath,'w');
            
            fprintf(fid, '# Generated with AD9361 Filter Design Wizard %s\r\n', get_version);
            fprintf(fid, '# MATLAB %s, %s\r\n', version(), datestr(now()));
            fprintf(fid, '# Inputs:\r\n');
            
            %FIXME
            %converter_rate = get_converter_clk(handles);
            %converter_rate = get_ADC_clk(handles);
            
%             PLL_rate = str2double(get(handles.Pll_rate, 'String')) * 1e6;
%             [rx_FIR_rate, rx_HB1_rate, rx_HB2_rate, rx_HB3_rate, rx_RFbw_hw, ...
%                 tx_FIR_rate, tx_HB1_rate, tx_HB2_rate, tx_HB3_rate, tx_RFbw_hw] = get_path_rates(handles);
            
            %fprintf(fid, '# PLL CLK Frequency = %f Hz\r\n', pll_rate);
            %fprintf(fid, '# Converter Sample Frequency = %f Hz\r\n', converter_rate);
            fprintf(fid, '# Data Sample Frequency = %.0f Hz\r\n', app.FilterDesigner.DataRate);
            if app.FilterDesigner.RxPhaseEQ
                fprintf(fid, '# RX Phase equalization = %f ns\r\n', app.FilterDesigner.TxTargetDelay);
            end
            if app.FilterDesigner.RxPhaseEQ
                fprintf(fid, '# TX Phase equalization = %f ns\r\n', app.FilterDesigner.RxTargetDelay);
            end
            fprintf(fid, 'TX 3 GAIN %d INT %d\r\n', handles.tx.gain, app.FilterDesigner.TxFIR);
            fprintf(fid, 'RX 3 GAIN %d DEC %d\r\n', handles.rx.gain, app.FilterDesigner.RxFIR);
            fprintf(fid, 'RTX %.0f %.0f %.0f %.0f %.0f %.0f\r\n', app.FilterDesigner.PLLRate, app.FilterDesigner.TxRateHB3,...
                app.FilterDesigner.TxRateHB2, app.FilterDesigner.TxRateHB1, app.FilterDesigner.TxRateFIR, app.FilterDesigner.DataRate);
            fprintf(fid, 'RRX %.0f %.0f %.0f %.0f %.0f %.0f\r\n', app.FilterDesigner.PLLRate, app.FilterDesigner.RxRateHB3,...
                app.FilterDesigner.RxRateHB2, app.FilterDesigner.RxRateHB1, app.FilterDesigner.RxRateFIR, app.FilterDesigner.DataRate);
            fprintf(fid, 'BWTX %.0f\r\n', app.FilterDesigner.AnalogCuttoffFrequencyTxHW);
            fprintf(fid, 'BWRX %.0f\r\n', app.FilterDesigner.AnalogCuttoffFrequencyRxHW);
            
            % concat and transform Rx and Tx coefficient matrices for output
            coefficients = flip(rot90(vertcat(handles.tfirtaps, handles.rfirtaps)));
            
            % output all non-zero coefficients since they're padded to 128 with zeros
            for i = 1:handles.nfirtaps
                fprintf(fid, '%d,%d\r\n', coefficients(i,:));
            end
            
            fclose(fid);
        end
        
        function t = createTable(~,data)
            fn = fieldnames(data);
            fd = [];
            for f = 1:length(fn)
                fd = [fd;data.(fn{f})];
            end
            t = table(fn,fd);
        end
        
        function SaveToNoOS(app)
            
            [filename,path] = uiputfile('*.c', 'Save coefficients as');
            if filename == 0
                error('No filename provided');
            else
                newpath = strcat(path, filename);
            end
            fid = fopen(newpath,'w');
            
            fprintf(fid, '// Generated with AD9361 Filter Design Wizard %s\n', get_version);
            fprintf(fid, '// MATLAB %s, %s\n', version(), datestr(now()));
            fprintf(fid, '// Inputs:\n');
            
            fprintf(fid, '// Data Sample Frequency = %.0f Hz\n', app.FilterDesigner.DataRate);
            if app.FilterDesigner.RxPhaseEQ
                fprintf(fid, '// RX Phase equalization = %f ns\n', app.FilterDesigner.RxPhaseEQ);
            end
            if app.FilterDesigner.TxPhaseEQ
                fprintf(fid, '// TX Phase equalization = %f ns\n', app.FilterDesigner.TxPhaseEQ);
            end
            
            % Rx
            fprintf(fid, '\nAD9361_RXFIRConfig rx_fir_config = {\n');
            fprintf(fid, '\t3, // rx\n');
            fprintf(fid, '\t%d, // rx_gain\n', handles.rx.gain);
            fprintf(fid, '\t%d, // rx_dec\n', app.FilterDesigner.RxFIR);
            
            coefficients = sprintf('%.0f,', flip(rot90(handles.rfirtaps)));
            coefficients = coefficients(1:end-1); % strip final comma
            fprintf(fid, '\t{%s}, // rx_coef[128]\n', coefficients);
            
            fprintf(fid, '\t%d, // rx_coef_size\n', handles.nfirtaps);
            fprintf(fid, '\t{%.0f,%.0f,%.0f,%.0f,%.0f,%.0f}, // rx_path_clks[6]\n', ...
                app.FilterDesigner.PLLRate, app.FilterDesigner.RxRateHB3, app.FilterDesigner.RxRateHB2,...
                app.FilterDesigner.RxRateHB1, app.FilterDesigner.RxRateFIR, app.FilterDesigner.DataRate);
            fprintf(fid, '\t%.0f // rx_bandwidth\n', app.FilterDesigner.RxAnalogCuttoffFrequencyActual);
            fprintf(fid, '};\n\n');
            
            % Tx
            fprintf(fid, 'AD9361_TXFIRConfig tx_fir_config = {\n');
            fprintf(fid, '\t3, // tx\n');
            fprintf(fid, '\t%d, // tx_gain\n', handles.tx.gain);
            fprintf(fid, '\t%d, // tx_int\n', app.FilterDesigner.TxFIR);
            coefficients = sprintf('%.0f,', flip(rot90(handles.tfirtaps)));
            coefficients = coefficients(1:end-1); % strip final comma
            fprintf(fid, '\t{%s}, // tx_coef[128]\n', coefficients);
            fprintf(fid, '\t%d, // tx_coef_size\n', handles.nfirtaps);
            fprintf(fid, '\t{%.0f,%.0f,%.0f,%.0f,%.0f,%.0f}, // tx_path_clks[6]\n', ...
                app.FilterDesigner.PLLRate, app.FilterDesigner.TxRateHB3, app.FilterDesigner.TxRateHB2,...
                app.FilterDesigner.TxRateHB1, app.FilterDesigner.TxRateFIR, app.FilterDesigner.DataRate);
            fprintf(fid, '\t%.0f // tx_bandwidth\n', app.FilterDesigner.TxAnalogCuttoffFrequencyActual);
            fprintf(fid, '};\n');
            
            fclose(fid);
            
        end
        
        function UpdateExportButtons(app)
            if app.FIRDesignComplete()
                app.ToWorkspaceButtonRX.Enable = true;
                app.ToWorkspaceButtonTX.Enable = true;
                app.ToftrFileButtonRX.Enable = true;
                app.ToftrFileButtonTX.Enable = true;
                app.ToNoOSFileButtonRX.Enable = true;
                app.ToNoOSFileButtonTX.Enable = true;
                app.ProgramBoardButtonRX.Enable = true;
                app.ProgramBoardButtonTX.Enable = true;
            else
                app.ToWorkspaceButtonRX.Enable = false;
                app.ToWorkspaceButtonTX.Enable = false;
                app.ToftrFileButtonRX.Enable = false;
                app.ToftrFileButtonTX.Enable = false;
                app.ToNoOSFileButtonRX.Enable = false;
                app.ToNoOSFileButtonTX.Enable = false;
                app.ProgramBoardButtonRX.Enable = false;
                app.ProgramBoardButtonTX.Enable = false;
            end
            
        end
        
        function FixDesignerInputs(app)
            
            if app.TxFstop.Value > app.DataRate.Value/2
                app.TxFstop.Value = app.DataRate.Value/2;
            end
            if app.TxFpass.Value > app.TxFstop.Value
                app.TxFpass.Value = app.TxFstop.Value;
            end
            
            if app.RxFstop.Value > app.DataRate.Value/2
                app.RxFstop.Value = app.DataRate.Value/2;
            end
            if app.RxFpass.Value > app.RxFstop.Value
                app.RxFpass.Value = app.RxFstop.Value;
            end
            
        end
    end
    
    methods (Access = public)
        
        function value = FIRDesignComplete(app)
            value = app.RxDesignComplete && app.TxDesignComplete;
        end
    end


    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            % Setup HB configuration            
            app.FilterDesigner = ad936xFilterDesigner();
            % Update rates on GUI based on filter designer
            app.updateGUIRates();
            % Update HBs on GUI based on filter designer
            app.updateGUIHBs();
            % Update data rate
            app.DataRate.Value = app.FilterDesigner.DataRate/1e6;
            % Update filter configs
            app.updateGUIFilterConfigs();
            % Update export buttons
            app.UpdateExportButtons();
        end

        % Value changed function: DACDivider, DataRate, PLLDivider, 
        % RxFIR, RxHB1, RxHB2, RxHB3, TxFIR, TxHB1, TxHB2, TxHB3
        function HBValueChanged(app, event)
            attributesToUpdate = {...
                'TxHB1','TxHB2','TxHB3','TxFIR',...
                'RxHB1','RxHB2','RxHB3','RxFIR',...
                'DACDivider','PLLDivider','DataRate'};
            app.RxDesignComplete = false;
            app.TxDesignComplete = false;
            app.UpdateExportButtons();
            % Pass new configuration into Filter Designer for parameter validation
            for attr = attributesToUpdate
                app.updateFilterDesigner(attr{:});
            end

            % Check if configuration is valid
            app.updateConfigStatus();
            app.updateGUIRates();   

        end

        % Value changed function: RxAnalogCuttoffFrequency, RxApass, 
        % RxAstop, RxFpass, RxFstop, TxAnalogCuttoffFrequency, 
        % TxApass, TxAstop, TxFpass, TxFstop
        function DesignConfigValueChanged(app, event)
            
            FixDesignerInputs(app);
            
            attributesToUpdate = {...
                'TxApass','TxAstop','TxFpass','TxFstop','TxAnalogCuttoffFrequency'...
                'RxApass','RxAstop','RxFpass','RxFstop','RxAnalogCuttoffFrequency'};
                
            % Pass new configuration into Filter Designer for parameter validation
            for attr = attributesToUpdate
                app.updateFilterDesigner(attr{:});
            end
            
            % Check if configuration is valid
            app.updateConfigStatus();
            
        end

        % Value changed function: EnTxphEQ
        function EnTxphEQValueChanged(app, event)
            value = app.EnTxphEQ.Value;
            if value
                app.TxphEQ.Editable = true;
                app.FilterDesigner.TxPhaseEQ = 1;
            else
                app.TxphEQ.Editable = false;
                app.FilterDesigner.TxPhaseEQ = 0;
            end                
        end

        % Value changed function: TxphEQ
        function TxphEQValueChanged(app, event)
            app.FilterDesigner.TxTargetDelay = app.TxphEQ.Value;
        end

        % Button pushed function: DesignTXFilterButton
        function DesignTXFilterButtonPushed(app, event)
            app.DesignTXFilterButton.BackgroundColor = app.Grey;
            pause(0.1);
            [pass,~,~,fpTaps,stats] = app.FilterDesigner.designFilter('Tx');
            if pass
                app.DesignTXFilterButton.BackgroundColor = app.Green;
            else
                app.DesignTXFilterButton.BackgroundColor = app.Red;
            end
            app.TxDesignComplete = true;
            % Update plot
            FIR = dsp.FIRInterpolator('Numerator',fpTaps,'InterpolationFactor',app.FilterDesigner.TxFIR);
            FilterStages = createObjects(app.FilterDesigner);
            addStage(FilterStages, FIR, 1);
            UpdatePlot(app,FilterStages,app.UIAxesTX,app.FilterDesigner.DACRate,app.FilterDesigner.TxApass,...
                app.FilterDesigner.TxAstop,app.FilterDesigner.TxFpass,app.FilterDesigner.TxFstop);
            % Update design state table
            app.DesignResultsTX.Data = createTable(app,stats);
        end

        % Button pushed function: AutoConfigureHBFIRButton
        function AutoConfigureHBFIRButtonPushed(app, event)
            app.FilterDesigner.AutoSetRates();
            % Set rates
            app.updateGUIRates;
            app.updateGUIHBs;
        end

        % Button pushed function: DesignRXFilterButton
        function DesignRXFilterButtonPushed(app, event)
            app.DesignRXFilterButton.BackgroundColor = app.Grey;
            pause(0.1);
            [pass,app.RxFIRTaps,~,fpTaps,stats] = app.FilterDesigner.designFilter('Rx');
            if pass
                app.DesignRXFilterButton.BackgroundColor = app.Green;
            else
                app.DesignRXFilterButton.BackgroundColor = app.Red;
            end
            app.RxDesignComplete = true;
            % Update plot
            FIR = dsp.FIRDecimator('Numerator',fpTaps,'DecimationFactor',app.FilterDesigner.RxFIR);
            [~,FilterStages] = createObjects(app.FilterDesigner);
            addStage(FilterStages, FIR);
            UpdatePlot(app,FilterStages,app.UIAxesRX,app.FilterDesigner.ADCRate,app.FilterDesigner.RxApass,...
                app.FilterDesigner.RxAstop,app.FilterDesigner.RxFpass,app.FilterDesigner.RxFstop);
            % Update design state table
            app.DesignResultsRX.Data = createTable(app,stats);
        end

        % Button pushed function: ToftrFileButtonRX, ToftrFileButtonTX
        function ToftrFileButton_2Pushed(app, event)
            SaveToFTR(app);
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ad936xFilterWizardUIFigure
            app.ad936xFilterWizardUIFigure = uifigure;
            app.ad936xFilterWizardUIFigure.Position = [100 100 1110 522];
            app.ad936xFilterWizardUIFigure.Name = 'ad936x Filter Wizard';

            % Create TXPathPanel
            app.TXPathPanel = uipanel(app.ad936xFilterWizardUIFigure);
            app.TXPathPanel.Title = 'TX Path';
            app.TXPathPanel.Position = [11 273 160 150];

            % Create HB1DropDownLabel
            app.HB1DropDownLabel = uilabel(app.TXPathPanel);
            app.HB1DropDownLabel.HorizontalAlignment = 'right';
            app.HB1DropDownLabel.Position = [11 71 29 22];
            app.HB1DropDownLabel.Text = 'HB1';

            % Create TxHB1
            app.TxHB1 = uidropdown(app.TXPathPanel);
            app.TxHB1.Items = {'1x', '2x'};
            app.TxHB1.ItemsData = {'1', '2'};
            app.TxHB1.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.TxHB1.Position = [55 71 39 22];
            app.TxHB1.Value = '1';

            % Create HB2DropDownLabel
            app.HB2DropDownLabel = uilabel(app.TXPathPanel);
            app.HB2DropDownLabel.HorizontalAlignment = 'right';
            app.HB2DropDownLabel.Position = [11 38 29 22];
            app.HB2DropDownLabel.Text = 'HB2';

            % Create TxHB2
            app.TxHB2 = uidropdown(app.TXPathPanel);
            app.TxHB2.Items = {'1x', '2x'};
            app.TxHB2.ItemsData = {'1', '2'};
            app.TxHB2.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.TxHB2.Position = [55 38 39 22];
            app.TxHB2.Value = '1';

            % Create HB3DropDownLabel
            app.HB3DropDownLabel = uilabel(app.TXPathPanel);
            app.HB3DropDownLabel.HorizontalAlignment = 'right';
            app.HB3DropDownLabel.Position = [11 6 29 22];
            app.HB3DropDownLabel.Text = 'HB3';

            % Create TxHB3
            app.TxHB3 = uidropdown(app.TXPathPanel);
            app.TxHB3.Items = {'1x', '2x', '3x'};
            app.TxHB3.ItemsData = {'1', '2', '3'};
            app.TxHB3.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.TxHB3.Position = [55 6 39 22];
            app.TxHB3.Value = '1';

            % Create FIRDropDownLabel
            app.FIRDropDownLabel = uilabel(app.TXPathPanel);
            app.FIRDropDownLabel.HorizontalAlignment = 'right';
            app.FIRDropDownLabel.Position = [11 104 25 22];
            app.FIRDropDownLabel.Text = 'FIR';

            % Create TxFIR
            app.TxFIR = uidropdown(app.TXPathPanel);
            app.TxFIR.Items = {'1x', '2x', '4x'};
            app.TxFIR.ItemsData = {'1', '2', '4'};
            app.TxFIR.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.TxFIR.Position = [55 104 39 22];
            app.TxFIR.Value = '1';

            % Create TxRateFIR
            app.TxRateFIR = uieditfield(app.TXPathPanel, 'numeric');
            app.TxRateFIR.Editable = 'off';
            app.TxRateFIR.Position = [103 106 50 20];

            % Create TxRateHB1
            app.TxRateHB1 = uieditfield(app.TXPathPanel, 'numeric');
            app.TxRateHB1.Editable = 'off';
            app.TxRateHB1.Position = [103 73 50 20];

            % Create TxRateHB2
            app.TxRateHB2 = uieditfield(app.TXPathPanel, 'numeric');
            app.TxRateHB2.Editable = 'off';
            app.TxRateHB2.Position = [103 40 50 20];

            % Create TxRateHB3
            app.TxRateHB3 = uieditfield(app.TXPathPanel, 'numeric');
            app.TxRateHB3.Editable = 'off';
            app.TxRateHB3.Position = [103 8 50 20];

            % Create RXPathPanel
            app.RXPathPanel = uipanel(app.ad936xFilterWizardUIFigure);
            app.RXPathPanel.Title = 'RX Path';
            app.RXPathPanel.Position = [172 273 160 150];

            % Create HB1DropDown_2Label
            app.HB1DropDown_2Label = uilabel(app.RXPathPanel);
            app.HB1DropDown_2Label.HorizontalAlignment = 'right';
            app.HB1DropDown_2Label.Position = [9 71 29 22];
            app.HB1DropDown_2Label.Text = 'HB1';

            % Create RxHB1
            app.RxHB1 = uidropdown(app.RXPathPanel);
            app.RxHB1.Items = {'1x', '2x'};
            app.RxHB1.ItemsData = {'1', '2'};
            app.RxHB1.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.RxHB1.Position = [53 71 37 22];
            app.RxHB1.Value = '1';

            % Create HB2DropDown_2Label
            app.HB2DropDown_2Label = uilabel(app.RXPathPanel);
            app.HB2DropDown_2Label.HorizontalAlignment = 'right';
            app.HB2DropDown_2Label.Position = [9 38 29 22];
            app.HB2DropDown_2Label.Text = 'HB2';

            % Create RxHB2
            app.RxHB2 = uidropdown(app.RXPathPanel);
            app.RxHB2.Items = {'1x', '2x'};
            app.RxHB2.ItemsData = {'1', '2'};
            app.RxHB2.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.RxHB2.Position = [53 38 37 22];
            app.RxHB2.Value = '1';

            % Create HB3DropDown_2Label
            app.HB3DropDown_2Label = uilabel(app.RXPathPanel);
            app.HB3DropDown_2Label.HorizontalAlignment = 'right';
            app.HB3DropDown_2Label.Position = [9 6 29 22];
            app.HB3DropDown_2Label.Text = 'HB3';

            % Create RxHB3
            app.RxHB3 = uidropdown(app.RXPathPanel);
            app.RxHB3.Items = {'1x', '2x', '3x'};
            app.RxHB3.ItemsData = {'1', '2', '3'};
            app.RxHB3.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.RxHB3.Position = [53 6 37 22];
            app.RxHB3.Value = '1';

            % Create FIRDropDown_2Label
            app.FIRDropDown_2Label = uilabel(app.RXPathPanel);
            app.FIRDropDown_2Label.HorizontalAlignment = 'right';
            app.FIRDropDown_2Label.Position = [9 104 25 22];
            app.FIRDropDown_2Label.Text = 'FIR';

            % Create RxFIR
            app.RxFIR = uidropdown(app.RXPathPanel);
            app.RxFIR.Items = {'1x', '2x', '4x'};
            app.RxFIR.ItemsData = {'1', '2', '4'};
            app.RxFIR.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.RxFIR.Position = [53 104 37 22];
            app.RxFIR.Value = '1';

            % Create RxRateFIR
            app.RxRateFIR = uieditfield(app.RXPathPanel, 'numeric');
            app.RxRateFIR.Editable = 'off';
            app.RxRateFIR.Position = [98 104 55 22];

            % Create RxRateHB1
            app.RxRateHB1 = uieditfield(app.RXPathPanel, 'numeric');
            app.RxRateHB1.Editable = 'off';
            app.RxRateHB1.Position = [98 71 55 22];

            % Create RxRateHB2
            app.RxRateHB2 = uieditfield(app.RXPathPanel, 'numeric');
            app.RxRateHB2.Editable = 'off';
            app.RxRateHB2.Position = [98 38 55 22];

            % Create RxRateHB3
            app.RxRateHB3 = uieditfield(app.RXPathPanel, 'numeric');
            app.RxRateHB3.Editable = 'off';
            app.RxRateHB3.Position = [98 6 55 22];

            % Create ClocksPanel
            app.ClocksPanel = uipanel(app.ad936xFilterWizardUIFigure);
            app.ClocksPanel.Title = 'Clocks';
            app.ClocksPanel.Position = [11 122 320 151];

            % Create DACDividerDropDownLabel
            app.DACDividerDropDownLabel = uilabel(app.ClocksPanel);
            app.DACDividerDropDownLabel.HorizontalAlignment = 'right';
            app.DACDividerDropDownLabel.Position = [8 108 72 22];
            app.DACDividerDropDownLabel.Text = 'DAC Divider';

            % Create DACDivider
            app.DACDivider = uidropdown(app.ClocksPanel);
            app.DACDivider.Items = {'1x', '2x'};
            app.DACDivider.ItemsData = {'1', '2'};
            app.DACDivider.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.DACDivider.Position = [91 108 50 22];
            app.DACDivider.Value = '1';

            % Create PLLDividerDropDownLabel
            app.PLLDividerDropDownLabel = uilabel(app.ClocksPanel);
            app.PLLDividerDropDownLabel.HorizontalAlignment = 'right';
            app.PLLDividerDropDownLabel.Position = [8 83 68 22];
            app.PLLDividerDropDownLabel.Text = 'PLL Divider';

            % Create PLLDivider
            app.PLLDivider = uidropdown(app.ClocksPanel);
            app.PLLDivider.Items = {'1x', '2x', '4x', '8x', '16x', '32x', '64x'};
            app.PLLDivider.ItemsData = {'1', '2', '4', '8', '16', '32', '64'};
            app.PLLDivider.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.PLLDivider.Position = [91 83 50 22];
            app.PLLDivider.Value = '1';

            % Create StatusLight
            app.StatusLight = uilamp(app.ClocksPanel);
            app.StatusLight.Position = [8 30 20 20];
            app.StatusLight.Color = [1 0 0];

            % Create ErrorStatus
            app.ErrorStatus = uilabel(app.ClocksPanel);
            app.ErrorStatus.FontWeight = 'bold';
            app.ErrorStatus.Position = [5 4 311 21];
            app.ErrorStatus.Text = '';

            % Create MHzLabel_13
            app.MHzLabel_13 = uilabel(app.ClocksPanel);
            app.MHzLabel_13.Position = [272 108 44 22];
            app.MHzLabel_13.Text = 'MHz';

            % Create DACRateEditField_2Label
            app.DACRateEditField_2Label = uilabel(app.ClocksPanel);
            app.DACRateEditField_2Label.HorizontalAlignment = 'right';
            app.DACRateEditField_2Label.Position = [157 108 60 22];
            app.DACRateEditField_2Label.Text = 'DAC Rate';

            % Create DACRate
            app.DACRate = uieditfield(app.ClocksPanel, 'numeric');
            app.DACRate.Editable = 'off';
            app.DACRate.Position = [227 108 42 22];

            % Create MHzLabel_14
            app.MHzLabel_14 = uilabel(app.ClocksPanel);
            app.MHzLabel_14.Position = [272 85 44 22];
            app.MHzLabel_14.Text = 'MHz';

            % Create ADCRateEditFieldLabel
            app.ADCRateEditFieldLabel = uilabel(app.ClocksPanel);
            app.ADCRateEditFieldLabel.HorizontalAlignment = 'right';
            app.ADCRateEditFieldLabel.Position = [157 84 60 22];
            app.ADCRateEditFieldLabel.Text = 'ADC Rate';

            % Create ADCRate
            app.ADCRate = uieditfield(app.ClocksPanel, 'numeric');
            app.ADCRate.Editable = 'off';
            app.ADCRate.Position = [227 84 42 22];

            % Create MHzLabel_15
            app.MHzLabel_15 = uilabel(app.ClocksPanel);
            app.MHzLabel_15.Position = [272 61 44 22];
            app.MHzLabel_15.Text = 'MHz';

            % Create PLLRateEditFieldLabel
            app.PLLRateEditFieldLabel = uilabel(app.ClocksPanel);
            app.PLLRateEditFieldLabel.HorizontalAlignment = 'right';
            app.PLLRateEditFieldLabel.Position = [161 59 56 22];
            app.PLLRateEditFieldLabel.Text = 'PLL Rate';

            % Create PLLRate
            app.PLLRate = uieditfield(app.ClocksPanel, 'numeric');
            app.PLLRate.Editable = 'off';
            app.PLLRate.Position = [227 59 42 22];

            % Create AvailableFIRTapsLabel
            app.AvailableFIRTapsLabel = uilabel(app.ClocksPanel);
            app.AvailableFIRTapsLabel.HorizontalAlignment = 'right';
            app.AvailableFIRTapsLabel.Position = [1 59 105 22];
            app.AvailableFIRTapsLabel.Text = 'Available FIR Taps';

            % Create AvailableFIRTaps
            app.AvailableFIRTaps = uieditfield(app.ClocksPanel, 'numeric');
            app.AvailableFIRTaps.Editable = 'off';
            app.AvailableFIRTaps.Position = [116 59 42 22];

            % Create DeviceSettingsPanel
            app.DeviceSettingsPanel = uipanel(app.ad936xFilterWizardUIFigure);
            app.DeviceSettingsPanel.Title = 'Device Settings';
            app.DeviceSettingsPanel.Position = [11 423 320 90];

            % Create DataRateEditFieldLabel
            app.DataRateEditFieldLabel = uilabel(app.DeviceSettingsPanel);
            app.DataRateEditFieldLabel.HorizontalAlignment = 'right';
            app.DataRateEditFieldLabel.Position = [2 8 56 22];
            app.DataRateEditFieldLabel.Text = 'DataRate';

            % Create DataRate
            app.DataRate = uieditfield(app.DeviceSettingsPanel, 'numeric');
            app.DataRate.Limits = [0.520833 61.44];
            app.DataRate.ValueDisplayFormat = '%5.3f';
            app.DataRate.ValueChangedFcn = createCallbackFcn(app, @HBValueChanged, true);
            app.DataRate.Position = [65 8 59 22];
            app.DataRate.Value = 8;

            % Create DeviceDropDownLabel
            app.DeviceDropDownLabel = uilabel(app.DeviceSettingsPanel);
            app.DeviceDropDownLabel.HorizontalAlignment = 'right';
            app.DeviceDropDownLabel.Position = [2 38 42 22];
            app.DeviceDropDownLabel.Text = 'Device';

            % Create DeviceDropDown
            app.DeviceDropDown = uidropdown(app.DeviceSettingsPanel);
            app.DeviceDropDown.Items = {'AD9361', 'AD9364', 'Pluto'};
            app.DeviceDropDown.Position = [59 38 82 22];
            app.DeviceDropDown.Value = 'AD9361';

            % Create ProfileDropDownLabel
            app.ProfileDropDownLabel = uilabel(app.DeviceSettingsPanel);
            app.ProfileDropDownLabel.HorizontalAlignment = 'right';
            app.ProfileDropDownLabel.Position = [154 38 40 22];
            app.ProfileDropDownLabel.Text = 'Profile';

            % Create ProfileDropDown
            app.ProfileDropDown = uidropdown(app.DeviceSettingsPanel);
            app.ProfileDropDown.Items = {'LTE5', 'LTE10', 'LTE15', 'LTE20', 'Custom'};
            app.ProfileDropDown.Position = [209 38 100 22];
            app.ProfileDropDown.Value = 'LTE5';

            % Create AutoConfigureHBFIRButton
            app.AutoConfigureHBFIRButton = uibutton(app.DeviceSettingsPanel, 'push');
            app.AutoConfigureHBFIRButton.ButtonPushedFcn = createCallbackFcn(app, @AutoConfigureHBFIRButtonPushed, true);
            app.AutoConfigureHBFIRButton.Position = [161 8 150 22];
            app.AutoConfigureHBFIRButton.Text = 'Auto Configure HB/FIR';

            % Create MHzLabel_21
            app.MHzLabel_21 = uilabel(app.DeviceSettingsPanel);
            app.MHzLabel_21.Position = [130 8 44 22];
            app.MHzLabel_21.Text = 'MHz';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.ad936xFilterWizardUIFigure);
            app.TabGroup.Position = [331 13 770 500];

            % Create TXFilterTab
            app.TXFilterTab = uitab(app.TabGroup);
            app.TXFilterTab.Title = 'TX Filter';

            % Create UIAxesTX
            app.UIAxesTX = uiaxes(app.TXFilterTab);
            xlabel(app.UIAxesTX, 'Frequency (MHz)')
            ylabel(app.UIAxesTX, 'Magnitude (dB)')
            app.UIAxesTX.Box = 'on';
            app.UIAxesTX.XGrid = 'on';
            app.UIAxesTX.YGrid = 'on';
            app.UIAxesTX.Position = [3 146 598 327];

            % Create FilterSpecificationsPanel
            app.FilterSpecificationsPanel = uipanel(app.TXFilterTab);
            app.FilterSpecificationsPanel.Title = 'Filter Specifications';
            app.FilterSpecificationsPanel.Position = [0 0 601 140];

            % Create FrequencyPanel
            app.FrequencyPanel = uipanel(app.FilterSpecificationsPanel);
            app.FrequencyPanel.Title = 'Frequency';
            app.FrequencyPanel.Position = [210 1 210 120];

            % Create FpassEditFieldLabel
            app.FpassEditFieldLabel = uilabel(app.FrequencyPanel);
            app.FpassEditFieldLabel.HorizontalAlignment = 'right';
            app.FpassEditFieldLabel.Position = [59 73 38 22];
            app.FpassEditFieldLabel.Text = 'Fpass';

            % Create TxFpass
            app.TxFpass = uieditfield(app.FrequencyPanel, 'numeric');
            app.TxFpass.Limits = [0 Inf];
            app.TxFpass.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.TxFpass.Position = [112 73 60 22];

            % Create FstopEditFieldLabel
            app.FstopEditFieldLabel = uilabel(app.FrequencyPanel);
            app.FstopEditFieldLabel.HorizontalAlignment = 'right';
            app.FstopEditFieldLabel.Position = [61 47 36 22];
            app.FstopEditFieldLabel.Text = 'Fstop';

            % Create TxFstop
            app.TxFstop = uieditfield(app.FrequencyPanel, 'numeric');
            app.TxFstop.Limits = [0 Inf];
            app.TxFstop.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.TxFstop.Position = [112 47 60 22];

            % Create AnalogFcutoffEditFieldLabel
            app.AnalogFcutoffEditFieldLabel = uilabel(app.FrequencyPanel);
            app.AnalogFcutoffEditFieldLabel.HorizontalAlignment = 'right';
            app.AnalogFcutoffEditFieldLabel.Position = [14 20 83 22];
            app.AnalogFcutoffEditFieldLabel.Text = 'Analog Fcutoff';

            % Create TxAnalogCuttoffFrequency
            app.TxAnalogCuttoffFrequency = uieditfield(app.FrequencyPanel, 'numeric');
            app.TxAnalogCuttoffFrequency.Limits = [0 Inf];
            app.TxAnalogCuttoffFrequency.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.TxAnalogCuttoffFrequency.Position = [112 20 60 22];

            % Create MHzLabel
            app.MHzLabel = uilabel(app.FrequencyPanel);
            app.MHzLabel.Position = [175 73 30 22];
            app.MHzLabel.Text = 'MHz';

            % Create MHzLabel_3
            app.MHzLabel_3 = uilabel(app.FrequencyPanel);
            app.MHzLabel_3.Position = [175 20 30 22];
            app.MHzLabel_3.Text = 'MHz';

            % Create MHzLabel_2
            app.MHzLabel_2 = uilabel(app.FrequencyPanel);
            app.MHzLabel_2.Position = [175 47 30 22];
            app.MHzLabel_2.Text = 'MHz';

            % Create MagnitudePanel
            app.MagnitudePanel = uipanel(app.FilterSpecificationsPanel);
            app.MagnitudePanel.Title = 'Magnitude';
            app.MagnitudePanel.Position = [0 1 210 120];

            % Create ApassEditFieldLabel
            app.ApassEditFieldLabel = uilabel(app.MagnitudePanel);
            app.ApassEditFieldLabel.HorizontalAlignment = 'right';
            app.ApassEditFieldLabel.Position = [35 73 39 22];
            app.ApassEditFieldLabel.Text = 'Apass';

            % Create TxApass
            app.TxApass = uieditfield(app.MagnitudePanel, 'numeric');
            app.TxApass.Limits = [0 Inf];
            app.TxApass.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.TxApass.Position = [89 73 60 22];

            % Create dBLabel
            app.dBLabel = uilabel(app.MagnitudePanel);
            app.dBLabel.Position = [152 73 25 22];
            app.dBLabel.Text = 'dB';

            % Create dBLabel_2
            app.dBLabel_2 = uilabel(app.MagnitudePanel);
            app.dBLabel_2.Position = [152 48 25 22];
            app.dBLabel_2.Text = 'dB';

            % Create AstopEditFieldLabel
            app.AstopEditFieldLabel = uilabel(app.MagnitudePanel);
            app.AstopEditFieldLabel.HorizontalAlignment = 'right';
            app.AstopEditFieldLabel.Position = [38 48 36 22];
            app.AstopEditFieldLabel.Text = 'Astop';

            % Create TxAstop
            app.TxAstop = uieditfield(app.MagnitudePanel, 'numeric');
            app.TxAstop.Limits = [0 Inf];
            app.TxAstop.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.TxAstop.Position = [89 48 60 22];

            % Create AdvancedPanel
            app.AdvancedPanel = uipanel(app.FilterSpecificationsPanel);
            app.AdvancedPanel.Title = 'Advanced';
            app.AdvancedPanel.Position = [420 1 181 120];

            % Create TargetDelayEditFieldLabel
            app.TargetDelayEditFieldLabel = uilabel(app.AdvancedPanel);
            app.TargetDelayEditFieldLabel.HorizontalAlignment = 'right';
            app.TargetDelayEditFieldLabel.Enable = 'off';
            app.TargetDelayEditFieldLabel.Position = [7 50 73 22];
            app.TargetDelayEditFieldLabel.Text = 'Target Delay';

            % Create TxphEQ
            app.TxphEQ = uieditfield(app.AdvancedPanel, 'numeric');
            app.TxphEQ.Limits = [0 Inf];
            app.TxphEQ.ValueChangedFcn = createCallbackFcn(app, @TxphEQValueChanged, true);
            app.TxphEQ.Enable = 'off';
            app.TxphEQ.Position = [95 50 62 22];

            % Create EnTxphEQ
            app.EnTxphEQ = uicheckbox(app.AdvancedPanel);
            app.EnTxphEQ.ValueChangedFcn = createCallbackFcn(app, @EnTxphEQValueChanged, true);
            app.EnTxphEQ.Enable = 'off';
            app.EnTxphEQ.Text = 'Phase Equalize';
            app.EnTxphEQ.Position = [41 73 105 22];

            % Create nsLabel
            app.nsLabel = uilabel(app.AdvancedPanel);
            app.nsLabel.Position = [160 50 25 22];
            app.nsLabel.Text = 'ns';

            % Create DesignResultsTX
            app.DesignResultsTX = uitable(app.TXFilterTab);
            app.DesignResultsTX.ColumnName = {'Design Results'};
            app.DesignResultsTX.RowName = {};
            app.DesignResultsTX.Position = [605 326 146 140];

            % Create DesignTXFilterButton
            app.DesignTXFilterButton = uibutton(app.TXFilterTab, 'push');
            app.DesignTXFilterButton.ButtonPushedFcn = createCallbackFcn(app, @DesignTXFilterButtonPushed, true);
            app.DesignTXFilterButton.BackgroundColor = [0 1 0];
            app.DesignTXFilterButton.Position = [606 305 144 20];
            app.DesignTXFilterButton.Text = 'Design TX Filter';

            % Create ExportDesignPanelTX
            app.ExportDesignPanelTX = uipanel(app.TXFilterTab);
            app.ExportDesignPanelTX.Title = 'Export Design';
            app.ExportDesignPanelTX.Position = [597 0 171 180];

            % Create ToWorkspaceButtonTX
            app.ToWorkspaceButtonTX = uibutton(app.ExportDesignPanelTX, 'push');
            app.ToWorkspaceButtonTX.Position = [35 128 100 22];
            app.ToWorkspaceButtonTX.Text = 'To Workspace';

            % Create ToftrFileButtonTX
            app.ToftrFileButtonTX = uibutton(app.ExportDesignPanelTX, 'push');
            app.ToftrFileButtonTX.ButtonPushedFcn = createCallbackFcn(app, @ToftrFileButton_2Pushed, true);
            app.ToftrFileButtonTX.Position = [35 98 100 22];
            app.ToftrFileButtonTX.Text = 'To ftr File';

            % Create ToNoOSFileButtonTX
            app.ToNoOSFileButtonTX = uibutton(app.ExportDesignPanelTX, 'push');
            app.ToNoOSFileButtonTX.Position = [35 68 100 22];
            app.ToNoOSFileButtonTX.Text = 'To No-OS File';

            % Create ProgramBoardButtonTX
            app.ProgramBoardButtonTX = uibutton(app.ExportDesignPanelTX, 'push');
            app.ProgramBoardButtonTX.Position = [35 38 100 22];
            app.ProgramBoardButtonTX.Text = 'Program Board';

            % Create IPAddressEditField_2Label
            app.IPAddressEditField_2Label = uilabel(app.ExportDesignPanelTX);
            app.IPAddressEditField_2Label.HorizontalAlignment = 'right';
            app.IPAddressEditField_2Label.Position = [5 8 68 22];
            app.IPAddressEditField_2Label.Text = 'IP Address:';

            % Create IPAddressEditField_2
            app.IPAddressEditField_2 = uieditfield(app.ExportDesignPanelTX, 'text');
            app.IPAddressEditField_2.HorizontalAlignment = 'right';
            app.IPAddressEditField_2.Position = [78 8 76 22];
            app.IPAddressEditField_2.Value = '192.168.3.2';

            % Create RXFilterTab
            app.RXFilterTab = uitab(app.TabGroup);
            app.RXFilterTab.Title = 'RX Filter';

            % Create UIAxesRX
            app.UIAxesRX = uiaxes(app.RXFilterTab);
            xlabel(app.UIAxesRX, 'Frequency (MHz)')
            ylabel(app.UIAxesRX, 'Magnitude (dB)')
            app.UIAxesRX.Box = 'on';
            app.UIAxesRX.XGrid = 'on';
            app.UIAxesRX.YGrid = 'on';
            app.UIAxesRX.Position = [3 146 598 327];

            % Create FilterSpecificationsPanel_2
            app.FilterSpecificationsPanel_2 = uipanel(app.RXFilterTab);
            app.FilterSpecificationsPanel_2.Title = 'Filter Specifications';
            app.FilterSpecificationsPanel_2.Position = [0 0 601 140];

            % Create FrequencyPanel_2
            app.FrequencyPanel_2 = uipanel(app.FilterSpecificationsPanel_2);
            app.FrequencyPanel_2.Title = 'Frequency';
            app.FrequencyPanel_2.Position = [210 1 210 120];

            % Create FpassEditField_2Label
            app.FpassEditField_2Label = uilabel(app.FrequencyPanel_2);
            app.FpassEditField_2Label.HorizontalAlignment = 'right';
            app.FpassEditField_2Label.Position = [59 73 38 22];
            app.FpassEditField_2Label.Text = 'Fpass';

            % Create RxFpass
            app.RxFpass = uieditfield(app.FrequencyPanel_2, 'numeric');
            app.RxFpass.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.RxFpass.Position = [112 73 60 22];

            % Create FstopEditField_2Label
            app.FstopEditField_2Label = uilabel(app.FrequencyPanel_2);
            app.FstopEditField_2Label.HorizontalAlignment = 'right';
            app.FstopEditField_2Label.Position = [61 47 36 22];
            app.FstopEditField_2Label.Text = 'Fstop';

            % Create RxFstop
            app.RxFstop = uieditfield(app.FrequencyPanel_2, 'numeric');
            app.RxFstop.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.RxFstop.Position = [112 47 60 22];

            % Create AnalogFcutoffEditField_2Label
            app.AnalogFcutoffEditField_2Label = uilabel(app.FrequencyPanel_2);
            app.AnalogFcutoffEditField_2Label.HorizontalAlignment = 'right';
            app.AnalogFcutoffEditField_2Label.Position = [14 20 83 22];
            app.AnalogFcutoffEditField_2Label.Text = 'Analog Fcutoff';

            % Create RxAnalogCuttoffFrequency
            app.RxAnalogCuttoffFrequency = uieditfield(app.FrequencyPanel_2, 'numeric');
            app.RxAnalogCuttoffFrequency.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.RxAnalogCuttoffFrequency.Position = [112 20 60 22];

            % Create MHzLabel_18
            app.MHzLabel_18 = uilabel(app.FrequencyPanel_2);
            app.MHzLabel_18.Position = [175 73 30 22];
            app.MHzLabel_18.Text = 'MHz';

            % Create MHzLabel_19
            app.MHzLabel_19 = uilabel(app.FrequencyPanel_2);
            app.MHzLabel_19.Position = [175 20 30 22];
            app.MHzLabel_19.Text = 'MHz';

            % Create MHzLabel_20
            app.MHzLabel_20 = uilabel(app.FrequencyPanel_2);
            app.MHzLabel_20.Position = [175 47 30 22];
            app.MHzLabel_20.Text = 'MHz';

            % Create MagnitudePanel_2
            app.MagnitudePanel_2 = uipanel(app.FilterSpecificationsPanel_2);
            app.MagnitudePanel_2.Title = 'Magnitude';
            app.MagnitudePanel_2.Position = [1 1 210 120];

            % Create ApassEditField_2Label
            app.ApassEditField_2Label = uilabel(app.MagnitudePanel_2);
            app.ApassEditField_2Label.HorizontalAlignment = 'right';
            app.ApassEditField_2Label.Position = [35 73 39 22];
            app.ApassEditField_2Label.Text = 'Apass';

            % Create RxApass
            app.RxApass = uieditfield(app.MagnitudePanel_2, 'numeric');
            app.RxApass.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.RxApass.Position = [89 73 60 22];

            % Create dBLabel_3
            app.dBLabel_3 = uilabel(app.MagnitudePanel_2);
            app.dBLabel_3.Position = [152 73 25 22];
            app.dBLabel_3.Text = 'dB';

            % Create dBLabel_4
            app.dBLabel_4 = uilabel(app.MagnitudePanel_2);
            app.dBLabel_4.Position = [152 48 25 22];
            app.dBLabel_4.Text = 'dB';

            % Create AstopEditField_2Label
            app.AstopEditField_2Label = uilabel(app.MagnitudePanel_2);
            app.AstopEditField_2Label.HorizontalAlignment = 'right';
            app.AstopEditField_2Label.Position = [38 48 36 22];
            app.AstopEditField_2Label.Text = 'Astop';

            % Create RxAstop
            app.RxAstop = uieditfield(app.MagnitudePanel_2, 'numeric');
            app.RxAstop.ValueChangedFcn = createCallbackFcn(app, @DesignConfigValueChanged, true);
            app.RxAstop.Position = [89 48 60 22];

            % Create AdvancedPanel_2
            app.AdvancedPanel_2 = uipanel(app.FilterSpecificationsPanel_2);
            app.AdvancedPanel_2.Title = 'Advanced';
            app.AdvancedPanel_2.Position = [420 1 180 120];

            % Create TargetDelayEditField_2Label
            app.TargetDelayEditField_2Label = uilabel(app.AdvancedPanel_2);
            app.TargetDelayEditField_2Label.HorizontalAlignment = 'right';
            app.TargetDelayEditField_2Label.Enable = 'off';
            app.TargetDelayEditField_2Label.Position = [7 50 73 22];
            app.TargetDelayEditField_2Label.Text = 'Target Delay';

            % Create TargetDelayEditField_2
            app.TargetDelayEditField_2 = uieditfield(app.AdvancedPanel_2, 'numeric');
            app.TargetDelayEditField_2.Enable = 'off';
            app.TargetDelayEditField_2.Position = [95 50 62 22];

            % Create PhaseEqualizeCheckBox_2
            app.PhaseEqualizeCheckBox_2 = uicheckbox(app.AdvancedPanel_2);
            app.PhaseEqualizeCheckBox_2.Enable = 'off';
            app.PhaseEqualizeCheckBox_2.Text = 'Phase Equalize';
            app.PhaseEqualizeCheckBox_2.Position = [41 73 105 22];

            % Create nsLabel_2
            app.nsLabel_2 = uilabel(app.AdvancedPanel_2);
            app.nsLabel_2.Position = [160 50 25 22];
            app.nsLabel_2.Text = 'ns';

            % Create DesignResultsRX
            app.DesignResultsRX = uitable(app.RXFilterTab);
            app.DesignResultsRX.ColumnName = {'Design Results'};
            app.DesignResultsRX.RowName = {};
            app.DesignResultsRX.Position = [605 326 146 140];

            % Create DesignRXFilterButton
            app.DesignRXFilterButton = uibutton(app.RXFilterTab, 'push');
            app.DesignRXFilterButton.ButtonPushedFcn = createCallbackFcn(app, @DesignRXFilterButtonPushed, true);
            app.DesignRXFilterButton.BackgroundColor = [0 1 0];
            app.DesignRXFilterButton.Position = [606 303 144 22];
            app.DesignRXFilterButton.Text = 'Design RX Filter';

            % Create ExportDesignPanelRX
            app.ExportDesignPanelRX = uipanel(app.RXFilterTab);
            app.ExportDesignPanelRX.Title = 'Export Design';
            app.ExportDesignPanelRX.Position = [597 0 171 180];

            % Create ToWorkspaceButtonRX
            app.ToWorkspaceButtonRX = uibutton(app.ExportDesignPanelRX, 'push');
            app.ToWorkspaceButtonRX.Position = [35 128 100 22];
            app.ToWorkspaceButtonRX.Text = 'To Workspace';

            % Create ToftrFileButtonRX
            app.ToftrFileButtonRX = uibutton(app.ExportDesignPanelRX, 'push');
            app.ToftrFileButtonRX.ButtonPushedFcn = createCallbackFcn(app, @ToftrFileButton_2Pushed, true);
            app.ToftrFileButtonRX.Position = [35 98 100 22];
            app.ToftrFileButtonRX.Text = 'To ftr File';

            % Create ToNoOSFileButtonRX
            app.ToNoOSFileButtonRX = uibutton(app.ExportDesignPanelRX, 'push');
            app.ToNoOSFileButtonRX.Position = [35 68 100 22];
            app.ToNoOSFileButtonRX.Text = 'To No-OS File';

            % Create ProgramBoardButtonRX
            app.ProgramBoardButtonRX = uibutton(app.ExportDesignPanelRX, 'push');
            app.ProgramBoardButtonRX.Position = [35 38 100 22];
            app.ProgramBoardButtonRX.Text = 'Program Board';

            % Create IPAddressEditFieldLabel
            app.IPAddressEditFieldLabel = uilabel(app.ExportDesignPanelRX);
            app.IPAddressEditFieldLabel.HorizontalAlignment = 'right';
            app.IPAddressEditFieldLabel.Position = [5 8 68 22];
            app.IPAddressEditFieldLabel.Text = 'IP Address:';

            % Create IPAddressEditField
            app.IPAddressEditField = uieditfield(app.ExportDesignPanelRX, 'text');
            app.IPAddressEditField.HorizontalAlignment = 'right';
            app.IPAddressEditField.Position = [78 8 76 22];
            app.IPAddressEditField.Value = '192.168.3.2';

            % Create MHzLabel_16
            app.MHzLabel_16 = uilabel(app.ad936xFilterWizardUIFigure);
            app.MHzLabel_16.Position = [127 401 44 22];
            app.MHzLabel_16.Text = 'MHz';

            % Create MHzLabel_17
            app.MHzLabel_17 = uilabel(app.ad936xFilterWizardUIFigure);
            app.MHzLabel_17.Position = [291 401 44 22];
            app.MHzLabel_17.Text = 'MHz';

            % Create Button
            app.Button = uibutton(app.ad936xFilterWizardUIFigure, 'push');
            app.Button.Icon = 'Analog_Devices_Logo.png';
            app.Button.Position = [11 13 320 110];
            app.Button.Text = '';
        end
    end

    methods (Access = public)

        % Construct app
        function app = DesignerGUI_exported(varargin)

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ad936xFilterWizardUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ad936xFilterWizardUIFigure)
        end
    end
end
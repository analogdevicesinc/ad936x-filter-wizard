classdef FilterDesignerTests < matlab.unittest.TestCase
    % Test ad936x Filter Designer
    %
    % This testing utilizes the provided 'ad9361_settings.mat' file to
    % generate test input vectors which are permuted to exercise additional
    % functionality of the designer.
    %
    % Currently these tests compare results of the generated (MEX/DLL) with
    % the existing implementation of the filter designer. The following
    % fields are compared:
    % - firtaps
    %
    % The DLL generated code is tested through an opaque call from a MEX
    % function allowing direct access to the DLL outputs in MATLAB
    %
    % Example call:
    %  test = FilterDesignerTests;
    %  test.run()
    
    properties
        args = '{input.Rdata, input.Fpass, input.Fstop, input.caldiv, input.FIR, input.HB1, input.PLL_mult, input.Apass, input.Astop, input.phEQ, input.HB2, input.HB3, input.Type, input.RxTx, input.RFbw, input.DAC_div, input.converter_rate, input.PLL_rate, input.Fcenter, input.wnom, input.FIRdBmin, input.int_FIR}';
        functionName = 'internal_design_filter_cg';
        passedCodegenMEX = false;
        passedCodegenDLL = false;
        settingsLoaded = false;
        ad9361_settings = [];
        TargetHWDeviceType = [];
    end
    
    methods(TestClassSetup)
        % Load example settings from mat file
        function loadTemplateSettings(testCase)
            a = load('ad9361_settings.mat');
            testCase.ad9361_settings = a.ad9361_settings;
            testCase.settingsLoaded = true;
        end
        % Set codegen configuration
        function setCodegenConfig(testCase)
           if ismac
               testCase.TargetHWDeviceType = 'Intel->x86-64 (Mac OS X)';
           elseif ispc
               testCase.TargetHWDeviceType = 'Intel->x86-64 (Windows64)';
           else
               testCase.TargetHWDeviceType = 'Intel->x86-64 (Linux 64)';
           end
        end
        % Test Codegen of MEX Target
        function testCodegenBuildMEX(testCase)
            %% Setup workspace
            % Get example struct
            inputVar = testCase.ad9361_settings.tx.LTE5;
            % Fill out necessary fields
            input = process_input(inputVar); %#ok<NASGU>
            %% Call codegen
            cfg = coder.config('mex');
            cfg.TargetLang='C++';
            result = codegen('-config','cfg',testCase.functionName,'-args',testCase.args);
            testCase.passedCodegenMEX = result.summary.passed;
            testCase.verifyTrue(result.summary.passed);
        end
        % Test Codegen of DLL Target
        function testCodegenBuildDLL(testCase)
            %% Setup workspace
            % Get example struct
            inputVar = testCase.ad9361_settings.tx.LTE5;
            % Fill out necessary fields
            input = process_input(inputVar); %#ok<NASGU>
            %% Call codegen
            cfg = coder.config('lib');
            cfg.TargetLang='C++';
            cfg.FilePartitionMethod='SingleFile';
            cfg.HardwareImplementation.TargetHWDeviceType = testCase.TargetHWDeviceType;
            result = codegen('-config','cfg',testCase.functionName,'-O ','disable:openmp','-args',testCase.args);
            testCase.passedCodegenDLL = result.summary.passed;
            testCase.verifyTrue(result.summary.passed);
            if testCase.passedCodegenDLL && ismac % Move library to root (rpath is buggy on mac)
                !mv codegen/lib/internal_design_filter_cg/internal_design_filter_cg.a .
            end
        end
        % Create MAT file for generated tests and generate tests
        function GenTestMATFiles(testCase)
            if ~testCase.passedCodegenDLL || ~testCase.settingsLoaded
                error('Must generate code first and load settings');
            end
            filename = 'ad9361_settings_processed_test';
            
            %% RX
            in = testCase.ad9361_settings.rx.LTE5;
            % Generate template expected results
            input = testCase.input_cooker(in);
            refResult = internal_design_filter(input); % reference
            firtaps = refResult.firtaps; %#ok<NASGU>
            % Save to file
            save(filename, 'firtaps', 'input');
            % Generate TestCase
            cfg = coder.config('mex');
            cfg.TargetLang='C++';
            if ismac
                cfg.CustomLibrary = [testCase.functionName,'.a'];
            elseif isunix
                cfg.CustomLibrary = [testCase.functionName,'.a'];
            else
                cfg.CustomLibrary = [testCase.functionName,'.lib'];
            end
            additionalSource = {[testCase.functionName,'.h']};
            cfg.CustomInclude = ['codegen/lib/',testCase.functionName,'/'];
            result = codegen('-config','cfg','TestToBeGenerated',...
                additionalSource{:},'-o','TestToBeGenerated_rx_mex');
            testCase.verifyTrue(result.summary.passed);
            %% TX
            in = testCase.ad9361_settings.tx.LTE5;
            % Generate template expected results
            input = testCase.input_cooker(in);
            refResult = internal_design_filter(input); % reference
            firtaps = refResult.firtaps; %#ok<NASGU>
            % Save to file
            save(filename, 'firtaps', 'input');
            % Generate TestCase
            cfg = coder.config('mex');
            cfg.TargetLang='C++';
            if ismac
                cfg.CustomLibrary = [testCase.functionName,'.a'];
            elseif isunix
                cfg.CustomLibrary = [testCase.functionName,'.a'];
            else
                cfg.CustomLibrary = [testCase.functionName,'.lib'];
            end
            additionalSource = {[testCase.functionName,'.h']};
            cfg.CustomInclude = ['codegen/lib/',testCase.functionName,'/'];
            result = codegen('-config','cfg','TestToBeGenerated',...
                additionalSource{:},'-o','TestToBeGenerated_tx_mex');
            testCase.verifyTrue(result.summary.passed);
        end
        
        
    end
    
    methods(TestClassTeardown)
        function removeCodegenFiles(testCase)
            % Remove generated code
            if testCase.passedCodegenDLL && testCase.passedCodegenMEX
                [~,~,~] = rmdir('codegen','s');
            end
            if ismac
                delete([testCase.functionName,'.a']);
            end
            delete('TestToBeGenerated_tx_mex.*')
            delete('TestToBeGenerated_rx_mex.*')
            delete('ad9361_settings_processed_test.mat')
        end
    end
    
    methods (Static)
        
        % Build input so all fields are filled
        function input = input_cooker(input)
            % support a simple data rate input otherwise it must be a structure
            if isfloat(input)
                input = struct('Rdata', input);
            end
            input = cook_input(input);
            
            % use the internal FIR if unspecified
            if ~isfield(input, 'int_FIR')
                input.int_FIR = 1;
            end
            
            % nominal frequency can't be zero
            if ~input.wnom
                input.wnom = (input.PLL_rate/input.caldiv)*(log(2)/(2*pi));
            end
        end
        
        % Modify input based on additional configuration
        function input = modifyInput(input,config)
            cFields = fields(config);
            for field = 1:length(cFields)
                if strcmp(cFields{field},'txrx')% This field must be set by the user first
                    continue;
                end
                input = setfield(input, cFields{field}, getfield(config,cFields{field})); %#ok<SFLD,GFLD>
            end
        end
        
    end
    
    methods % Non-Static Test Scaffolding
        
        % This will test the mexed designer with int_FIR=1
        function testFunctionGeneral(testCase,config)
            if ~testCase.passedCodegenMEX || ~testCase.settingsLoaded
                error('Must generate code first and load settings');
            end
            % Get settings
            if strcmp(config.txrx,'tx')
                txrx = testCase.ad9361_settings.tx;
            else
                txrx = testCase.ad9361_settings.rx;
            end
            frt = fields(txrx);
            % Test all configurations LTE5-20
            for s = 1:length(fields(txrx))
                % Build input
                str = char(frt{s});
                in = getfield(txrx,str); %#ok<GFLD>
                input = testCase.input_cooker(in);
                % Update settings based on config
                input = testCase.modifyInput(input,config);
                % Test
                cgResultFirtaps = call_filter_designer_cg(input,true); % codegen mex
                refResult = internal_design_filter(input); % reference
                % Evaluate errors
                testCase.verifyEqual(double(cgResultFirtaps),refResult.firtaps,'AbsTol',2);
            end
            
        end
        % This will test the mexed designer with int_FIR=0
        function testFunctionGeneralLengthCheck(testCase,config)
            if ~testCase.passedCodegenMEX || ~testCase.settingsLoaded
                error('Must generate code first and load settings');
            end
            % Get settings
            if strcmp(config.txrx,'tx')
                txrx = testCase.ad9361_settings.tx;
            else
                txrx = testCase.ad9361_settings.rx;
            end
            frt = fields(txrx);
            % Test all configurations LTE5-20
            for s = 1:length(fields(txrx))
                % Build input
                str = char(frt{s});
                in = getfield(txrx,str); %#ok<GFLD>
                input = testCase.input_cooker(in);
                % Update settings based on config
                input = testCase.modifyInput(input,config);
                % Test
                cgResultFirtaps = call_filter_designer_cg(input,true); % codegen mex
                refResult = internal_design_filter(input); % reference
                % Evaluate errors
                testCase.verifyEqual(length(cgResultFirtaps),length(refResult.firtaps));
                testCase.verifyEqual(double(cgResultFirtaps),double(int16(refResult.firtaps)),'AbsTol',2);
            end
            
        end
        % This will test the generated DLL file with a mex wrapper
        function testGeneratedFunctionGeneral(testCase,config)
            if ~testCase.passedCodegenMEX || ~testCase.settingsLoaded
                error('Must generate code first and load settings');
            end
            % Get settings
            if strcmp(config.txrx,'tx')
                txrx = testCase.ad9361_settings.tx;
            else
                txrx = testCase.ad9361_settings.rx;
            end
            frt = fields(txrx);
            % Test all configurations LTE5-20
            for s = 1:length(fields(txrx))
                % Build input
                str = char(frt{s});
                in = getfield(txrx,str); %#ok<GFLD>
                input = testCase.input_cooker(in);
                % Update settings based on config
                input = testCase.modifyInput(input,config);
                % Save test data to file
                filename = 'ad9361_settings_processed_test';
                refResult = internal_design_filter(input); % reference
                firtaps = refResult.firtaps; %#ok<NASGU>
                % Save to file
                save(filename, 'firtaps', 'input');
                % Test
                if strcmp(config.txrx,'tx')
                    firtaps = TestToBeGenerated_tx_mex();
                else
                    firtaps = TestToBeGenerated_rx_mex();
                end
                testCase.verifyEqual(firtaps,int16(refResult.firtaps));
            end
            
        end
        % This will test the generated DLL file with a mex wrapper with
        % int_FIR=0
        function testGeneratedFunctionGeneralLengthCheck(testCase,config)
            if ~testCase.passedCodegenMEX || ~testCase.settingsLoaded
                error('Must generate code first and load settings');
            end
            % Get settings
            if strcmp(config.txrx,'tx')
                txrx = testCase.ad9361_settings.tx;
            else
                txrx = testCase.ad9361_settings.rx;
            end
            frt = fields(txrx);
            % Test all configurations LTE5-20
            for s = 1:length(fields(txrx))
                % Build input
                str = char(frt{s});
                in = getfield(txrx,str); %#ok<GFLD>
                input = testCase.input_cooker(in);
                % Update settings based on config
                input = testCase.modifyInput(input,config);
                % Save test data to file
                filename = 'ad9361_settings_processed_test';
                refResult = internal_design_filter(input); % reference
                firtaps = refResult.firtaps; %#ok<NASGU>
                % Save to file
                save(filename, 'firtaps', 'input');
                % Test
                if strcmp(config.txrx,'tx')
                    firtaps = TestToBeGenerated_tx_mex();
                else
                    firtaps = TestToBeGenerated_rx_mex();
                end
                % Evaluate errors
                testCase.verifyEqual(length(firtaps),length(refResult.firtaps));
                testCase.verifyEqual(double(firtaps),double(int16(refResult.firtaps)),'AbsTol',2);
            end
            
        end
        
        
    end
    
    methods (Test)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Tests
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%% RX
        
        % Test MEX results with standard setting
        function testRXMEX(testCase)
            config = struct;
            config.txrx = 'rx';
            testCase.testFunctionGeneral(config);
        end
        % Test MEX results with reference without EQ
        function testRXEQMEX(testCase)
            config = struct;
            config.txrx = 'rx';
            config.phEQ = 1;
            testCase.testFunctionGeneral(config);
        end
        % Test MEX results with reference without EQ
        function testRXNonStandardFIRMEX(testCase)
            config = struct;
            config.txrx = 'rx';
            config.int_FIR = 0;
            testCase.testFunctionGeneralLengthCheck(config);
        end
        % Test MEX results with reference without EQ
        function testRXGenericMEX(testCase)
            config = cook_input(struct('Rdata',61.44e6,'RxTx','Rx'));
            config.txrx = 'rx';
            testCase.testFunctionGeneral(config);
        end
        % Test MEX results with reference without EQ at many different
        % rates
        function testRXManyRatesMEX(testCase)
            rates = [1,10:10:60].*1e6;
            for rate = rates
                config = cook_input(struct('Rdata',rate,'RxTx','Rx'));
                config.txrx = 'rx';
                testCase.testFunctionGeneral(config);
            end
        end
        % Test DLL results with standard setting
        function testRXDLL(testCase)
            config = struct;
            config.txrx = 'rx';
            testCase.testGeneratedFunctionGeneral(config);
        end
        % Test DLL results with best length FIR
        function testRXNonStandardFIRDLL(testCase)
            config = struct;
            config.txrx = 'rx';
            config.int_FIR = 0;
            testCase.testGeneratedFunctionGeneralLengthCheck(config);
        end

        %%%% TX
        
        % Test MEX results with reference without EQ
        function testTXMEX(testCase)
            config = struct;
            config.txrx = 'tx';
            testCase.testFunctionGeneral(config);
        end
        % Test MEX results with reference without EQ
        function testTXEQMEX(testCase)
            config = struct;
            config.txrx = 'tx';
            config.phEQ = 1;
            testCase.testFunctionGeneral(config);
        end
        % Test MEX results with reference without EQ
        function testTXNonStandardFIRMEX(testCase)
            config = struct;
            config.txrx = 'tx';
            config.int_FIR = 0;
            testCase.testFunctionGeneralLengthCheck(config);
        end
        % Test MEX results with reference without EQ
        function testTXGenericMEX(testCase)
            config = cook_input(struct('Rdata',61.44e6,'RxTx','Tx'));
            config.txrx = 'tx';
            testCase.testFunctionGeneral(config);
        end
        % Test MEX results with reference without EQ at many different
        % rates
        function testTXManyRatesMEX(testCase)
            rates = [1,10:10:60].*1e6;
            for rate = rates
                config = cook_input(struct('Rdata',rate,'RxTx','Tx'));
                config.txrx = 'tx';
                testCase.testFunctionGeneral(config);
            end
        end
        % Test DLL results with standard setting
        function testTXDLL(testCase)
            config = struct;
            config.txrx = 'tx';
            testCase.testGeneratedFunctionGeneral(config);
        end
        % Test DLL results with best length FIR
        function testTXNonStandardFIRDLL(testCase)
            config = struct;
            config.txrx = 'tx';
            config.int_FIR = 0;
            testCase.testGeneratedFunctionGeneralLengthCheck(config);
        end
    end
    
end
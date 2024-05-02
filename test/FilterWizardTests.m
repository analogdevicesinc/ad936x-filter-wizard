classdef FilterWizardTests < matlab.unittest.TestCase
    
    properties
        SampleRates = [0.52083333,0.6:0.1:61.44].*1e6;
        LTEModes = {'LTE5','LTE10','LTE15','LTE20'};
        MaxRippleDB = 1;
    end
    
    methods(TestClassSetup)
    end
    
    methods (Static)
    end
    
    methods (Test)
        
        function testAutoGenerationRipple(testCase)
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            sr = testCase.SampleRates;
            limit = testCase.MaxRippleDB;
            % Test ripple of generated filters           
            % parfor r = 1:length(sr)
            %     out = internal_design_filter_opt_ripple(sr(r));
            %     verifyThat(testCase, out.Apass_actual, IsLessThanOrEqualTo(limit), ...
            %         sprintf('Generated filter for rate %d with ripple %f (Limit %f)',...
            %         sr(r),out.Apass_actual,limit))
            % end
        end
        
        function testLTEFilterGeneration(testCase)
            import matlab.unittest.constraints.IsEqualTo
            % Verify LTE filters are created correctly
            for m = testCase.LTEModes
                i = load([m{:},'_input.mat']);
                o = load([m{:},'_output.mat']);
                output = internal_design_filter(i.input);
                testCase.verifyThat(o.output.firtaps,IsEqualTo(output.firtaps));
            end
        end
        
        function testOFFilterGeneration(testCase)
            import matlab.unittest.constraints.IsEqualTo
            % Verify EZ OR filter
            i = load('ez555448.mat');
            output = internal_design_filter(i.input);
            testCase.verifyThat(...
                double(int16(output.firtaps)),...
                IsEqualTo(...
                double(int32(output.firtaps))));
        end
        
    end
end
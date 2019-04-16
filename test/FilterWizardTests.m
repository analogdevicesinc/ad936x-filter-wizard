classdef FilterWizardTests < matlab.unittest.TestCase
    
    properties
        SampleRates = [0.52083333,0.6:0.1:61.44].*1e6;
    end
    
    methods(TestClassSetup)
    end
    
    methods (Static)
    end
    
    methods (Test)
        
        function testAutoGenerationRipple(testCase)
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            sr = testCase.SampleRates;
            % Test ripple of generated filters           
            parfor r = 1:length(sr)
                out = internal_design_filter_opt_ripple(sr(r));
                verifyThat(testCase, out.Apass_actual, IsLessThanOrEqualTo(out.Apass), ...
                    sprintf('Generated filter for rate %d failed to meet ripple target',sr(r)))
            end
        end
        
    end
end
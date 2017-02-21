import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.TestRunProgressPlugin

% Build up test suite
suite = TestSuite.fromClass(?FilterDesignerTests);
% Outline runner
runner = TestRunner.withNoPlugins;
% Add progress verbosity
p = TestRunProgressPlugin.withVerbosity(4);
runner.addPlugin(p);
% Run and display info
results = runner.run(suite);
disp(results.table);

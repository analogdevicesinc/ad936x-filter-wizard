function [found,h] = testCodegen(rate)


fd = ad936xFilterDesigner;
fd.DataRate = rate;
fd.AutoSetRates();
[found,h] = fd.designFilter('Tx');

end

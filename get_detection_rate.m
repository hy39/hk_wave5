% get_detection_rate
% To calcualate the PCR detection rate
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk) 
function [ det_rate ] = get_detection_rate(IFD50)

Iprev = 1:5000000;
det_rate = 0.7*exp(-log(2).*(Iprev./IFD50))+0.1;

end

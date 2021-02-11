
function oe = mee2oe(mee)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION THAT CONVERTS FROM THE SET OF MODIFIED EQUINOCIAL ELEMENTS TO 
%  THE SET OF CLASSICAL ORBITAL ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%  mee(1) = semilatus rectum of orbit (adimen)
%  mee(2) = f equinoctial element
%  mee(3) = g equinoctial element
%  mee(4) = h equinoctial element
%  mee(5) = k equinoctial element
%  mee(6) = true longitude (radians)
%
% OUPUT
% oe(1) = semimajor-axis
% oe(2) = eccentricity
% oe(3) = inclination
% oe(4) = RAN
% oe(5) = argument of the perigee
% oe(6) = true anomaly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a  = mee(:,1)./(1-mee(:,2).^2-mee(:,3).^2);
e  = sqrt(mee(:,2).^2 + mee(:,3).^2);
i  = atan2(2*sqrt(mee(:,4).^2+mee(:,5).^2), 1-mee(:,4).^2-mee(:,5).^2);
Om = atan2(mee(:,5),mee(:,4));
w  = atan2(mee(:,3).*mee(:,4)-mee(:,2).*mee(:,5),mee(:,2).*mee(:,4) + mee(:,3).*mee(:,5));
theta = mee(:,6) - Om -w;

oe = [a e i Om w theta];
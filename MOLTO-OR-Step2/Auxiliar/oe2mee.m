function mee = oe2mee(oe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION THAT CONVERTS FROM THE SET OF CLASSICAL ORBITAL ELEMENTS TO 
%  THE SET OF MODIFIED EQUINOCCIAL ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% oe(1) = semimajor-axis
% oe(2) = eccentricity
% oe(3) = inclination
% oe(4) = RAN
% oe(5) = argument of the perigee
% oe(6) = true anomaly
% OUTPUT
%  mee(1) = semilatus rectum of orbit (adimen)
%  mee(2) = f equinoctial element
%  mee(3) = g equinoctial element
%  mee(4) = h equinoctial element
%  mee(5) = k equinoctial element
%  mee(6) = true longitude (radians)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = oe(:,1).*(1-oe(:,2).^2);
f = oe(:,2).*cos(oe(:,5) + oe(:,4));
g = oe(:,2).*sin(oe(:,5) + oe(:,4));
h = tan(oe(:,3)/2).*cos(oe(:,4));
k = tan(oe(:,3)/2).*sin(oe(:,4));
L = oe(:,4) + oe(:,5) + oe(:,6);

mee = [p f g h k L];
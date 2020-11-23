function eps_pathlength = ioi_epsilon_pathlength(whichCurve,...
    baseline_hbt, baseline_hbo, baseline_hbr, debug)
%	This function estimates epsilon * D, it takes into account the camera
%	response, the leds spectra and uses a pathlength factor either set from Kohl
%	or Dunn in the literature.
%   This module is dependent on this file which contains all hardware info for
%   the setup, needs to specify the leds and the camera response. we are still a
%   bit dependent on the RGY but we could abstract this (however the lambdas
%   would need to be registered to specific hardware still
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

load('SysSpect.mat')
Optics =  ReadConfigFile();

% Rough baseline concentrations (in uM) : 100 uM (in the brain)
c_tot = baseline_hbt*1e-6; %100e-6;

lambda_vec= linspace(400, 700, 301);

%Camera spectrum:
switch Optics.Camera
    case 'PF1024'
        c_camera = PF1024;
    case 'PF1312'
        c_camera = PF1312;
    case 'BFly' 
        c_camera = BFly;
    case 'CS2100M'
        c_camera = ThorQLux;
    otherwise
        c_camera = PF1024;
end

% Red spectrum:
switch Optics.RIllum
    case 'Red'
        c_led(1,:) = Red;
    otherwise 
        c_led(1,:) = zeros(1,301);
end
switch Optics.RFilter_illum
    case {'FF01496LP', 'FF0161850'}
        eval(['c_led(1,:) = c_led(1,:).*' Optics.RFilter_illum ';']);
    otherwise 
        c_led(1,:) = c_led(1,:);
end
switch Optics.RFilter_detec
    case {'FF0161850', 'FF01496LP'}
        eval(['c_led(1,:) = c_led(1,:).*' Optics.RFilter_detec ';']);
    otherwise 
        c_led(1,:) = c_led(1,:);
end

% Green spectrum:
switch Optics.GIllum
    case 'Green'
        c_led(2,:) = Green;
    otherwise 
        c_led(2,:) = zeros(1,301);
end
switch Optics.GFilter_illum
    case {'FF01496LP'}
        eval(['c_led(2,:) = c_led(2,:).*' Optics.GFilter_illum ';']);
    otherwise 
        c_led(2,:) = c_led(2,:);
end
switch Optics.GFilter_detec
    case {'FF0151444', 'FF01496LP'}
        eval(['c_led(2,:) = c_led(2,:).*' Optics.GFilter_detec ';']);
    otherwise 
        c_led(2,:) = c_led(2,:);
end

% Yellow spectrum:
switch Optics.YIllum
    case 'Yellow'
        c_led(3,:) = Yellow;
    otherwise 
        c_led(3,:) = zeros(1,301);
end
switch Optics.YFilter_illum
    case {'FF01496LP', 'FF0161850'}
        eval(['c_led(3,:) = c_led(3,:).*' Optics.YFilter_illum ';']);
    otherwise 
        c_led(3,:) = c_led(3,:);
end
switch Optics.YFilter_detec
    case {'FF0161850', 'FF01496LP'}
        eval(['c_led(3,:) = c_led(3,:).*' Optics.YFilter_detec ';']);
    otherwise 
        c_led(3,:) = c_led(3,:);
end

c_pathlength = ioi_path_length_factor(400, 700, 301, c_tot*1000, whichCurve);
[c_ext_hbo,c_ext_hbr] = ioi_get_extinctions(400, 700, 301);

if( debug==1 )
    figure;
    subplot(2,2,1)
    plot(lambda_vec,c_led(1,:),'r')
    hold on
    plot(lambda_vec,c_led(2,:),'g')
    hold on
    plot(lambda_vec,c_led(3,:),'y')
end

%Input basehbt1 not used...

% Create vectors of values for the fits
CHbO = baseline_hbo/baseline_hbt*c_tot*(.5:.1:1.5);
CHbR = baseline_hbr/baseline_hbt*c_tot*(.5:.1:1.5);

% In this computation below we neglect the fact that pathlength changes
% with total concentration (it is fixed for a Ctot of 100e-6)
eps_pathlength=zeros(3,2);

% Preallocating memory for speed increase //EGC
IHbO = zeros(size(CHbO));
IHbR = zeros(size(CHbR));

for iled=1:3
    for iconc = 1:length(CHbO)
        IHbO(iconc) = sum(c_camera .* c_led(iled,:) .* exp(-c_ext_hbo .* c_pathlength * CHbO(iconc)),2) ; %	Measured intensity for different concentrations
        IHbR(iconc) = sum(c_camera .* c_led(iled,:) .* exp(-c_ext_hbr .* c_pathlength * CHbR(iconc)),2) ;
    end
    IHbO = IHbO/median(IHbO);
    IHbR = IHbR/median(IHbR);

    % Compute effective eps
    p1 = polyfit(CHbO,-log(IHbO),1);
    p2 = polyfit(CHbR,-log(IHbR),1);
    HbRL = p2(1); %epsilon*D HbR effectif
    HbOL = p1(1);%epsilon*D HbO effectif
    eps_pathlength(iled,1)=HbOL;
    eps_pathlength(iled,2)=HbRL;
end



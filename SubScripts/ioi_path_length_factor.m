function pathlength = ioi_path_length_factor(L1, L2, N, cHBT, whichCurve)
% Compute the pathlength factor in OD x mm(-1)

% Values taken from: 
%   Kohl, Matthias,Ute Lindauer, Georg Royl, Marc K�uhl, Lorenz Gold,
%   Arno Villringer and Ulrich Dirnagl Physical model for the spectroscopic analysis of cortical
%   intrinsic optical signals,Phys. Med. Biol. 45 (2000) 3749�3764
%   figure 8a) S=50%, Tot hb 0.02
%   Values between 400-470 are set to the same value as 480.
%   Values between 690-800 are set ot the same value as 680

%  de la figure 8a) est ecrit en mm-1, mais les vraies unit�s sont des mm
%   Lambda    differential path length
%   colonne 2 OD*mm(converti par la suite en cm)

%   L1 Start wavelength
%   L2 End wavelength
%   N  Number of points
%   cHBT : Total hemoglobin concentration (default 0.02mM)
%   whichCurve: 'Kohl' or 'Dunn'
%   deltaL: Step by half a step to evaluate values at mid point

Ext.Kohl=   ...
   [350	0.2     0.2 % longueur d'onde, sat=50% , sat=70%
    400	0.2     0.2
    410	0.2     0.2
    420	0.2     0.2
    430	0.2     0.2
    440	0.2     0.2
    450	0.2     0.2
    460	0.2     0.2
    470	0.2     0.2
    480	0.2     0.2
    490	0.2     0.2
    500	0.2     0.2
    510	0.176	0.176
    520	0.161	0.161
    530	0.125	0.125
    540	0.11	0.11
    550	0.11	0.11
    560	0.118	0.118
    570	0.118	0.118
    580	0.118	0.12272
    590	0.2     0.216
    600	0.363	0.40656
    610	0.5     0.58
    620	0.663	0.7956
    630	0.76	0.912
    640	0.89	1.068
    650	0.94	1.128
    660	1.05	1.26
    670	1.18	1.416
    680	1.25	1.5
    690	1.25	1.5
    700	1.25	1.5
    710	1.25	1.5
    720	1.25	1.5
    730	1.25	1.5
    740	1.25	1.5
    750	1.25	1.5
    760	1.25	1.5
    770	1.25	1.5
    780	1.25	1.5
    790	1.25	1.5
    820	1.25	1.5
    ];

Ext.Dunn=[ 450 0.52     %extrapol� Simon
        560 0.52    %Donn�es de Dunn2005
        570 0.5     %Donn�es de Dunn2005
        580 0.5     %Donn�es de Dunn2005
        590 .93     %Donn�es de Dunn2005
        600 1.58    %Donn�es de Dunn2005
        610 2.05    %Donn�es de Dunn2005
        620 2.7     %extrapol� � l'aide de Kohl
        630 3.1     %extrapol� � l'aide de Kohl
        640 3.64    %extrapol� � l'aide de Kohl
        650 3.85    %extrapol� � l'aide de Kohl
        750 3.85];  %extrapol� Simon

if strcmp(whichCurve,'Kohl')
    E=Ext.Kohl;
    E(:,2:end)=E(:,2:end)/10; % In centimeters
elseif strcmp(whichCurve,'Dunn')
    E=Ext.Dunn;
    E(:,2:end)=E(:,2:end)/10; % In centimeters
else
    disp('ioi_path_length_factor: Error, no curve specified, taking data from Dunn')
    E=Ext.Dunn;
    E(:,2:end)=E(:,2:end)/10; % In centimeters
end

% Add boundaries if start and end wavelengths are outside our reach
if E(1,1)>L1
    E=[ L1 0 ;
        E(1,1)*.9999 0;
        E(:,:) ];
end

if E(end,1)<L2
    E=[ E(:,:);
        E(end,1)*1.00001 0;
        L2 0 ];
end

xi = linspace(L1,L2,N);
x=E(:,1); % avant x = linspace(E(1,1),E(end,1),size(E,1)); %remplac� fev 09
pathlength = E(:,2);
pathlength = interp1(x,pathlength,xi);

% trouve L pour des concentration diff�rentes de 0.02 (� partir de la
% figure 8 de kohl). On multiplie la courbe � 0.02 par un facteur de
% correction   
x=[0.01 0.02 0.04 0.06 0.08 0.10 0.14 0.18]; %concentration d'h�moglobine
y=[  .33 .21 .115 .09 .075 0.062 .05 .04 ]; % valeur Da=[OD*mm] � 480 nm
correction=interp1(x,y,cHBT)/interp1(x,y,.02);
pathlength=pathlength*correction;

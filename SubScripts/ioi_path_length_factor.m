function pathlength = ioi_path_length_factor(L1, L2, N, cHBT, whichCurve)
% Compute the pathlength factor in OD x mm
%
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
    
Ext.Silico=[  400.0000    0.0899;  401.0000    0.0919;  402.0000    0.0865;
              403.0000    0.0826;  404.0000    0.0771;  405.0000    0.0749;
              406.0000    0.0695;  407.0000    0.0641;  408.0000    0.0577;
              409.0000    0.0588;  410.0000    0.0544;  411.0000    0.0532;
              412.0000    0.0503;  413.0000    0.0507;  414.0000    0.0480;
              415.0000    0.0494;  416.0000    0.0485;  417.0000    0.0481;
              418.0000    0.0469;  419.0000    0.0495;  420.0000    0.0489;
              421.0000    0.0493;  422.0000    0.0525;  423.0000    0.0514;
              424.0000    0.0534;  425.0000    0.0544;  426.0000    0.0569;
              427.0000    0.0591;  428.0000    0.0602;  429.0000    0.0594;
              430.0000    0.0599;  431.0000    0.0585;  432.0000    0.0620;
              433.0000    0.0649;  434.0000    0.0683;  435.0000    0.0735;  
              436.0000    0.0743;  437.0000    0.0777;  438.0000    0.0840;
              439.0000    0.0920;  440.0000    0.1011;  441.0000    0.1082;
              442.0000    0.1121;  443.0000    0.1276;  444.0000    0.1393;
              445.0000    0.1501;  446.0000    0.1602;  447.0000    0.1853;
              448.0000    0.2110;  449.0000    0.2396;  450.0000    0.2755;
              451.0000    0.3112;  452.0000    0.3479;  453.0000    0.3924;
              454.0000    0.4300;  455.0000    0.4517;  456.0000    0.4688;
              457.0000    0.4759;  458.0000    0.4951;  459.0000    0.5040;
              460.0000    0.5294;  461.0000    0.5415;  462.0000    0.5608;
              463.0000    0.5653;  464.0000    0.5780;  465.0000    0.5905;
              466.0000    0.6049;  467.0000    0.6231;  468.0000    0.6315;
              469.0000    0.6514;  470.0000    0.6653;  471.0000    0.6664;
              472.0000    0.6796;  473.0000    0.6895;  474.0000    0.7047;
              475.0000    0.7116;  476.0000    0.7105;  477.0000    0.7282;
              478.0000    0.7271;  479.0000    0.7448;  480.0000    0.7498;
              481.0000    0.7568;  482.0000    0.7621;  483.0000    0.7658;
              484.0000    0.7623;  485.0000    0.7745;  486.0000    0.7714;
              487.0000    0.7704;  488.0000    0.7776;  489.0000    0.7769;
              490.0000    0.7731;  491.0000    0.7795;  492.0000    0.7815;
              493.0000    0.7780;  494.0000    0.7856;  495.0000    0.7834;
              496.0000    0.7860;  497.0000    0.7826;  498.0000    0.7878;
              499.0000    0.7901;  500.0000    0.7807;  501.0000    0.7846;
              502.0000    0.7843;  503.0000    0.7769;  504.0000    0.7742;
              505.0000    0.7733;  506.0000    0.7781;  507.0000    0.7703;
              508.0000    0.7696;  509.0000    0.7690;  510.0000    0.7564;
              511.0000    0.7480;  512.0000    0.7391;  513.0000    0.7434;
              514.0000    0.7313;  515.0000    0.7274;  516.0000    0.7193;
              517.0000    0.7076;  518.0000    0.6921;  519.0000    0.6766;
              520.0000    0.6718;  521.0000    0.6542;  522.0000    0.6408;
              523.0000    0.6221;  524.0000    0.6136;  525.0000    0.5913;
              526.0000    0.5770;  527.0000    0.5643;  528.0000    0.5444;
              529.0000    0.5335;  530.0000    0.5138;  531.0000    0.5032;
              532.0000    0.4908;  533.0000    0.4767;  534.0000    0.4640;
              535.0000    0.4604;  536.0000    0.4482;  537.0000    0.4403;
              538.0000    0.4360;  539.0000    0.4290;  540.0000    0.4244;
              541.0000    0.4299;  542.0000    0.4122;  543.0000    0.4287;
              544.0000    0.4302;  545.0000    0.4318;  546.0000    0.4266;
              547.0000    0.4365;  548.0000    0.4426;  549.0000    0.4391;
              550.0000    0.4504;  551.0000    0.4601;  552.0000    0.4671;
              553.0000    0.4753;  554.0000    0.4763;  555.0000    0.4851;
              556.0000    0.5006;  557.0000    0.4998;  558.0000    0.5038;
              559.0000    0.4982;  560.0000    0.5060;  561.0000    0.5097;
              562.0000    0.5115;  563.0000    0.5100;  564.0000    0.5113;
              565.0000    0.5151;  566.0000    0.5059;  567.0000    0.4973;
              568.0000    0.4916;  569.0000    0.4930;  570.0000    0.4813;
              571.0000    0.4715;  572.0000    0.4652;  573.0000    0.4594;
              574.0000    0.4494;  575.0000    0.4441;  576.0000    0.4377;
              577.0000    0.4428;  578.0000    0.4477;  579.0000    0.4600;
              580.0000    0.4713;  581.0000    0.5017;  582.0000    0.5244;
              583.0000    0.5573;  584.0000    0.5872;  585.0000    0.6208;
              586.0000    0.6719;  587.0000    0.7162;  588.0000    0.7447;
              589.0000    0.8076;  590.0000    0.8653;  591.0000    0.9152;
              592.0000    0.9774;  593.0000    1.0215;  594.0000    1.0827;
              595.0000    1.1396;  596.0000    1.1990;  597.0000    1.2383;
              598.0000    1.3008;  599.0000    1.3749;  600.0000    1.4345;
              601.0000    1.4507;  602.0000    1.4921;  603.0000    1.5145;
              604.0000    1.5511;  605.0000    1.6079;  606.0000    1.6005;
              607.0000    1.6308;  608.0000    1.6634;  609.0000    1.7148;
              610.0000    1.7336;  611.0000    1.7498;  612.0000    1.7823;
              613.0000    1.7995;  614.0000    1.8219;  615.0000    1.8488;
              616.0000    1.8656;  617.0000    1.8671;  618.0000    1.8817;
              619.0000    1.9012;  620.0000    1.9292;  621.0000    1.9491; 
              622.0000    1.9696;  623.0000    1.9676;  624.0000    1.9818;
              625.0000    2.0029;  626.0000    2.0102;  627.0000    2.0126;
              628.0000    2.0095;  629.0000    2.0221;  630.0000    2.0414;
              631.0000    2.0690;  632.0000    2.0347;  633.0000    2.0811;
              634.0000    2.0720;  635.0000    2.0835;  636.0000    2.0946;
              637.0000    2.0809;  638.0000    2.1022;  639.0000    2.1079;
              640.0000    2.1160;  641.0000    2.1332;  642.0000    2.1058;
              643.0000    2.1403;  644.0000    2.1447;  645.0000    2.1412;
              646.0000    2.1357;  647.0000    2.1657;  648.0000    2.1699;
              649.0000    2.1933;  650.0000    2.1874;  651.0000    2.1985;
              652.0000    2.1716;  653.0000    2.1855;  654.0000    2.1854;
              655.0000    2.1775;  656.0000    2.2080;  657.0000    2.2116;
              658.0000    2.2192;  659.0000    2.2191;  660.0000    2.2303;
              661.0000    2.2363;  662.0000    2.2349;  663.0000    2.2369;
              664.0000    2.2460;  665.0000    2.2533;  666.0000    2.2734;
              667.0000    2.2596;  668.0000    2.2652;  669.0000    2.2529;
              670.0000    2.2632;  671.0000    2.2764;  672.0000    2.2658;
              673.0000    2.2895;  674.0000    2.3142;  675.0000    2.2744;
              676.0000    2.2977;  677.0000    2.2977;  678.0000    2.3037;
              679.0000    2.3076;  680.0000    2.2961;  681.0000    2.2959;
              682.0000    2.3043;  683.0000    2.3246;  684.0000    2.2999;
              685.0000    2.3268;  686.0000    2.3174;  687.0000    2.3113;
              688.0000    2.3398;  689.0000    2.3388;  690.0000    2.3332;
              691.0000    2.3209;  692.0000    2.3543;  693.0000    2.3142;
              694.0000    2.3198;  695.0000    2.3554;  696.0000    2.3697;
              697.0000    2.3239;  698.0000    2.3485;  699.0000    2.3541;
              700.0000    2.3827];
          
if strcmp(whichCurve,'Kohl')
    E=Ext.Kohl;
    E(:,2:end)=E(:,2:end)/10; % In centimeters
elseif strcmp(whichCurve,'Dunn')
    E=Ext.Dunn;
    E(:,2:end)=E(:,2:end)/10; % In centimeters
elseif strcmp(whichCurve, 'Silico')
    E = Ext.Silico;
    E(:,2) = E(:,2)/10;  %In CM
elseif strcmp(whichCurve, 'Approx')
    [eHbO, eHbR] = ioi_get_extinctions(L1, L2, N);
    bCHbO = 60e-6;
    bCHbR = 40e-6;
    W = linspace(L1, L2, N);
    musp = 12.0*(W/560).^(-2);
    E(:,1) = W;
    E(:,2) = 0.5*sqrt(3*musp./(eHbO*bCHbO + eHbR*bCHbR))/10;
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
%%% Modifications Samuel Belanger (2020/03/17): pas vrai pour Dunn. Voir article (Dunn
%%% 2005). Concentration de base Hbo = 60uM/L et HbR = 40uM/L...
if strcmp(whichCurve,'Kohl')
    x=[0.01 0.02 0.04 0.06 0.08 0.10 0.14 0.18]; %concentration d'h�moglobine
    y=[  .33 .21 .115 .09 .075 0.062 .05 .04 ]; % valeur Da=[OD*mm] � 480 nm
    correction=interp1(x,y,cHBT)/interp1(x,y,.02);
    pathlength=pathlength*correction;
end

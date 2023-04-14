clc; clear; close all;
load('./0.5_0.1_0.98175_65_Underbrink.mat');
N = 50;
z0 = 2;
f = 2000;
SNR = 40;
source = [10,2;20,2];

[DAS_result, a, CSM] = DAS(N, z0, f, coordinates, source, SNR);

figure(1);
contourf(abs(DAS_result));

DAMAS_result = MYDAMAS(DAS_result, a, 100);
figure(2);
contourf(abs(DAMAS_result));

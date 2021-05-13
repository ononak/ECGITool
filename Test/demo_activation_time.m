clc
clear
addpath('../ActivationRecoveryTime');
addpath('../Filter');
addpath('../Regularization');
load('utah.mat') 

[ lat] = LocalActTime(ep);

figure(1);
patch('Vertices',epigeom.pts','Faces',epigeom.fac','FaceVertexCData',lat,'FaceColor','interp');
colorbar
title('Local Activation Time Method')

[scat] = SpCoherentActTime(ep, epigeom);

figure(2);
patch('Vertices',epigeom.pts','Faces',epigeom.fac','FaceVertexCData',scat,'FaceColor','interp');
colorbar
title('Spationaly Coherent Activation Time Method')



[ L ] = SurfaceLaplacian(epigeom);
[stat] = SpatioTemporalAT(ep, L);
figure(3);
patch('Vertices',epigeom.pts','Faces',epigeom.fac','FaceVertexCData',stat,'FaceColor','interp');
colorbar
title('Spatio-Temporal Activation Time Method')
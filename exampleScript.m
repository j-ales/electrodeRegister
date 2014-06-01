%%%%%%%%%% EXAMPLE SCRIPT
%%%% Load data
load('exampleData.mat')

%% Fit to a hi resolution scalp using electrodes and extra digitized scalp points


[rotationMtx, movedPoints, movedHead, maxError, movedElecLo, movedElecHi] ...
    = registerElectrodesToScalp(hiresScalp,elecFiducials,mriFiducials,electrodes,headShapePoints);

figure(1)
clf;
patch('faces',hiresScalp.faces,'vertices',hiresScalp.vertices,'linestyle','none','facecolor',[.8 .7 .6]);
setPlotOptions;
hold on;
scatter3( movedHead(:,1), movedHead(:,2),movedHead(:,3), 80, [.7 .7 0], 'filled' ),

for iPar =1:size(movedElecLo,3),
    
    for iElec = 1:size(movedPoints,1),
        
        line2plot = [ movedElecLo(iElec,:,iPar); ...
            movedElecHi(iElec,:,iPar) ];
        line(line2plot(:,1),line2plot(:,2),line2plot(:,3),'linewidth',2,'color',[.8 .4 0]);
    end

    
end


%% Fit to a hi resolution scalp using electrodes without extra points on scalp
[rotationMtx, movedPoints, movedHead, maxError, movedElecLo, movedElecHi] ...
    = registerElectrodesToScalp(hiresScalp,elecFiducials,mriFiducials,electrodes,[]);

figure(2)
clf;
patch('faces',hiresScalp.faces,'vertices',hiresScalp.vertices,'linestyle','none','facecolor',[.8 .7 .6]);
setPlotOptions;
hold on;

for iPar =1:size(movedElecLo,3),
     
    for iElec = 1:size(movedPoints,1),
        
        line2plot = [ movedElecLo(iElec,:,iPar); ...
            movedElecHi(iElec,:,iPar) ];           
        line(line2plot(:,1),line2plot(:,2),line2plot(:,3),'linewidth',2,'color',[.8 .4 0]);
    end

    
end

 
%% Fit to a low resolution scalp using electrodes and extra points on scalp
[rotationMtx, movedPoints, movedHead, maxError, movedElecLo, movedElecHi] ...
    = registerElectrodesToScalp(loresScalp,elecFiducials,mriFiducials,electrodes,headShapePoints);

figure(3)
clf;

patch('faces',loresScalp.faces,'vertices',loresScalp.vertices,'linestyle','none','facecolor',[.8 .7 .6]);
setPlotOptions;
hold on;
scatter3( movedHead(:,1), movedHead(:,2),movedHead(:,3), 80, [.7 .7 0], 'filled' )

for iPar =1:size(movedElecLo,3),
     
    for iElec = 1:size(movedPoints,1),
        
        line2plot = [ movedElecLo(iElec,:,iPar); ...
            movedElecHi(iElec,:,iPar) ];           
        line(line2plot(:,1),line2plot(:,2),line2plot(:,3),'linewidth',2,'color',[.8 .4 0]);
    end

    
end

%% Fit to a low resolution scalp using just electrodes and fiducials
[rotationMtx, movedPoints, movedHead, maxError, movedElecLo, movedElecHi] ...
    = registerElectrodesToScalp(loresScalp,elecFiducials,mriFiducials,electrodes,[]);

figure(4)
clf;
patch('faces',loresScalp.faces,'vertices',loresScalp.vertices,'linestyle','none','facecolor',[.8 .7 .6]);
setPlotOptions;

for iPar =1:size(movedElecLo,3),
     
    for iElec = 1:size(movedPoints,1),
        
        line2plot = [ movedElecLo(iElec,:,iPar); ...
            movedElecHi(iElec,:,iPar) ];           
        line(line2plot(:,1),line2plot(:,2),line2plot(:,3),'linewidth',2,'color',[.8 .4 0]);
    end

    
end


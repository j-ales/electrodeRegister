%%%%%%%%%% EXAMPLE SCRIPT


%%%% Load data
load('exampleData.mat')

%% Fit to a hi resolution scalp using electrodes and extra digitized scalp points

%hiresScalp: Scalp surface found using FreeSurfer
%elecFiducials: Digitized fiducial points: Left ear (tragus), Right Ear (tragus), Nasion
%mriFiducials: Corresponding fiducial points as they are in the MRI frame
%electrodes: digitized locations for electrodes
%headShapePoints: digitized points from around the scalp.

[rotationMtx, movedPoints, movedHead, maxError, movedElecLo, movedElecHi] ...
    = registerElectrodesToScalp(hiresScalp,elecFiducials,mriFiducials,electrodes,headShapePoints);

figure(1)
clf;
patch('faces',hiresScalp.faces,'vertices',hiresScalp.vertices,'linestyle','none','facecolor',[.8 .7 .6]);
setPlotOptions;
hold on;
hs=scatter3( movedHead(:,1), movedHead(:,2),movedHead(:,3), 80, [.7 .7 0], 'filled' );

for iPar =1:size(movedElecLo,3),
    
    for iElec = 1:size(movedPoints,1),
        
        line2plot = [ movedElecLo(iElec,:,iPar); ...
            movedElecHi(iElec,:,iPar) ];
        hl=line(line2plot(:,1),line2plot(:,2),line2plot(:,3),'linewidth',2,'color',[.8 .4 0]);
    end

    
end

legend([hs;hl],'Headshape points','Electrode Confidence Interval');


%% Fit to a hi resolution scalp using electrodes without extra points on scalp
%hiresScalp: Scalp surface found using FreeSurfer
%elecFiducials: Digitized fiducial points: Left ear (tragus), Right Ear (tragus), Nasion
%mriFiducials: Corresponding fiducial points as they are in the MRI frame
%electrodes: digitized locations for electrodes
%headShapePoints: NOT USED.

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
       hl= line(line2plot(:,1),line2plot(:,2),line2plot(:,3),'linewidth',2,'color',[.8 .4 0]);
    end

    
end
legend([hl],'Electrode Confidence Interval');

 
%% Fit to a low resolution scalp using electrodes and extra points on scalp
%loresScalp: Scalp surface found using FSL BET
%elecFiducials: Digitized fiducial points: Left ear (tragus), Right Ear (tragus), Nasion
%mriFiducials: Corresponding fiducial points as they are in the MRI frame
%electrodes: digitized locations for electrodes
%headShapePoints: Points collectred from the scalp surface.

[rotationMtx, movedPoints, movedHead, maxError, movedElecLo, movedElecHi] ...
    = registerElectrodesToScalp(loresScalp,elecFiducials,mriFiducials,electrodes,headShapePoints);

figure(3)
clf;

patch('faces',loresScalp.faces,'vertices',loresScalp.vertices,'linestyle','none','facecolor',[.8 .7 .6]);
setPlotOptions;
hold on;
hs=scatter3( movedHead(:,1), movedHead(:,2),movedHead(:,3), 80, [.7 .7 0], 'filled' );

for iPar =1:size(movedElecLo,3),
     
    for iElec = 1:size(movedPoints,1),
        
        line2plot = [ movedElecLo(iElec,:,iPar); ...
            movedElecHi(iElec,:,iPar) ];           
        hl=line(line2plot(:,1),line2plot(:,2),line2plot(:,3),'linewidth',2,'color',[.8 .4 0]);
    end

    
end
legend([hs;hl],'Headshape points','Electrode Confidence Interval');

%% Fit to a low resolution scalp using just electrodes and fiducials
[rotationMtx, movedPoints, movedHead, maxError, movedElecLo, movedElecHi] ...
    = registerElectrodesToScalp(loresScalp,elecFiducials,mriFiducials,electrodes,[]);
%loresScalp: Scalp surface found using FSL BET
%elecFiducials: Digitized fiducial points: Left ear (tragus), Right Ear (tragus), Nasion
%mriFiducials: Corresponding fiducial points as they are in the MRI frame
%electrodes: digitized locations for electrodes
%headShapePoints: NOT USED.

figure(4)
clf;
patch('faces',loresScalp.faces,'vertices',loresScalp.vertices,'linestyle','none','facecolor',[.8 .7 .6]);
setPlotOptions;

for iPar =1:size(movedElecLo,3),
     
    for iElec = 1:size(movedPoints,1),
        
        line2plot = [ movedElecLo(iElec,:,iPar); ...
            movedElecHi(iElec,:,iPar) ];           
        hl=line(line2plot(:,1),line2plot(:,2),line2plot(:,3),'linewidth',2,'color',[.8 .4 0]);
    end

    
end

legend([hl],'Electrode Confidence Interval');

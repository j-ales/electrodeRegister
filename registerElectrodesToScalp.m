function [transformMtx, movedElec, movedHeadshape, maxError, movedElecLo, movedElecHi] = registerElectrodesToScalp(scalp,elecFiducials,mriFiducials,electrodes,headshape)

%fitPointsToScalp  - Find a rigid transformation the best fits 2 sets of points
%function [transformMtx, movedElec, movedHeadshape, maxError, movedElecLo, movedElecHi] = ... 
%    registerElectrodesToScalp(scalp,fiducials,electrodes,[headshape])
% 
% This script does a 6 degree of freedom fit of electrodes to a scalp
% surface. 
%
% All units expected to be in mm
%
% INPUTS:
% scalp = Matlab structure with 2 fields: faces & vertices
% fiducials = 3x3 set of fiducial
% electrodes = nx3 electrode coordinates
% optional:
% headshape  = nx3 set of points taken from the scalp 
%
% OUTPUTS:
% transformMtx  (4x4)  - The transformation matrix that transforms the electrode
%                        coordinates to the scalp
% movedElec     (nx3)  - The transformed electrodes
% movedHeadshape (nx3) - The transformed headshape points
% maxError    (scalar) - The maximum electrode distance error given 95% confidence
%                        on the final fit.
% movedElecLo (nx3x6)  -
% movedElecHi (nx3x6)  - These two matrices contain the electrodes
%                        transformed by parameters corresponding to the the range 
%                        for the 95% confidence there are 6 parameters: 3 translations 
%                        3 rotations. Each (nx3x1) points corresponds to
%                        changing one paramter while holding the other parameters 
%                        fixed at their optimum values
%
%Depends on nearpoints from vistasoft: https://github.com/vistalab/vistasoft
%Requires the statistics and optimization toolboxes from matlab

% get all installed toolbox names
v = ver;
% collect the names in a cell array
[installedToolboxes{1:length(v)}] = deal(v.Name);

% check 
tf = all(ismember('Optimization Toolbox',installedToolboxes));

if tf ==false
   error('This function requires the Optimization Toolbox, which is not installed') 
end

tf = all(ismember('Statistics Toolbox',installedToolboxes));

if tf ==false
   error('This function requires the Statistics Toolbox, which is not installed') 
end


%Set some options for the search
optns = optimset(@lsqnonlin);
optns = optimset(optns,'display', 'off', 'maxfunevals', 1000, 'MaxIter', 100,...
    'FinDiffType','central','LargeScale','on','TolX',1e-7,'TolFun',1e-7);

sf = 1/100; % Scale factor for making mm units similar to radians, this improves the nonlinear search
            % probably not needed in modern versions of matlab, but old
            % habits. JMA
stationaryPoints = scalp.vertices*sf;
elecFiducials = elecFiducials*sf;
mriFiducials = mriFiducials*sf;
electrodes = electrodes*sf;

%Checking to see if extra scalp points should be used in fitting
if nargin<4 || isempty(headshape),
    headshape = [];
    movedHeadshape = [];
end
headshape = headshape*sf;

%Find intitial registration using the known fiducial points.
[tOrig rOrig] = alignFiducials(elecFiducials,mriFiducials);
transOrig = [ rOrig, tOrig'; 0 0 0 1];
transFid = transOrig*[elecFiducials, [1; 1; 1;]]';
transFid = [transFid(1:3,:)]';
t = tOrig;
r = rOrig;
%Transformation matrix
trans = [ r, t'; 0 0 0 1];
    
%Fiducials, electrodes and Headshape points translated to
%initial conditions found above using just the mri and elec fiducials.
electrodes = applyTransform(trans,electrodes);
elecFiducials = applyTransform(trans,elecFiducials);
initialTransform = trans;

            
          
% If there are scalp points apply the initial transform
if ~isempty(headshape)
    headshape = applyTransform(trans,headshape);
    
    %Do an intial QA check for scalp points that aren't close.
    %Sometimes people pull the stylus off the scalp before discontinuing
    %digitization
    [K,D] = nearpoints(headshape', stationaryPoints');
    
    %throw out points farther than 3 cm from scalp
    npointsRem = sum(full(sqrt(D)>(30*sf)));
    disp([num2str(npointsRem) ' removed because they are further than 3 cm from inital scalp']); 
    headshape = headshape(sqrt(D)<(sf*30),:);
    

end




%This code is a quick way to calculate surface normals.
%Surface normals are used to be able to calculate if electrodes are inside
%the head.
%
scalpCenter = mean(scalp.vertices);
fig = figure;

N = get(patch('vertices',scalp.vertices,'faces',scalp.faces(:,[3 2 1])),'vertexnormals')';
N = N ./ repmat(sqrt(sum(N.^2)),3,1);
signedDist = dot(  bsxfun(@minus,scalp.vertices,scalpCenter)', N);

%We must also check that the structure given to us uses the appropriate
%convention for the winding of faces.  We want positive surface normals to
%point outwards. So we assume the midpoint of the scalp surface will be
%more or less inside and check if most normals point away from the center.
%
if mean(sign(signedDist))<0 %Checking that more than 50% normals point in the "correct" directions
    disp('Surface normals found with opposite assumed winding, flipping normals')
    N = get(patch('vertices',scalp.vertices,'faces',scalp.faces(:,[1 2 3])),'vertexnormals')';
    N = N ./ repmat(sqrt(sum(N.^2)),3,1);
end
    



close(fig)

initial = [0 0 0 0 0 0];

%[params fval] = fminsearch(@translate,[0 0 0 0],optns,v1fMRI,VEP);
%[params fval EXITFLAG, OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(@rotcostfunc2,initial,[],[],[],[],lowlim,uplim,[],optns, ...
%						 stationaryPoints,electrodes,headshape,N);


% %Constrained, uses large scale algorithm
%Some mild constraints on the search space.
% lowlim = [ -100*sf -100*sf -100*sf -pi/2 -pi/2 -pi/2];
% uplim  = [  100*sf  100*sf  100*sf  pi/2  pi/2  pi/2];
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(@rotcostfunclsq,initial,lowlim,uplim,optns, ...
% 						 stationaryPoints,electrodes,headshape,N);

%Unconstrained uses line-search
disp('Registering Electrodes.. . .   .       .');
[X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(@rotcostfunclsq,initial,[],[],optns, ...
						 stationaryPoints,electrodes,headshape,N);
disp('Registration Done');

params = X;

sigmaHat = RESNORM*length(RESIDUAL); %Sigma^2 estimator of measurement noise from data

%Not used currently, but useful to know.  The following code will calculate
%the covariance matrix for the parameter fits.
covMat = sigmaHat*inv(JACOBIAN'*JACOBIAN); %Estimated Covariance Matrix;

%Easy way is to use the matlab function nlparci to interpret to outputs of
%lsqnonlin
CI = nlparci(X,RESIDUAL,'jacobian',JACOBIAN);

%Must undo the scaling used for the fitting procedure.
if ~isempty(headshape)
    movedHeadshape = rigidRotate(params,headshape)./sf;
end

movedElec = rigidRotate(params,electrodes)./sf;
params(1:3) = params(1:3)./sf;

transformMtx = makeRotMtx(params);

%Total transform is the initial transform plus the final transform.
transformMtx = transformMtx*initialTransform;


%Know let's find how much we can shift electrodes and still have the cost 
%be essentially the same. We're going to use the 95% confidence intervals
%calculated above. We will change each parameter individual to the extreme
%of the 95% ranges and see how much that moves the electrodes.
%Then we'll take the largest movement over all electrodes and parameters
%and call that our margin of error.
%
%Now, this ignores possible covariance between parameters, but it's a
%fairly good metric of how well the cost function pins down the
%registration.

for iPar =1:length(CI),
    
    thisParLo = CI(iPar,1);
    thisParHi = CI(iPar,2);
    %    parSamp(iPar,:) = linspace(CI(iPar,1),CI(iPar,2),nP);

    
    thisXLo = X;
    thisXHi = X;

    thisXLo(iPar) = thisParLo;
    thisXHi(iPar) = thisParHi;

    movedElecLo(:,:,iPar) = rigidRotate(thisXLo,electrodes)./sf;
    movedElecHi(:,:,iPar) = rigidRotate(thisXHi,electrodes)./sf;

%    dist(iPar) = max(sqrt(sum((movedElecLo(:,:,iPar)-movedElecHi(:,:,iPar)).^2,2)));
    
    dist1(iPar) = max(sqrt(sum((movedElecLo(:,:,iPar)-movedElec(:,:)).^2,2)));
    dist2(iPar) = max(sqrt(sum((movedElec(:,:)-movedElecHi(:,:,iPar)).^2,2)));
    
    dist(iPar) = max(dist1(iPar),dist2(iPar));
    
end

%axLim = axis;
maxError = max(dist);

%title(['Electrodes fit to within: ' num2str(max(dist/2)) ' mm'])
%text(axLim(2),axLim(3),axLim(6),['Electrodes fit to within: ' num2str(max(dist/2),2) ' mm']);
disp(['Electrodes fit to within: ' num2str(maxError,3) ' mm'])


end 


function dist = rotcostfunc(params,stationaryPoints,movablePoints)

	movedPoints = rigidRotate(params,movablePoints);
    %[K,D] = nearpoints(src, dest) 
    
    [K,D] = nearpoints(movedPoints', stationaryPoints'); 
  % [K,D] = nearpoints(stationaryPoints',movedPoints'); 
  
 
  %[K,D] = dsearchn(stationaryPoints,movedPoints);
    
    
    %dist = sum(D(sqrt(D)<20));
    dist = sum(D);
    
	%dist = sum(sum((stationaryPoints - movedPoints).^2));
end

    
function dist = rotcostfunc2(params,stationaryPoints,electrodes,headshape,N)

movedElec = rigidRotate(params,electrodes);


[kE,dE] = nearpoints(movedElec', stationaryPoints');

signedDist = dot(  (movedElec - stationaryPoints(kE,:))', N(:,kE));
%    SSE = sum(signedDist-mean(signedDist)).^2;


dE = sqrt(dE);
dE(dE<30) = 30;
dE = (dE-30).^2; 


if ~isempty(headshape)
    movedHead = rigidRotate(params,headshape);
    [kH,dH] = nearpoints(movedHead', stationaryPoints');
   % [mean(dH) var(signedDist)]
    dist = mean(dH)+var(signedDist);
else
    
    dist = var(signedDist) +sum(dE);
end



	%dist = sum(sum((stationaryPoints - movedPoints).^2));

	
%{
	hold off;
	clf;
    scatter3(stationaryPoints(:,1),stationaryPoints(:,2),stationaryPoints(:,3),'r')
	hold on
	scatter3(movedPoints(:,1),movedPoints(:,2),movedPoints(:,3),'k')
	scatter3(movablePoints(:,1),movablePoints(:,2),movablePoints(:,3),'b')
	drawnow;
pause

%}
end 
    
    
        
function [dist meanDist] = rotcostfunclsq(params,stationaryPoints,electrodes,headshape,N)

movedElec = rigidRotate(params,electrodes);


[kE,dE] = nearpoints(movedElec', stationaryPoints');

signedDist = dot(  (movedElec - stationaryPoints(kE,:))', N(:,kE));
%    SSE = sum(signedDist-mean(signedDist)).^2;

meanDist = mean(signedDist);

% figure(1);
% hold on;
% plot(signedDist)
% 
% figure(2);
% hold on;
% plot(dE)
% 
% figure(3)
% hold on;
% scatter3(movedElec(:,1),movedElec(:,2),movedElec(:,3),'k')
%[K,D] = dsearchn(stationaryPoints,movedPoints);

if ~isempty(headshape)
    movedHead = rigidRotate(params,headshape);
    [kH,dH] = nearpoints(movedHead', stationaryPoints');
   % [mean(dH) var(signedDist)]
    dist = [dH [signedDist-mean(signedDist)]];
else
    dE = sqrt(dE);
    dE(dE<30) = 30;
    dE = (dE-30);
    dist = [dE + [signedDist - meanDist]];
end

end


function rotMtx = makeRotMtx(params)

 xShift = params(1);
 yShift = params(2);
 zShift = params(3);
 a   = params(4);
 b   = params(5);
 g   = params(6);

 Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a) ];
 Ry = [ cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b) ];
 Rz = [ cos(g) -sin(g) 0;  sin(g) cos(g) 0;  0 0 1 ];


 newR = Rx*Ry*Rz;


rotMtx = [newR, [xShift yShift zShift]'; 0 0 0 1];

end


function transPoints = applyTransform(trans,points)

points = [points ones(length(points),1)];
transPoints = trans*points';
transPoints = transPoints(1:3,:)';

end




function movedPoints = rigidRotate(params,movablePoints)
%rigidRotate - performs a rigid body transformation
%movedPoints = rigidRotate(params,movablePoints)
%
%
%This function that performs a rigid body transformation on movablePoints
%that is specified as follows:
%
%   movablePoints - Nx3 matrix
%
%	xShift = params(1);
%	yShift = params(2);
%	zShift = params(3);
%	Rx   = params(4); rotation around x axis
%	Ry   = params(5); rotation around y
%	Rz   = params(6); rotation around z

xShift = params(1);
yShift = params(2);
zShift = params(3);
ang1   = params(4);
ang2   = params(5);
ang3   = params(6);


movedPoints =  Rx(movablePoints,ang1);
movedPoints =  Ry(movedPoints,ang2);
movedPoints =  Rz(movedPoints,ang3);

movedPoints = 	[movedPoints(:,1)+xShift, ...
    movedPoints(:,2)+yShift, ...
    movedPoints(:,3)+zShift];

end

% R? functions:
% Rx - Rotate 3D Cartesian coordinates around the X axis
% Ry - Rotate 3D Cartesian coordinates around the Y axis
% Rz - Rotate 3D Cartesian coordinates around the Z axis
%
% Useage: [XYZ] = rx(XYZ,alpha,units)
%
% XYZ is a [3,N] or [N,3] matrix of 3D Cartesian coordinates
%
% 'alpha' - angle of rotation about the axis
% 'units' - angle is either 'degrees' or 'radians'
%           the default is alpha in radians
%
% If input XYZ = eye(3), the XYZ returned is
% the rotation matrix.
% $Revision: 1.10 $ $Date: 2004/04/16 18:49:10 $
% Licence:  GNU GPL, no express or implied warranties
% History:  04/2002, Darren.Weber_at_radiology.ucsf.edu
%                    Developed after example 3.1 of
%                    Mathews & Fink (1999), Numerical
%                    Methods Using Matlab. Prentice Hall: NY.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [XYZ] = Rx(XYZ,a,units)
        
        
        
        if ~exist('units','var'), units = 'radians'; end
        
        % convert degrees to radians
        if isequal(units,'degrees'),
            a = a*pi/180;
        end
        
        Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a) ];
        
        if isequal(size(XYZ,1),3),
            XYZ = Rx * XYZ;
        else
            XYZ = XYZ';
            if isequal(size(XYZ,1),3),
                XYZ = [Rx * XYZ]';
            else
                error('Rx: Input XYZ must be [N,3] or [3,N] matrix.\n');
            end
        end
        
    end


    function [XYZ] = Ry(XYZ,b,units)
        
        if ~exist('units','var'), units = 'radians'; end
        
        % convert degrees to radians
        if isequal(units,'degrees'),
            b = b*pi/180;
        end
        
        Ry = [ cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b) ];
        
        if isequal(size(XYZ,1),3),
            XYZ = Ry * XYZ;
        else
            XYZ = XYZ';
            if isequal(size(XYZ,1),3),
                XYZ = [Ry * XYZ]';
            else
                error('Ry: Input XYZ must be [N,3] or [3,N] matrix.\n');
            end
        end
        
    end

    
    function [XYZ] = Rz(XYZ,g,units)
        
        
        if ~exist('units','var'), units = 'radians'; end
        
        % convert degrees to radians
        if isequal(units,'degrees'),
            g = g*pi/180;
        end
        
        Rz = [ cos(g) -sin(g) 0;  sin(g) cos(g) 0;  0 0 1 ];
        
        if isequal(size(XYZ,1),3),
            XYZ = Rz * XYZ;
        else
            XYZ = XYZ';
            if isequal(size(XYZ,1),3),
                XYZ = [Rz * XYZ]';
            else
                error('Rz: Input XYZ must be [N,3] or [3,N] matrix.\n');
            end
        end
        
    end
    
    

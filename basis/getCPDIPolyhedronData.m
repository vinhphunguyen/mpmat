function [res,V, facesData] = getCPDIPolyhedronData(faces,xp,particle,res)
% Compute CPDI-ngon data at particle "pid".
%
% Inputs:
%
% face:     corner indices of the face
% center:   centroid of the polyhedron
% particle: particle mesh
%
% Output:
%
% data.wf:    function weights
% data.wf:    gradient weights
%
% Call functions:
%   getCPDIFaceData
%
% VP Nguyen
% 20 April 2016
% Monash University

faceCount = length(faces);

V = 0;

facesData.coord = zeros(faceCount,3); % centers of faces
facesData.wf    = zeros(faceCount,1); % weights of face centers
facesData.wg    = zeros(faceCount,3); % gradient weights of face centers

% loop over the faces
for f=1:faceCount
    face            = faces{f};
    [res,vol,fRes]  = getCPDIFaceData(face,xp,particle,res);
    V               = V + vol;
    facesData.coord(f,:) = fRes.coord;
    facesData.wf(f)      = fRes.wf;
    facesData.wg(f,:)    = fRes.wg;
end




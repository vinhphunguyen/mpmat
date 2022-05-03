function [res,V,fRes] = getCPDIFaceData(face,center,particle,res)
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
%   getCPDITet4Data()
%
% VP Nguyen
% 20 April 2016
% Monash University

% corners of the particle domain

corners   = particle.node(face,:);
nodeCount = length(face);
faceCenter= mean(corners);

V = 0;

fRes.wf    = 0;
fRes.wg    = [0 0 0];
fRes.coord = faceCenter;

% loop over the edges of the face
for i = 1:nodeCount
  j    = rem(i,nodeCount) + 1;

  % form a sub-tetrahedron by connecting edge ij with the face center and
  % the centroid
  cornersI = [corners(i,1) corners(i,2) corners(i,3);...
              corners(j,1) corners(j,2) corners(j,3);...
              faceCenter;              
              center];
          
  data = getCPDITet4Data(cornersI);
  
  first = 1; secon = 2;
  
  if ( data.vol < 0 )
      cornersI = [corners(j,1) corners(j,2) corners(j,3);...
                  corners(i,1) corners(i,2) corners(i,3);...
                  faceCenter;
                  center];
      
      data = getCPDITet4Data(cornersI);
      first = 2; secon = 1;
  end
  if ( data.vol < 0 )
      error('Impossible');
  end
  V    = V + data.vol;

  
  % global indices of corners i and j
  ig          = face(i); 
  jg          = face(j); 
  
  res.wf(ig)   = res.wf(ig)   + data.wf(first);
  res.wf(jg)   = res.wf(jg)   + data.wf(secon);
  res.wf(end)  = res.wf(end)  + data.wf(4); % for the particle centroid
  fRes.wf      = fRes.wf      + data.wf(3); % for the face center
  
  res.wg(ig,:)   = res.wg(ig,:)   + data.wg(first,:);
  res.wg(jg,:)   = res.wg(jg,:)   + data.wg(secon,:);
  res.wg(end,:)  = res.wg(end,:)  + data.wg(4,:); % for the particle centroid
  fRes.wg        = fRes.wg        + data.wg(3,:); % for the face center
end




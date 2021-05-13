function [ D ] = SurfaceGradient(geom)
% SURFACE_GRADIENT  Calculates the surface Gradient operator for 3D mesh
% surface
%
% Usage:   D = SurfaceGradient(geom)
% 
% Input:
%   geom    Geometry of the given surface
%       geom.pts is points (3 x nLeads), 
%       geom.fac are the triangles (size: (3 x numOfTri); 
%               node indexing should start with 1, not 0)
%
% Output:
%   D       Surface Gradient matrix (operator acting on the x vector)
%
% Reference: A Comparison of Gradient Estimation Methods for Volume Rendering on Unstructured Meshes
%                Carlos D. Correa, Robert Hero, and Kwan-Liu Ma,
%               DOI: 10.1109/TVCG.2009.105
%
% Author: Onder Nazim Onak ononak@gmail.com
%


%initialize
neighbours = NeighbourList(geom.fac);
D = zeros(size(neighbours,1));

    for leadNo = 1:size(neighbours,1)       
        currentLeadCoordinate = geom.pts(:,leadNo);
        den = zeros(1,neighbours(leadNo,2));
       
        % for each neighbour node
        for nbrIdx = 1:neighbours(leadNo,2)
            nbrVtx = neighbours(leadNo, nbrIdx +2);
            nbrCoordinate = geom.pts(:,nbrVtx);
            den(nbrIdx) = norm(nbrCoordinate - currentLeadCoordinate);%^2; % weight
        end
        
        den_normalized = den/sum(den);

        for nbrIdx = 1:neighbours(leadNo,2)
            nbrVtx = neighbours(leadNo, nbrIdx +2);
            D(leadNo,nbrVtx) = den_normalized(nbrIdx);
            D(leadNo,leadNo) = D(leadNo,leadNo) - den_normalized(nbrIdx);
            %D(leadNo,leadNo) = D(leadNo,leadNo) - den(nbrIdx);
        end
    end

end


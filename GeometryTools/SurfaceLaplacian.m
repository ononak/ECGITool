function [ L ] = SurfaceLaplacian(geom)
% SURFACE_LAPLACIAN  Calculates the surface Laplacian operator for 3D mesh
%
% Usage:   L = surface_laplacian(geom)
% 
% Input:
%   geom    Geometry of the given surface
%       geom.pts is points (3 x nLeads), 
%       geom.fac are the triangles (size: (3 x numOfTri); 
%               node indexing should start with 1, not 0)
%
% Output:
%   L       Surface Gradient matrix (operator acting on the x vector)
%
% Reference: Gabriel Peyré. Numerical Mesh Processing. 2008. hal-00365931
%
% Author: Onder Nazim Onak ononak@gmail.com
%
%initialize
neighbours = NeighbourList(geom.fac);
D = zeros(size(neighbours,1));
W = zeros(size(neighbours,1));

    for leadNo = 1:size(neighbours,1)       
        currentLeadCoordinate = geom.pts(:,leadNo);
        den = zeros(1,neighbours(leadNo,2));
        
        % for each neighbour node
        for nbrIdx = 1:neighbours(leadNo,2)
            nbrVtx = neighbours(leadNo, nbrIdx +2);
            nbrCoordinate = geom.pts(:,nbrVtx);
            den(nbrIdx) = 1/norm(nbrCoordinate - currentLeadCoordinate) ^2; % weight
        end
        
         D(leadNo,leadNo) = sum(den);
         
        for nbrIdx = 1:neighbours(leadNo,2)
            nbrVtx = neighbours(leadNo, nbrIdx +2);
            W(leadNo,nbrVtx) = den(nbrIdx);
        end
    end

    L = D - W;
end


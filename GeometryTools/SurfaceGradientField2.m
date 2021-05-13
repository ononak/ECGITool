function [GF, GVal] = SurfaceGradientField2(geom, valueVector)
% SURFACEGRADIENTFIELD2  Calculates the surface Gradient field for 3D mesh
%
% Usage:   GF = SurfaceGradientField(geom)
% 
% Input:
%   geom    Geometry of the surface
%       geom.pts is points (3 x nLeads), 
%       geom.fac are the triangles (size: (3 x numOfTri); 
%               node indexing should start with 1, not 0)
%   activationTimeVec  Contains the activation time of the each lead
%   located on the vtx
%
% Output:
%   GF       Surface Gradient field matrix
%
% Reference: A Comparison of Gradient Estimation Methods for Volume Rendering on Unstructured Meshes
%                Carlos D. Correa, Robert Hero, and Kwan-Liu Ma,
%               DOI: 10.1109/TVCG.2009.105
%
% Author: Onder Nazim Onak ononak@gmail.com
%
%initialize
neighbours = NeighbourList(geom.fac);
GF = zeros(size(neighbours,1),3);
GVal = zeros(size(neighbours,1),1);

    for leadNo = 1:size(neighbours,1)     
        numberofNbrVtx = neighbours(leadNo,2);
        currentLeadCoordinate = geom.pts(:,leadNo);
        w = zeros(1,numberofNbrVtx);
        VecDiff = zeros(numberofNbrVtx,3);
        ActDiff = zeros(numberofNbrVtx,1);
        
        % for each neighbour node
        for nbrIdx = 1:neighbours(leadNo,2)
            nbrVtx = neighbours(leadNo, nbrIdx +2);
            nbrCoordinate = geom.pts(:,nbrVtx);
            w(nbrIdx) = 1/norm(nbrCoordinate - currentLeadCoordinate)^2; % weight
            VecDiff(nbrIdx,:) = currentLeadCoordinate' - nbrCoordinate';
            ActDiff(nbrIdx,:) = valueVector(leadNo) - valueVector(nbrIdx);
        end
        
        W = diag(w);
        GF(leadNo,:) = ((VecDiff'*W^2*VecDiff)\(VecDiff'*W^2*ActDiff))';
        
        GVal(leadNo,1) = norm(GF(leadNo,:));
        GF(leadNo,:) = GF(leadNo,:)/norm(GF(leadNo,:));
    end

end


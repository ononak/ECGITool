function [GradientVector] = SurfaceGradientField(geom, funcValue)
% SURFACEGRADIENTFIELD  Generates the gradient field vector for the function defined on
% the 3D mesh surface
%
% 
% Input:
%   geom    Geometry of the given surface
%       geom.pts is points (3 x nLeads), 
%       geom.fac are the triangles (size: (3 x numOfTri); 
%               node indexing should start with 1, not 0)
%  funcValue  surface function values calculated at points geom.pts (example: epicardial potential values)
%             
%
% Output:
%  GradientVector   Surface Gradient matrix applied to funcValue
%
% Reference: A Comparison of Gradient Estimation Methods for Volume Rendering on Unstructured Meshes
%                Carlos D. Correa, Robert Hero, and Kwan-Liu Ma,
%               DOI: 10.1109/TVCG.2009.105
%
% Author: Onder Nazim Onak ononak@gmail.com
%


%initialize
number_of_vtx = size(geom.pts,2);
neighbours = NeighbourList(geom.fac);

GradientVector = zeros(3,number_of_vtx);

    for vertex = 1:number_of_vtx
        
        number_of_nb = neighbours(vertex,2);
        D = zeros(number_of_nb, 3);
        fD = zeros(number_of_nb,1);
        currentLeadCoordinate = geom.pts(:,vertex);
       
        % for each neighbour node
        for nbrIdx = 1:number_of_nb 
            nbrVtx = neighbours(vertex, nbrIdx +2);
            nbrCoordinate = geom.pts(:,nbrVtx);
            D(nbrIdx,:) = (nbrCoordinate - currentLeadCoordinate)';
            fD(nbrIdx,:) = funcValue(nbrVtx) - funcValue(vertex);
        end
        vertex
        GradientVector(:,vertex) = D\fD;
    end

end
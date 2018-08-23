%
% Contents file for package: htm
% Created: 29-Dec-2015
%---------
% gc_mid_section.m :  Given two points on a sphere, find the central point lying on the the shortest great circle section connecting the two points.
% htm_build.m :  Build Hierarchical Triangular Mesh (HTM) structure. This structure can be use for fast searches of data in catalogs on a sphere.
% htm_build_son.m :  An auxilary function for htm_build.m for building Hierarchical Triangular Mesh (HTM) structure.
% htm_search_point.m :  Search for a single point-like coordinate in an HTM tree.
% in_halfspace.m :  Given a unit vector R and half space (N,C) test if the point is contained inside the half space (N dot R > C). A halfspace is a plane that splits the sphere in two. It is defined by a direction vector (N) and a signed scalar (C), measured along the normal vector from the origin of the sphere.
% in_polysphere.m :  Check if a list of positions are found within a convex polygon on the celestial sphere in which its sides are great circles. The polygon should be defined according to the right-hand rule.
% search_htm_coocat.m :  
% tree_collect_leafs.m :  Given a tree and a pointer to a node, collect all the leafs found at the last level below this node.

function [Ptr]=tree_collect_leafs(Tree,StartPtr,Collect)
% Collect leafs in a tree
% Package: celestial.htm
% Description: Given a tree and a pointer to a node, collect all the leafs
%              found at the last level below this node.
% Input  : - A tree structure array (e.g., htm_build.m).
%          - Pointer in the tree below to collect all leafs.
%          - Vector of collected pointers used for recusion.
%            For internal use. Default is empty matrix.
% Output : - Vector of pointers to all tree elements at the last level
%            below the input node (i.e., all leafs that belong to a node).
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Aug 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [HTM,LevList]=celestial.htm.htm_build(4);
%          [Ptr]=tree_collect_leafs(HTM,100)
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin==2),
    Collect = [];
end

Nleaf = numel(Tree(StartPtr).son);
if (Nleaf==0),
   % final leaf found
   Collect = [Collect, StartPtr];
else
   for Ileaf=1:1:Nleaf,
       Collect = celestial.htm.tree_collect_leafs(Tree,Tree(StartPtr).son(Ileaf),Collect);
   end
end

Ptr = Collect;

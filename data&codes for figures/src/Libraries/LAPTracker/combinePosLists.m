function [ posListAll ] = combinePosLists( varargin )
% [ posListAll ] = combinePosLists( posList1, ..., posListN )
%
% Concatenates posLists in the order given,
% updating frame numbers accordlying.
% 
% PARAM
% posList1, ..., posListN posLists to concatonate
%
%RETURNS
% posListAll concatonated pos list with frames updated
%
% Author: $Author: cotter $
% Revision: $Revision: 1.1 $

%$Log: combinePosLists.m,v $
%Revision 1.1  2014/02/11 17:19:11  cotter
%_
%
    names = fieldnames(varargin{1});
    for i = names'
        posListAll.(i{1}) = [];
    end
    
    maxFrame = 0;
    for i = 1:length(varargin)
        posList = varargin{i};
        
        for n = names'
            if(strcmp(n,'frame'))
                posListAll.frame = [posListAll.frame; posList.frame + maxFrame];
            else
                posListAll.(n{1}) = [posListAll.(n{1}); posList.(n{1})];
            end
        end

        maxFrame = max(posList.frame);
    end
end
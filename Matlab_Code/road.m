classdef road < handle
    properties
        starting_point
        ending_point
        length
        width
    end
    methods
        function obj = road(starting_point,ending_point,length,width)
            obj.starting_point = starting_point;
            obj.ending_point = ending_point;
            obj.length = length;
            obj.width = width;     
        end     
    end
end

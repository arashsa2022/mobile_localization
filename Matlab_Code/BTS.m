classdef BTS < handle
    properties
        id
        coordinate
    end
    methods
        function obj = BTS(id, coordinate)
            global BTS_xy;
            obj.id = id;
            obj.coordinate = coordinate;
            BTS_xy(obj.id,:) = coordinate;
        end     
    end
end

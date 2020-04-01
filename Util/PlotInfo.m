classdef PlotInfo
    % Holds information about current function for ploting
    
    properties
        Limits, Minimum, Contours
    end
    
    methods
        function obj = PlotInfo(Limits, Minimum, Contours)
            %PLOTINFO Construct an instance of this class
            obj.Limits = Limits;
            obj.Minimum = Minimum;
            obj.Contours = Contours;
        end
    end
end


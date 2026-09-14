classdef BadModule < seqcraft.Module
    %BADMODULE Fixture proving that Module rejects non-LogicBlock output.

    methods
        function obj = BadModule(opts)
            obj@seqcraft.Module(opts);
        end
    end

    methods (Access = protected)
        function value = buildImplicit(~, varargin) %#ok<INUSD>
            value = mr.makeDelay(1e-3);
        end
    end
end

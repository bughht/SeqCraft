classdef TestReadoutModule < seqcraft.Module
    %TESTREADOUTMODULE Minimal fixture for the MATLAB Module contract.

    properties
        gradient (1, 1) struct
    end

    methods
        function obj = TestReadoutModule(opts)
            obj@seqcraft.Module(opts);
            obj.gradient = mr.makeTrapezoid("x", opts, "Area", 20);
        end
    end

    methods (Access = protected)
        function block = buildImplicit(obj, varargin) %#ok<INUSD>
            block = seqcraft.LogicBlock();
            block.add(0, obj.gradient);
        end
    end
end

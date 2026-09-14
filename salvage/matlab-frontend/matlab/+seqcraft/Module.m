classdef (Abstract) Module < handle
    %MODULE A reusable component that builds one LogicBlock.

    properties
        opts (1, 1) struct
        tag (1, 1) string = ""
    end

    methods
        function obj = Module(opts, options)
            arguments
                opts (1, 1) struct
                options.Tag (1, 1) string = ""
            end
            obj.opts = opts;
            obj.tag = options.Tag;
        end

        function block = build(obj, varargin)
            %BUILD Assemble, validate, and name this module's LogicBlock.
            block = obj.buildImplicit(varargin{:});
            if ~isa(block, "seqcraft.LogicBlock")
                error("seqcraft:InvalidModuleOutput", ...
                    "%s.buildImplicit returned %s, not a seqcraft.LogicBlock.", ...
                    class(obj), class(block));
            end
            if strlength(block.tag) == 0
                if strlength(obj.tag) > 0
                    block.tag = obj.tag;
                else
                    parts = split(string(class(obj)), ".");
                    block.tag = parts(end);
                end
            end
        end
    end

    methods (Abstract, Access = protected)
        block = buildImplicit(obj, varargin)
    end
end

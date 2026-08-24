classdef LogicBlock < handle
    %LOGICBLOCK A tree of Pulseq events with relative start times.

    properties
        tag (1, 1) string = ""
        nodes (1, :) cell = {}
    end

    properties (Dependent, SetAccess = private)
        duration (1, 1) double
    end

    methods
        function obj = LogicBlock(tag)
            if nargin > 0
                obj.tag = string(tag);
            end
        end

        function obj = add(obj, startSeconds, varargin)
            %ADD Append Pulseq events or nested blocks without sorting them.
            validateattributes(startSeconds, {'numeric'}, {'scalar', 'real', 'finite'});
            for index = 1:numel(varargin)
                item = varargin{index};
                isEvent = isstruct(item) && isscalar(item) && isfield(item, "type");
                if ~isa(item, "seqcraft.LogicBlock") && ~isEvent
                    error("seqcraft:InvalidItem", ...
                        "LogicBlock.add accepts mr.make* events or nested LogicBlocks, got %s.", ...
                        class(item));
                end
                obj.nodes{end + 1} = struct("start", double(startSeconds), "item", item);
            end
        end

        function value = get.duration(obj)
            if isempty(obj.nodes)
                value = 0;
                return
            end
            ends = zeros(1, numel(obj.nodes));
            for index = 1:numel(obj.nodes)
                node = obj.nodes{index};
                if isa(node.item, "seqcraft.LogicBlock")
                    itemDuration = node.item.duration;
                elseif string(node.item.type) == "seqcraft_barrier"
                    itemDuration = 0;
                else
                    itemDuration = mr.calcDuration(node.item);
                end
                ends(index) = node.start + itemDuration;
            end
            value = max(ends);
        end

        function out = copy(obj)
            %COPY Copy this block and its node container, sharing child items.
            out = seqcraft.LogicBlock(obj.tag);
            out.nodes = obj.nodes;
        end
    end
end

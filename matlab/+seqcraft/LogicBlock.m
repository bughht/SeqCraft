classdef LogicBlock
    %LOGICBLOCK MATLAB value-object builder for an LBTX tree.

    properties (SetAccess = private)
        Tag (1, 1) string = ""
        Nodes (1, :) cell = {}
    end

    methods
        function obj = LogicBlock(tag)
            if nargin > 0
                obj.Tag = string(tag);
            end
        end

        function obj = add(obj, startSeconds, varargin)
            %ADD Append items at one start time, preserving insertion order.
            start = seqcraft.decimal(startSeconds);
            for index = 1:numel(varargin)
                node = struct();
                node.start_s = start;
                node.item = seqcraft.encodeItem(varargin{index});
                obj.Nodes{end + 1} = node;
            end
        end

        function value = toStruct(obj)
            %TOSTRUCT Return the JSON-ready LBTX block payload.
            value = struct();
            value.tag = char(obj.Tag);
            value.nodes = obj.Nodes;
        end
    end

    methods (Static)
        function obj = fromStruct(value)
            %FROMSTRUCT Reconstruct a builder from a decoded LBTX block.
            obj = seqcraft.LogicBlock(string(value.tag));
            nodes = value.nodes;
            if isempty(nodes)
                return
            end
            if isstruct(nodes)
                nodes = num2cell(nodes);
            end
            for index = 1:numel(nodes)
                node = nodes{index};
                item = node.item;
                if string(item.kind) == "block"
                    decoded = seqcraft.LogicBlock.fromStruct(item.block);
                else
                    decoded = item;
                end
                obj = obj.add(str2double(string(node.start_s)), decoded);
            end
        end
    end
end

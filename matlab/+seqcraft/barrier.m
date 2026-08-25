function value = barrier(tag)
%BARRIER Return a zero-duration SeqCraft block-boundary marker.
arguments
    tag (1, 1) string = "barrier"
end
value = struct("type", "seqcraft_barrier", "tag", tag, "delay", 0);
end

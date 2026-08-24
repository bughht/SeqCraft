function value = barrier(tag)
%BARRIER Construct an LBTX SeqCraft boundary marker.
arguments
    tag (1, 1) string = "barrier"
end
value = seqcraft.event("seqcraft_barrier", struct("tag", tag));
end

function value = delay(durationSeconds)
%DELAY Construct an LBTX delay event.
arguments
    durationSeconds (1, 1) double {mustBeNonnegative}
end
value = seqcraft.event("delay", struct("duration_s", seqcraft.decimal(durationSeconds)));
end

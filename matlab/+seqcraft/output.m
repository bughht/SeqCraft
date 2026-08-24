function value = output(channel, durationSeconds, options)
%OUTPUT Construct an LBTX digital output event.
arguments
    channel (1, 1) string
    durationSeconds (1, 1) double {mustBeNonnegative}
    options.DelaySeconds (1, 1) double {mustBeNonnegative} = 0
end
payload = struct("channel", channel, "delay_s", seqcraft.decimal(options.DelaySeconds), ...
    "duration_s", seqcraft.decimal(durationSeconds));
value = seqcraft.event("output", payload);
end

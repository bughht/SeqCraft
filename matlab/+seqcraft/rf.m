function value = rf(signalHz, sampleTimesSeconds, options)
%RF Construct an LBTX RF event with explicit complex samples.
arguments
    signalHz (1, :) double {mustBeFinite}
    sampleTimesSeconds (1, :) double {mustBeNonnegative}
    options.ShapeDurationSeconds (1, 1) double {mustBeNonnegative}
    options.FrequencyOffsetHz (1, 1) double {mustBeFinite} = 0
    options.PhaseOffsetRad (1, 1) double {mustBeFinite} = 0
    options.FrequencyOffsetPpm (1, 1) double {mustBeFinite} = 0
    options.PhaseOffsetPpm (1, 1) double {mustBeFinite} = 0
    options.DeadTimeSeconds (1, 1) double {mustBeNonnegative} = 0
    options.RingdownTimeSeconds (1, 1) double {mustBeNonnegative} = 0
    options.DelaySeconds (1, 1) double {mustBeNonnegative} = 0
    options.CenterSeconds (1, 1) double {mustBeNonnegative}
    options.Use = []
end
if numel(signalHz) ~= numel(sampleTimesSeconds) || isempty(signalHz)
    error("seqcraft:InvalidRF", "signalHz and sampleTimesSeconds must have equal non-zero lengths.");
end
use = options.Use;
if isempty(use)
    use = string(missing);
else
    use = char(string(use));
end
payload = struct( ...
    "signal_hz", struct("real", real(signalHz), "imag", imag(signalHz)), ...
    "sample_times_s", string(arrayfun(@seqcraft.decimal, sampleTimesSeconds, UniformOutput=false)), ...
    "shape_duration_s", seqcraft.decimal(options.ShapeDurationSeconds), ...
    "frequency_offset_hz", options.FrequencyOffsetHz, ...
    "phase_offset_rad", options.PhaseOffsetRad, ...
    "frequency_offset_ppm", options.FrequencyOffsetPpm, ...
    "phase_offset_ppm", options.PhaseOffsetPpm, ...
    "dead_time_s", seqcraft.decimal(options.DeadTimeSeconds), ...
    "ringdown_time_s", seqcraft.decimal(options.RingdownTimeSeconds), ...
    "delay_s", seqcraft.decimal(options.DelaySeconds), ...
    "center_s", seqcraft.decimal(options.CenterSeconds), ...
    "use", use);
value = seqcraft.event("rf", payload);
end

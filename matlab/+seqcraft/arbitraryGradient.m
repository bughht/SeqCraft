function value = arbitraryGradient(channel, waveformHzPerM, sampleTimesSeconds, options)
%ARBITRARYGRADIENT Construct an arbitrary or extended LBTX gradient.
arguments
    channel (1, 1) string {mustBeMember(channel, ["x", "y", "z"])}
    waveformHzPerM (1, :) double {mustBeFinite}
    sampleTimesSeconds (1, :) double {mustBeNonnegative}
    options.DelaySeconds (1, 1) double {mustBeNonnegative} = 0
    options.ShapeDurationSeconds (1, 1) double = NaN
    options.AreaPerM (1, 1) double = NaN
    options.FirstHzPerM (1, 1) double {mustBeFinite} = 0
    options.LastHzPerM (1, 1) double {mustBeFinite} = 0
end
if numel(waveformHzPerM) ~= numel(sampleTimesSeconds) || isempty(waveformHzPerM)
    error("seqcraft:InvalidGradient", ...
        "waveformHzPerM and sampleTimesSeconds must have equal non-zero lengths.");
end
shapeDuration = options.ShapeDurationSeconds;
if isnan(shapeDuration)
    if numel(sampleTimesSeconds) == 1
        error("seqcraft:InvalidGradient", "ShapeDurationSeconds is required for one sample.");
    end
    shapeDuration = sampleTimesSeconds(end) + (sampleTimesSeconds(end) - sampleTimesSeconds(end-1)) / 2;
elseif ~isfinite(shapeDuration) || shapeDuration < 0
    error("seqcraft:InvalidGradient", "ShapeDurationSeconds must be finite and non-negative.");
end
area = options.AreaPerM;
if isnan(area)
    knotsT = [0, sampleTimesSeconds, shapeDuration];
    knotsG = [options.FirstHzPerM, waveformHzPerM, options.LastHzPerM];
    area = trapz(knotsT, knotsG);
elseif ~isfinite(area)
    error("seqcraft:InvalidGradient", "AreaPerM must be finite.");
end
payload = struct( ...
    "channel", channel, ...
    "waveform_hz_per_m", waveformHzPerM, ...
    "sample_times_s", string(arrayfun(@seqcraft.decimal, sampleTimesSeconds, UniformOutput=false)), ...
    "delay_s", seqcraft.decimal(options.DelaySeconds), ...
    "shape_duration_s", seqcraft.decimal(shapeDuration), ...
    "area_per_m", area, ...
    "first_hz_per_m", options.FirstHzPerM, ...
    "last_hz_per_m", options.LastHzPerM);
value = seqcraft.event("grad", payload);
end

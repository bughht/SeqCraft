function value = trapezoid(channel, amplitudeHzPerM, riseSeconds, flatSeconds, fallSeconds, options)
%TRAPEZOID Construct an LBTX trapezoid from canonical SI fields.
arguments
    channel (1, 1) string {mustBeMember(channel, ["x", "y", "z"])}
    amplitudeHzPerM (1, 1) double {mustBeFinite}
    riseSeconds (1, 1) double {mustBeNonnegative}
    flatSeconds (1, 1) double {mustBeNonnegative}
    fallSeconds (1, 1) double {mustBeNonnegative}
    options.DelaySeconds (1, 1) double {mustBeNonnegative} = 0
    options.FirstHzPerM (1, 1) double {mustBeFinite} = 0
    options.LastHzPerM (1, 1) double {mustBeFinite} = 0
    options.AreaPerM (1, 1) double = NaN
    options.FlatAreaPerM (1, 1) double = NaN
end
area = options.AreaPerM;
if isnan(area)
    area = amplitudeHzPerM * (flatSeconds + (riseSeconds + fallSeconds) / 2);
elseif ~isfinite(area)
    error("seqcraft:InvalidGradient", "AreaPerM must be finite.");
end
flatArea = options.FlatAreaPerM;
if isnan(flatArea)
    flatArea = amplitudeHzPerM * flatSeconds;
elseif ~isfinite(flatArea)
    error("seqcraft:InvalidGradient", "FlatAreaPerM must be finite.");
end
payload = struct( ...
    "channel", channel, ...
    "amplitude_hz_per_m", amplitudeHzPerM, ...
    "rise_time_s", seqcraft.decimal(riseSeconds), ...
    "flat_time_s", seqcraft.decimal(flatSeconds), ...
    "fall_time_s", seqcraft.decimal(fallSeconds), ...
    "area_per_m", area, ...
    "flat_area_per_m", flatArea, ...
    "delay_s", seqcraft.decimal(options.DelaySeconds), ...
    "first_hz_per_m", options.FirstHzPerM, ...
    "last_hz_per_m", options.LastHzPerM);
value = seqcraft.event("trap", payload);
end

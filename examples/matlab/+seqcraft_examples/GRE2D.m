classdef GRE2D < seqcraft.Module
    %GRE2D Small example module built entirely from native Pulseq events.

    properties (SetAccess = private)
        matrix (1, 2) double
        fovM (1, 2) double
        sliceThicknessM (1, 1) double
        trSeconds (1, 1) double
        rf (1, 1) struct
        gz (1, 1) struct
        gzRephase (1, 1) struct
        gxPre (1, 1) struct
        phaseEncodes (1, :) cell
        lineLabels (1, :) cell
        gx (1, 1) struct
        adc (1, 1) struct
        gxSpoil (1, 1) struct
        gzSpoil (1, 1) struct
        trDelay (1, 1) struct
        excitationSeconds (1, 1) double
        readoutStart (1, 1) double
        spoilStart (1, 1) double
        contentSeconds (1, 1) double
    end

    methods
        function obj = GRE2D(opts, options)
            arguments
                opts (1, 1) struct
                options.Matrix (1, 2) double = [16 8]
                options.FOVM (1, 2) double = [0.22 0.22]
                options.SliceThicknessM (1, 1) double = 5e-3
                options.FlipAngleRad (1, 1) double = deg2rad(10)
                options.TRSeconds (1, 1) double = 20e-3
                options.DwellSeconds (1, 1) double = 20e-6
                options.Tag (1, 1) string = "gre_2d"
            end
            obj@seqcraft.Module(opts, Tag=options.Tag);

            validateattributes(options.Matrix, {'numeric'}, ...
                {'integer', 'positive', 'finite'});
            validateattributes(options.FOVM, {'numeric'}, {'positive', 'finite'});
            validateattributes(options.SliceThicknessM, {'numeric'}, {'positive', 'finite'});
            validateattributes(options.TRSeconds, {'numeric'}, {'positive', 'finite'});
            validateattributes(options.DwellSeconds, {'numeric'}, {'positive', 'finite'});

            obj.matrix = options.Matrix;
            obj.fovM = options.FOVM;
            obj.sliceThicknessM = options.SliceThicknessM;
            obj.trSeconds = options.TRSeconds;

            % The constructor designs events; buildImplicit only assembles them.
            [obj.rf, obj.gz, gzRephaseMinimum] = mr.makeSincPulse( ...
                options.FlipAngleRad, opts, ...
                "Duration", 1e-3, ...
                "SliceThickness", obj.sliceThicknessM, ...
                "use", "excitation");

            readoutArea = obj.matrix(1) / obj.fovM(1);
            obj.gx = mr.makeTrapezoid( ...
                "x", opts, ...
                "FlatArea", readoutArea, ...
                "FlatTime", obj.matrix(1) * options.DwellSeconds);
            obj.adc = mr.makeAdc( ...
                obj.matrix(1), opts, ...
                "Dwell", options.DwellSeconds, ...
                "Delay", obj.gx.riseTime);

            winderMaxGrad = opts.maxGrad / sqrt(3);
            winderMaxSlew = opts.maxSlew / sqrt(3);
            gzRephaseConstrainedMinimum = mr.makeTrapezoid( ...
                "z", opts, "Area", gzRephaseMinimum.area, ...
                "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
            gxPreMinimum = mr.makeTrapezoid( ...
                "x", opts, "Area", -obj.gx.area / 2, ...
                "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
            centerLine = floor(obj.matrix(2) / 2);
            gyMaximum = mr.makeTrapezoid( ...
                "y", opts, "Area", centerLine / obj.fovM(2), ...
                "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
            winderSeconds = max([ ...
                mr.calcDuration(gzRephaseConstrainedMinimum), ...
                mr.calcDuration(gxPreMinimum), ...
                mr.calcDuration(gyMaximum)]);
            winderSeconds = ceil(winderSeconds / opts.gradRasterTime) * opts.gradRasterTime;

            obj.gzRephase = mr.makeTrapezoid( ...
                "z", opts, ...
                "Area", gzRephaseMinimum.area, ...
                "Duration", winderSeconds, ...
                "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
            obj.gxPre = mr.makeTrapezoid( ...
                "x", opts, "Area", -obj.gx.area / 2, "Duration", winderSeconds, ...
                "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);

            obj.phaseEncodes = cell(1, obj.matrix(2));
            obj.lineLabels = cell(1, obj.matrix(2));
            for lineIndex = 0:(obj.matrix(2) - 1)
                obj.phaseEncodes{lineIndex + 1} = mr.makeTrapezoid( ...
                    "y", opts, ...
                    "Area", (lineIndex - centerLine) / obj.fovM(2), ...
                    "Duration", winderSeconds, ...
                    "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
                obj.lineLabels{lineIndex + 1} = mr.makeLabel("SET", "LIN", lineIndex);
            end

            spoilMaxGrad = opts.maxGrad / sqrt(2);
            spoilMaxSlew = opts.maxSlew / sqrt(2);
            gxSpoilMinimum = mr.makeTrapezoid( ...
                "x", opts, "Area", 2 * readoutArea, ...
                "maxGrad", spoilMaxGrad, "maxSlew", spoilMaxSlew);
            gzSpoilMinimum = mr.makeTrapezoid( ...
                "z", opts, "Area", 2 / obj.sliceThicknessM, ...
                "maxGrad", spoilMaxGrad, "maxSlew", spoilMaxSlew);
            spoilSeconds = max( ...
                mr.calcDuration(gxSpoilMinimum), mr.calcDuration(gzSpoilMinimum));
            obj.gxSpoil = mr.makeTrapezoid( ...
                "x", opts, "Area", 2 * readoutArea, "Duration", spoilSeconds, ...
                "maxGrad", spoilMaxGrad, "maxSlew", spoilMaxSlew);
            obj.gzSpoil = mr.makeTrapezoid( ...
                "z", opts, "Area", 2 / obj.sliceThicknessM, "Duration", spoilSeconds, ...
                "maxGrad", spoilMaxGrad, "maxSlew", spoilMaxSlew);

            obj.excitationSeconds = mr.calcDuration(obj.rf, obj.gz);
            readoutSeconds = mr.calcDuration(obj.gx, obj.adc);
            obj.readoutStart = obj.excitationSeconds + winderSeconds;
            obj.spoilStart = obj.readoutStart + readoutSeconds;
            obj.contentSeconds = obj.spoilStart + spoilSeconds;
            if obj.contentSeconds >= obj.trSeconds
                error("seqcraft_examples:TRTooShort", ...
                    "GRE line duration %.3f ms does not fit in TR %.3f ms.", ...
                    obj.contentSeconds * 1e3, obj.trSeconds * 1e3);
            end
            obj.trDelay = mr.makeDelay(obj.trSeconds - obj.contentSeconds);
        end
    end

    methods (Access = protected)
        function root = buildImplicit(obj, varargin)
            root = seqcraft.LogicBlock();
            for lineIndex = 0:(obj.matrix(2) - 1)
                line = seqcraft.LogicBlock("line_" + lineIndex);
                line.add(0, obj.rf, obj.gz);
                line.add( ...
                    obj.excitationSeconds, ...
                    obj.gzRephase, obj.gxPre, obj.phaseEncodes{lineIndex + 1});
                line.add( ...
                    obj.readoutStart, ...
                    obj.gx, obj.adc, obj.lineLabels{lineIndex + 1});
                line.add(obj.spoilStart, obj.gxSpoil, obj.gzSpoil);
                line.add(obj.contentSeconds, obj.trDelay);
                root.add(lineIndex * obj.trSeconds, line);
            end
        end
    end
end

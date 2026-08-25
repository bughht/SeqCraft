classdef GRE2DTR < seqcraft.Module
    %GRE2DTR One example GRE 2D repetition built from native Pulseq events.

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
        phaseRewinds (1, :) cell
        lineLabels (1, :) cell
        gx (1, 1) struct
        adc (1, 1) struct
        gxSpoil (1, 1) struct
        gzSpoil (1, 1) struct
        trDelay (1, 1) struct
        excitationSeconds (1, 1) double
        readoutStart (1, 1) double
        tailStart (1, 1) double
        contentSeconds (1, 1) double
    end

    methods
        function obj = GRE2DTR(opts, options)
            arguments
                opts (1, 1) struct
                options.Matrix (1, 2) double = [16 8]
                options.FOVM (1, 2) double = [0.22 0.22]
                options.SliceThicknessM (1, 1) double = 5e-3
                options.FlipAngleRad (1, 1) double = deg2rad(10)
                options.TRSeconds (1, 1) double = 20e-3
                options.DwellSeconds (1, 1) double = 20e-6
                options.Tag (1, 1) string = ""
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

            % The constructor designs events; buildImplicit only assembles one requested TR.
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

            preEchoSample = floor(obj.matrix(1) / 2);
            flatBeforeEcho = ...
                obj.adc.delay - obj.gx.riseTime + ...
                (preEchoSample + 0.5) * options.DwellSeconds;
            prephaserArea = -obj.gx.amplitude * ...
                (0.5 * obj.gx.riseTime + flatBeforeEcho);

            winderMaxGrad = opts.maxGrad / sqrt(3);
            winderMaxSlew = opts.maxSlew / sqrt(3);
            gzRephaseConstrainedMinimum = mr.makeTrapezoid( ...
                "z", opts, "Area", gzRephaseMinimum.area, ...
                "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
            gxPreMinimum = mr.makeTrapezoid( ...
                "x", opts, "Area", prephaserArea, ...
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
                "z", opts, "Area", gzRephaseMinimum.area, "Duration", winderSeconds, ...
                "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
            obj.gxPre = mr.makeTrapezoid( ...
                "x", opts, "Area", prephaserArea, "Duration", winderSeconds, ...
                "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);

            tailMaxGrad = opts.maxGrad / sqrt(3);
            tailMaxSlew = opts.maxSlew / sqrt(3);
            gyRewindMaximum = mr.makeTrapezoid( ...
                "y", opts, "Area", centerLine / obj.fovM(2), ...
                "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
            gxSpoilMinimum = mr.makeTrapezoid( ...
                "x", opts, "Area", 2 * readoutArea, ...
                "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
            gzSpoilMinimum = mr.makeTrapezoid( ...
                "z", opts, "Area", 2 / obj.sliceThicknessM, ...
                "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
            tailSeconds = max([ ...
                mr.calcDuration(gyRewindMaximum), ...
                mr.calcDuration(gxSpoilMinimum), ...
                mr.calcDuration(gzSpoilMinimum)]);
            tailSeconds = ceil(tailSeconds / opts.gradRasterTime) * opts.gradRasterTime;

            obj.phaseEncodes = cell(1, obj.matrix(2));
            obj.phaseRewinds = cell(1, obj.matrix(2));
            obj.lineLabels = cell(1, obj.matrix(2));
            for lineIndex = 0:(obj.matrix(2) - 1)
                phaseArea = (lineIndex - centerLine) / obj.fovM(2);
                obj.phaseEncodes{lineIndex + 1} = mr.makeTrapezoid( ...
                    "y", opts, "Area", phaseArea, "Duration", winderSeconds, ...
                    "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
                obj.phaseRewinds{lineIndex + 1} = mr.makeTrapezoid( ...
                    "y", opts, "Area", -phaseArea, "Duration", tailSeconds, ...
                    "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
                obj.lineLabels{lineIndex + 1} = mr.makeLabel("SET", "LIN", lineIndex);
            end

            obj.gxSpoil = mr.makeTrapezoid( ...
                "x", opts, "Area", 2 * readoutArea, "Duration", tailSeconds, ...
                "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
            obj.gzSpoil = mr.makeTrapezoid( ...
                "z", opts, "Area", 2 / obj.sliceThicknessM, "Duration", tailSeconds, ...
                "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);

            obj.excitationSeconds = mr.calcDuration(obj.rf, obj.gz);
            obj.readoutStart = obj.excitationSeconds + winderSeconds;
            obj.tailStart = obj.readoutStart + mr.calcDuration(obj.gx, obj.adc);
            obj.contentSeconds = obj.tailStart + tailSeconds;
            if obj.contentSeconds >= obj.trSeconds
                error("seqcraft_examples:TRTooShort", ...
                    "GRE repetition %.3f ms does not fit in TR %.3f ms.", ...
                    obj.contentSeconds * 1e3, obj.trSeconds * 1e3);
            end
            obj.trDelay = mr.makeDelay(obj.trSeconds - obj.contentSeconds);
        end
    end

    methods (Access = protected)
        function line = buildImplicit(obj, lineIndex)
            validateattributes(lineIndex, {'numeric'}, {'scalar', 'integer', 'finite'});
            if lineIndex < 0 || lineIndex >= obj.matrix(2)
                error("seqcraft_examples:InvalidLine", ...
                    "lineIndex must be between 0 and %d, got %d.", ...
                    obj.matrix(2) - 1, lineIndex);
            end

            offset = lineIndex + 1;
            line = seqcraft.LogicBlock();
            line.add(0, obj.rf, obj.gz);
            line.add( ...
                obj.excitationSeconds, ...
                obj.gzRephase, obj.gxPre, obj.phaseEncodes{offset});
            line.add(obj.readoutStart, obj.gx, obj.adc, obj.lineLabels{offset});
            line.add(obj.tailStart, obj.phaseRewinds{offset}, obj.gxSpoil, obj.gzSpoil);
            line.add(obj.contentSeconds, obj.trDelay);
        end
    end
end

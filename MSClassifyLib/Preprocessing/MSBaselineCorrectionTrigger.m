classdef MSBaselineCorrectionTrigger<MSPreprocessingTrigger
    methods
        function apply(~,maldiData)
            % Applies baseline correction to the maldiData
            data = abs(msbackadj(maldiData.mzVector', maldiData.data'))';
            normalizat=maldiData.normalization;
            if ~strcmp(normalizat,'raw');
                maldiData.setNormalization('raw')
            end
            maldiData.data=data;
            % recalculate dataNorm
            maldiData.initNormalization;
            if ~strcmp(normalizat,'raw');
                maldiData.setNormalization(normalizat);
            end
        end
    end
end

function paramTranslated = MSTranslateDecEquivInputs(param)

  paramTranslated = param;
  if isfield(param,'rescaling')
    paramTranslated.rescaling = MSDecomposer.translateRescalingInputToChar(param.rescaling);
  end

end

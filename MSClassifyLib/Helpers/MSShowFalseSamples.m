function MSShowFalseSamples (prediction, reference, samples)
  % Show misclassified items in misclassified samples
  % Display an image where only samples misclassified after aggregation are
  % shown. Within these samples, correctly classified items are shown in
  % gray, whereas misclassified items are shown in color.

  % Compute prediction and reference on sample level
  samplePredict = prediction.aggregate(samples);
  sampleReference = reference.aggregate(samples, 'unique');
  % Identify misclassified samples
  sampleMask = samplePredict.data > 0 & samplePredict.data ~= sampleReference.data;
  % Generate predicted item labels for misclassified samples
  sampleItemsMask = samples.data > 0;
  falseSampleItems = zeros(prediction.numItems,1);
  falseSampleItems(sampleItemsMask) = prediction.data(sampleItemsMask) .* ...
                                     sampleMask(samples.data(sampleItemsMask));
  % Add a label for correctly classified items and replace resp. label values
  correctItems = falseSampleItems > 0 & falseSampleItems == reference.data;
  falseSampleItems(correctItems) = prediction.numLabels+1;
  % Show items of misclassified samples
  im = prediction.positions.encube(falseSampleItems);
  figure
  imagesc(flipud(im))
  colormap([ones(1,3); linspecer(prediction.numLabels); 0.5*ones(1,3)])
end


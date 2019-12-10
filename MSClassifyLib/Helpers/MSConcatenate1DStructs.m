function sout = MSConcatenate1DStructs(cellArrayStructs)

if length(cellArrayStructs) == 1;
  sout = cellArrayStructs{1};
  return;
end

sM = cellArrayStructs{1};
for k = 2:length(cellArrayStructs)
  sM = [sM cellArrayStructs{k}];
end

sout = sM(1);
fn = fieldnames(sM);
for l = 1:length(fn)
  for k = 1:length(sM)
    sout.(fn{l}) = [sout.(fn{l}) sM(k).(fn{l})];
  end
end

end
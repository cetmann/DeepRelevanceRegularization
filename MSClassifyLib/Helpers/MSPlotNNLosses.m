function MSPlotNNLosses(nnModel)
    [NLL,l1,l2,l1l2,logitSens,logitSqSens,l1l2LogitSens,probSens,lossSens] = MSReadNNLogFile(nnModel);
    figure
    plot(NLL)
    legend('NLL')
    
    figure
    plot(l1)
    hold on
    plot(l2)
    legend('l1 (weights)','l2^2 (weights)')
    
    figure
    plot(l1l2)
    legend('l1 / l2')
    
    figure
    plot(logitSens)
    hold on
    plot(logitSqSens)
    hold on
    plot(probSens)
    hold on
    plot(lossSens)
    legend('logit saliency (l1)', 'logit saliency (l2^2)',...
        'probability saliency','NLL saliency')

    figure
    plot(l1l2LogitSens)
    legend('logit saliency (l1/l2)')
end

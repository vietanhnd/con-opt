% Plotting iterations on himmelblau 
% fmincon-iter

% HimmelblauPlotIter('data/fmincon-iter/results_sqp.csv', false)
% % save as png to figures folder
% saveas(gcf, 'figures/hb-fmincon-iter.png')

% % 'local' SQP

% HimmelblauPlotIter('data/Himmelblau/eqcon/SQP_HB_BFGS_local.csv', true)
% % % save as png to figures folder
% saveas(gcf, 'figures/Himmelblau/hb-SQP-local.png')

% HimmelblauPlotIter('data/Himmelblau/noeqcon/SQP_HB_BFGS_local.csv', false)
% % % save as png to figures folder
% saveas(gcf, 'figures/Himmelblau/hb-SQP-local-noeqcon.png')

% 'line-search' SQP

HimmelblauPlotIter('data/Himmelblau/eqcon/SQP_HB_BFGS_line-search.csv', true)
% % save as png to figures folder
saveas(gcf, 'figures/Himmelblau/hb-SQP-line-search.png')
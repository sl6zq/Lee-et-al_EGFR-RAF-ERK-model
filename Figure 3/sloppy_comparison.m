function [norm_loadings, data_matlab, data_sloppy_sameorder] = sloppy_comparison(number, sampling, ncomp, select_ncomp, wanted_param, Xloadings, Q2Y, PLSR_cat)

wanted_sloppy = load('V','V'); %Select the sign of sloppiness analysis eigenvalues (V preferred)
data_sloppy = struct2array(wanted_sloppy);

% SA labels
mylabels;

% Sloppiness labels
labels = readtable('param-order.txt','Format','auto');
labels2 = table2array(labels); 
labels3 = string(labels2); %for reference

% Make sloppycell order follow my order
[data_sloppy_sameorder] = order_data(data_sloppy);

% Scale loadings from -0.6 to 0.6
norm_loadings = zeros(length(wanted_param), ncomp);
for m=1:ncomp
    norm_loadings(:,m) = rescale(Xloadings(:,m),-0.6,0.6);
end
%{
figure
h = heatmap(select_ncomp,PLSR_cat,norm_loadings(:,1:6),'CellLabelColor','none','GridVisible','off');
bwr = @(n)interp1([1 2 3], [0 0 1; 1 1 1; 1 0 0], linspace(1, 3, n), 'linear');
colormap(bwr(200));
h.xlabel('Principal Component');
ax = gca;
set(gca,'FontSize',10);
%}
normloadings_save = norm_loadings(:,1:6);
data_matlab = normloadings_save;

combine_matlab_sloppy = [norm_loadings(:,1:6) data_sloppy_sameorder];
combine_xvalues1 = {['1']; ['(Q^2Y = ' num2str(round(Q2Y(:,1),2)) ')']};
combine_xvalues2 = {['2']; ['(Q^2Y = ' num2str(round(Q2Y(:,2),2)) ')']};
combine_xvalues3 = {['3']; ['(Q^2Y = ' num2str(round(Q2Y(:,3),2)) ')']};
combine_xvalues4 = {['4']; ['(Q^2Y = ' num2str(round(Q2Y(:,4),2)) ')']};
combine_xvalues5 = {['5']; ['(Q^2Y = ' num2str(round(Q2Y(:,5),2)) ')']};
combine_xvalues6 = {['6']; ['(Q^2Y = ' num2str(round(Q2Y(:,6),2)) ')']};

combine_xvalues1 = ['1','(Q^2Y = ' num2str(round(Q2Y(:,1),2)) ')'];
combine_xvalues2 = ['2','(Q^2Y = ' num2str(round(Q2Y(:,2),2)) ')'];
combine_xvalues3 = ['3','(Q^2Y = ' num2str(round(Q2Y(:,3),2)) ')'];
combine_xvalues4 = ['4','(Q^2Y = ' num2str(round(Q2Y(:,4),2)) ')'];
combine_xvalues5 = ['5','(Q^2Y = ' num2str(round(Q2Y(:,5),2)) ')'];
combine_xvalues6 = ['6','(Q^2Y = ' num2str(round(Q2Y(:,6),2)) ')'];



combine_xvalues = {'1','2','3','4','5','6','0','1','2','3','4','5'};

%combine_xvalues = vertcat(combine_xvalues1, combine_xvalues2, combine_xvalues3, combine_xvalues4, combine_xvalues5, combine_xvalues6);

h1 = figure;
h2 = heatmap(combine_matlab_sloppy,'CellLabelColor','none','GridVisible','off');
h2.XDisplayLabels = combine_xvalues;
h2.YDisplayLabels = PLSR_cat;
h2.XLabel = 'PC                                                                     Mode';
h1.Position = [100 30 900 600];
bwr = @(n)interp1([1 2 3], [0 0 1; 1 1 1; 1 0 0], linspace(1, 3, n), 'linear');
colormap(bwr(200)); 
set(gca,'FontSize',9);
h2.FontSize = 6;
t1 = {['Parameter SA using ', num2str(number), ' ', sampling, ' parameter sets and Sloppiness Analysis']}; 
filename_combine = ['matlab_sloppy_heatmap ', num2str(number), sampling, '.fig'];
filename_combine2 = ['matlab_sloppy_heatmap ', num2str(number), sampling '.jpeg'];
%saveas(h1, filename_combine);
%saveas(h1, filename_combine2);


h4 = figure;
h3 = heatmap(combine_matlab_sloppy', 'CellLabelColor','none','GridVisible','off');
h3.XDisplayLabels = PLSR_cat;
h3.YDisplayLabels = combine_xvalues;
%h2.XLabel = 'PC                                                                     Mode';
h4.Position = [100 100 900 400];
bwr = @(n)interp1([1 2 3], [0 0 1; 1 1 1; 1 0 0], linspace(1, 3, n), 'linear');
colormap(bwr(200)); 
set(gca,'FontSize',8);
h3.FontSize = 8;
t1 = {['Parameter SA using ', num2str(number), ' ', sampling, ' parameter sets and Sloppiness Analysis']}; 
filename_combine = ['matlab_sloppy_heatmap_flip ', num2str(number), sampling, '.fig'];
filename_combine2 = ['matlab_sloppy_heatmap_flip ', num2str(number), sampling '.jpeg'];




end

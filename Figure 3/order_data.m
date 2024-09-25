function [data_sloppy_sameorder] = order_data(data_sloppy, sloppy_labels, PLSRSA_labels)
% Puts the parameters from Sloppiness analysis in the same order as PLSR SA
data_sloppy_sameorder          = zeros(size(data_sloppy,1), size(data_sloppy,2));    

for i=1:40
    [~, idx]                   = ismember(PLSRSA_labels, sloppy_labels);
    idx(i)
    data_sloppy_sameorder(i,:) = data_sloppy(idx(i),:);
end

end

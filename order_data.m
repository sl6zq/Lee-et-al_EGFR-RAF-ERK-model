function [data_sloppy_sameorder] = order_data(data_sloppy, sloppy_labels, PLSRSA_labels)

%     data_sloppy_sameorder(1,:) = data_sloppy(22,:);
%     data_sloppy_sameorder(2,:) = data_sloppy(26,:);
%     data_sloppy_sameorder(3,:) = data_sloppy(13,:);
%     data_sloppy_sameorder(4,:) = data_sloppy(36,:);
%     data_sloppy_sameorder(5,:) = data_sloppy(8,:);
%     data_sloppy_sameorder(6,:) = data_sloppy(7,:);
%     data_sloppy_sameorder(7,:) = data_sloppy(15,:);
%     data_sloppy_sameorder(8,:) = data_sloppy(6,:);
%     data_sloppy_sameorder(9,:) = data_sloppy(5,:);
%     data_sloppy_sameorder(10,:) = data_sloppy(21,:);
%     data_sloppy_sameorder(11,:) = data_sloppy(3,:);
%     data_sloppy_sameorder(12,:) = data_sloppy(34,:);
%     data_sloppy_sameorder(13,:) = data_sloppy(33,:);
%     data_sloppy_sameorder(14,:) = data_sloppy(28,:);
%     data_sloppy_sameorder(15,:) = data_sloppy(1,:);
%     data_sloppy_sameorder(16,:) = data_sloppy(38,:);
%     data_sloppy_sameorder(17,:) = data_sloppy(37,:);
%     data_sloppy_sameorder(18,:) = data_sloppy(12,:);
%     data_sloppy_sameorder(19,:) = data_sloppy(11,:);
%     data_sloppy_sameorder(20,:) = data_sloppy(2,:);
%     data_sloppy_sameorder(21,:) = data_sloppy(14,:);
%     data_sloppy_sameorder(22,:) = data_sloppy(30,:);
%     data_sloppy_sameorder(23,:) = data_sloppy(29,:);
%     data_sloppy_sameorder(24,:) = data_sloppy(16,:);
%     data_sloppy_sameorder(25,:) = data_sloppy(39,:);
%     data_sloppy_sameorder(26,:) = data_sloppy(9,:);
%     data_sloppy_sameorder(27,:) = data_sloppy(18,:);
%     data_sloppy_sameorder(28,:) = data_sloppy(4,:);
%     data_sloppy_sameorder(29,:) = data_sloppy(17,:);
%     data_sloppy_sameorder(30,:) = data_sloppy(20,:);
%     data_sloppy_sameorder(31,:) = data_sloppy(35,:);
%     data_sloppy_sameorder(32,:) = data_sloppy(27,:);
%     data_sloppy_sameorder(33,:) = data_sloppy(25,:);
%     data_sloppy_sameorder(34,:) = data_sloppy(24,:);
%     data_sloppy_sameorder(35,:) = data_sloppy(10,:);
%     data_sloppy_sameorder(36,:) = data_sloppy(32,:);
%     data_sloppy_sameorder(37,:) = data_sloppy(31,:);
%     data_sloppy_sameorder(38,:) = data_sloppy(40,:);
%     data_sloppy_sameorder(39,:) = data_sloppy(19,:);
%     data_sloppy_sameorder(40,:) = data_sloppy(23,:);
data_sloppy_sameorder          = zeros(size(data_sloppy,1), size(data_sloppy,2));    

for i=1:40
    [~, idx]                   = ismember(PLSRSA_labels, sloppy_labels);
    idx(i)
    data_sloppy_sameorder(i,:) = data_sloppy(idx(i),:);
end

end

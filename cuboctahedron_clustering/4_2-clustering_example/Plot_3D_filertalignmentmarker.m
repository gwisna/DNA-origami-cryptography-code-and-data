clc;
clear all;
close all;
data = readtable('decoded_data_global.csv');

decoded_data=table2array([data(:,1) data(:,2) data(:,3)]);
decoded = decoded_data(decoded_data(:,3)>=3,1:3)
z=zeros(8,16);
for i=1:1:8
    for j=1:1:16
        for k=1:1:size(decoded,1)
            if decoded(k,1)==j-1 && decoded(k,2)==i-1
                z(i,j)=z(i,j)+1
            end
        end
    end
end

b = bar3(z);
set(gca,'YTickLabel',[0 1 2 3 4 5 6 7])
set(gca,'XTickLabel',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15])
xlabel('Numbers');
ylabel('Position');
zlabel('Counts');
for i = 1:length(b)
  index = logical(kron(z(:, i) == 0, ones(6, 1)));
  zData = get(b(i), 'ZData');
  zData(index, :) = nan;
  set(b(i), 'ZData', zData);
end
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
colorbar

sorted_descend=sort((reshape(z.',1,[])),'descend');
maximum_SN_ratio=sorted_descend(1)/sorted_descend(2)
average_accuracy=sorted_descend(1)/sum(sorted_descend)
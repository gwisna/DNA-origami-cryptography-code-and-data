clear all;
clc;
close all;
%opening file
locs1 = h5read('4_2_cuboctahedron_DNA-PAINT.hdf5', '/locs');

%input
frame=locs1.frame;
x=117*locs1.x;
y=117*locs1.y;
photons=locs1.photons;
sx=locs1.sx;
sy=locs1.sy;
bg=locs1.bg;
lpx=locs1.lpx;
lpy=locs1.lpy;
ellipticity=locs1.ellipticity;
net_gradient=locs1.net_gradient;
z=locs1.z;
d_zcalib=locs1.d_zcalib;
group=locs1.group+1;  
for group_id=1:1:max(group)
pick_number=group_id;%input for pick number to be processed<----------------
%making folder based on pick number
folder_name=sprintf('Pick # %.f', pick_number);
mkdir (folder_name)
%creating proper matrix data
X=[x y z double(group)];
Y = X(X(:,4)==pick_number,1:3)
header_1 = {'x(nm)' 'y(nm)' 'z(nm)'};
Preprocessed_data=[header_1;num2cell([abs(-Y(:,1)) abs(-Y(:,2)) -Y(:,3)])]; %mirroring with respect to 0,0,0 and keeping x and y positive
%plotting pre-processed experimental data #Figure 1 and save it
fig1_name=sprintf('Pre-processed Data of Pick # %.f', pick_number);
fig1=figure('Name',fig1_name);
set(fig1,'visible','off');
scatter3(abs(-Y(:,1)), abs(-Y(:,2)), -Y(:,3), 15, 'filled');
title(fig1_name);
title(fig1_name)
xlabel('x');
ylabel('y');
zlabel('z');
pbaspect([1 1 1]);
%saving
saveas(fig1,fullfile(folder_name,fig1_name),'fig');
saveas(fig1,fullfile(folder_name,fig1_name),'pdf');
Preprocessed_data_name=sprintf('Preprocessed_data_pick_num%.f.csv', pick_number);
writecell(Preprocessed_data,fullfile(folder_name,Preprocessed_data_name));

cluster_trial=20%number of cluster to be tried for elbow method k-means
%clean up noise using dbscan
D = pdist2(Y,Y);
[idx_dbscan, corepts] = dbscan(D,7,4,'Distance','precomputed');
Y_new=Y(idx_dbscan~=-1,1:3);
idx_dbscan_new=idx_dbscan(idx_dbscan~=-1,1);
%creating proper matrix data
header_2 = {'x(nm)' 'y(nm)' 'z(nm)' 'cluster label (-1 is noise)'};
dbscan_data_combined_with_label=[header_2;num2cell([abs(-Y(:,1)) abs(-Y(:,2)) -Y(:,3) idx_dbscan])];

%Plotting dbscan result with noise and cluster data #Figure 2
fig2_name=sprintf('Density Based Clustering Data of Pick # %.f', pick_number);
fig2=figure('Name',fig2_name);
set(fig2,'visible','off');
scatter3(Y_new(:,1), Y_new(:,2), -Y_new(:,3), 15, 'filled');
hold on
scatter3(Y(:,1), Y(:,2), -Y(:,3), 15, idx_dbscan==-1, 'filled');
title(fig2_name);
xlabel('x');
ylabel('y');
zlabel('z');
pbaspect([1 1 1])
%saving
saveas(fig2,fullfile(folder_name,fig2_name),'fig');
saveas(fig2,fullfile(folder_name,fig2_name),'pdf');
dbscan_data_combined_with_label_name=sprintf('dbscan_data_combined_with_label_pick_num%.f.csv', pick_number);
writecell(dbscan_data_combined_with_label,fullfile(folder_name,dbscan_data_combined_with_label_name));

%clustering using k-means procedure start
%calculating inertia (elbow method)
  for i=1:cluster_trial
   opts = statset('Display','final');
   [idx,C,SSE] = kmeans(Y_new,i,'Maxiter',1000,'Replicates',10,'Options',opts);
   inertia_value(i)=sum(SSE);
   inertia_index(i)=i;
  end
%creating proper matrix data 
header_3 = {'num cluster' 'inertia'};
inertia_data=[header_3;num2cell([inertia_index' inertia_value'])];
%plotting inertia
fig3_name=sprintf('Inertia of Pick # %.f', pick_number);
fig3=figure('Name',fig3_name);
set(fig3,'visible','off');
plot(inertia_index,inertia_value);
title(fig3_name);
xlabel('Number of cluster');
ylabel('Inertia');
%saving
saveas(fig3,fullfile(folder_name,fig3_name),'fig');
saveas(fig3,fullfile(folder_name,fig3_name),'pdf');
inertia_data_name=sprintf('inertia_data_pick_num%.f.csv', pick_number);
writecell(inertia_data,fullfile(folder_name,inertia_data_name));


 %calculating gradient of the inertia
  for j=1:cluster_trial-1
       gradient_value(j)=inertia_value(j+1)-inertia_value(j);
       gradient_index(j)=j;
  end
%normalizing gradient
gradient_value_normalized=gradient_value/abs(min(gradient_value));

%creating proper matrix data 
header_4 = {'num cluster-1' 'normalized gradient'};
normalized_gradient_data=[header_4;num2cell([gradient_index' gradient_value_normalized'])];
%plotting gradient
fig4_name=sprintf('Normalized Gradient of Pick # %.f', pick_number);
fig4=figure('Name',fig4_name);
set(fig4,'visible','off');
plot(gradient_index,gradient_value_normalized)
title(fig4_name);
xlabel('Number of cluster-1');
ylabel('Normalized Gradient');
%saving
saveas(fig4,fullfile(folder_name,fig4_name),'fig');
saveas(fig4,fullfile(folder_name,fig4_name),'pdf');
normalized_gradient_data_name=sprintf('normalized_gradient_data_pick_num%.f.csv', pick_number);
writecell(normalized_gradient_data,fullfile(folder_name,normalized_gradient_data_name));

%picking optimum cluster
limit_gradient=-0.04 %<-------------------------------limit for gradient of elbow

 for j=1:cluster_trial-1
       if gradient_value_normalized(j)>limit_gradient
           list_cluster(j)=j;
       else
           list_cluster(j)=nan;
       end
       n_cluster=min(list_cluster)+1;
 end
%k-means clustering with optimum cluster
opts = statset('Display','final');
[idx,C,SSE] = kmeans(Y_new,n_cluster,'Maxiter',1000,'Replicates',50,'Options',opts);

%creating proper matrix data
header_5a = {'x(nm)' 'y(nm)' 'z(nm)' 'cluster label'};
kmeans_data_combined_with_label=[header_5a;num2cell([Y_new(:,1) Y_new(:,2) -Y_new(:,3) idx])];

header_5b = {'x(nm)' 'y(nm)' 'z(nm)'};
kmeans_centroid_data=[header_5b;num2cell([C(:,1) C(:,2) -C(:,3)])];

%plotting k-means result
fig5_name=sprintf('K-means Clustering of Pick # %.f', pick_number);
fig5=figure('Name',fig5_name);
set(fig5,'visible','off');
scatter3(Y_new(:,1), Y_new(:,2), -Y_new(:,3), 15, idx, 'filled');
hold on
scatter3(C(:,1), C(:,2), -C(:,3), 250,'kx');
title(fig5_name);
xlabel('x');
ylabel('y');
zlabel('z');
pbaspect([1 1 1])
%saving
saveas(fig5,fullfile(folder_name,fig5_name),'fig');
saveas(fig5,fullfile(folder_name,fig5_name),'pdf');
kmeans_data_combined_with_label_name=sprintf('kmeans_data_combined_with_label_pick_num%.f.csv', pick_number);
kmeans_centroid_data_name=sprintf('kmeans_centroid_data_pick_num%.f.csv', pick_number);
writecell(kmeans_data_combined_with_label,fullfile(folder_name,kmeans_data_combined_with_label_name));
writecell(kmeans_centroid_data,fullfile(folder_name,kmeans_centroid_data_name));


%Start alignment--------------------------------------

%finding center of mass
x_center_mass=mean(C(:,1));
y_center_mass=mean(C(:,2));
z_center_mass=mean(C(:,3));
%translation to zero center of mass
points=[C(:,1)-x_center_mass C(:,2)-y_center_mass -(C(:,3)-z_center_mass)];

%creating template cuboctahedron
factor=(max(points(:,3))-min(points(:,3)))/70 %z shrinking factor
Vcuboc=35*[1 1 0;1 0 -1*factor;0 1 -1*factor;-1 0 -1*factor;-1 1 0;-1 0 1*factor;-1 -1 0;0 1 1*factor;0 -1 1*factor;1 -1 0;1 0 1*factor;0 -1 -1*factor];
%VO=35*[1 1 0;0 1 -1*factor;-1 0 -1*factor;-1 1 0;-1 0 1*factor;-1 -1 0;0 1
%1*factor;1 -1 0;1 0 1*factor;0 -1 -1*factor];%initial pattern encryption
VO=35*[0 1 1*factor; 0 1 -1*factor; 0 -1 -1*factor; -1 -1 0; -1 1 0; 1 1 0; 1 -1 0;-1 0 1*factor;-1 0 -1*factor; 1 0 1*factor];%re-ordered to match the pattern encryption

%creating proper matrix data
header_6 = {'x(nm)' 'y(nm)' 'z(nm)'};
centroid_with_zero_cm_data=[header_6;num2cell(points)];

%plot new positions
fig6_name=sprintf('Clusters Centroid and Cuboctahedron Template of Pick # %.f', pick_number);
fig6=figure('Name',fig6_name);
set(fig6,'visible','off');
scatter3(points(:,1),points(:,2),points(:,3));
hold on;
pbaspect([1 1 1])
scatter3(VO(:,1), VO(:,2), VO(:,3), 250,'ko', 'filled', 'MarkerFaceColor','red');
scatter3(Vcuboc(:,1), Vcuboc(:,2), Vcuboc(:,3), 250,'ko', 'MarkerEdgeColor','red');
title(fig6_name);
xlabel('x');
ylabel('y');
zlabel('z');
pbaspect([1 1 1]);

%saving
saveas(fig6,fullfile(folder_name,fig6_name),'fig');
saveas(fig6,fullfile(folder_name,fig6_name),'pdf');
centroid_with_zero_cm_data_name=sprintf('centroid_with_zero_cm_data_pick_num%.f.csv', pick_number);
writecell(centroid_with_zero_cm_data,fullfile(folder_name,centroid_with_zero_cm_data_name));

%start translation alignment
%alignment input start------------------------------
n=4;%degree step
x_shift_initial=-20;%x translation initial
y_shift_initial=-20;%y translation initial
x_step=4;%x translation step
y_step=4;%y translation step
%alignment input end------------------------------
points_translate=points;
for i=1:1:-2*(x_shift_initial)/x_step 
    points_translate(:,1)=points(:,1)+x_shift_initial+x_step*(i-1);
        for j=1:1:-2*(y_shift_initial)/y_step 
            points_translate(:,2)=points(:,2)+y_shift_initial+y_step*(j-1);
            for k=1:1:360/n   %rotation alignment
                k
                R=rotz(n*k);
                points_rotate=points_translate*R';
                [index_weight,D] = knnsearch(VO,points_rotate);
                %-----add weight to orientation marker
                %------------------please continue working here, not
                %done
                
                for l=size(index_weight,1)
                    for m=1:1:7
                       if index_weight(l)==m
                        D(l)=2*D(l);%weight factor is introduced
                      end
                    end
                end
               
                squared_distance(k,j,i)=sum(D.^2);
            end
        end
end
                
%optimum alignment finalization
[squareddistance_min,index_min]=min(squared_distance(:));
[kk,jj,ii] = ind2sub(size(squared_distance),index_min);
points_temporary=[(points(:,1)+x_shift_initial+x_step*(ii-1)) (points(:,2)+y_shift_initial+y_step*(jj-1)) (points(:,3))];
points_new=points_temporary*(rotz(n*kk))';
%nearest neighbor projection after alignment
[Idx_final,D] = knnsearch(VO,points_new);
for j=1:1:n_cluster
    points_moved(j,:)=VO(Idx_final(j),:);
end

%creating proper matrix data and saving
header_7 = {'x(nm)' 'y(nm)' 'z(nm)'};
aligned_centroid_data=[header_7;num2cell(points_new)];
aligned_centroid_data_name=sprintf('aligned_centroid_data_pick_num%.f.csv', pick_number);
writecell(aligned_centroid_data,fullfile(folder_name,aligned_centroid_data_name));

header_8 = {'x(nm)' 'y(nm)' 'z(nm)'};
detected_points_attemplate_data=[header_8;num2cell(points_moved)];
detected_points_attemplate_data_name=sprintf('detected_points_attemplate_data_pick_num%.f.csv', pick_number);
writecell(detected_points_attemplate_data,fullfile(folder_name,detected_points_attemplate_data_name));

header_9 = {'Idx_final' 'Distance'};
Final_nearest_neighbor_data=[header_9;num2cell([Idx_final D])];
Final_nearest_neighbor_data_name=sprintf('Final_nearest_neighbor_data_pick_num%.f.csv', pick_number);
writecell(Final_nearest_neighbor_data,fullfile(folder_name,Final_nearest_neighbor_data_name));



%----------Decoder for cuboctahedron template start--------------
%VO index has been used as a template for decoder, the ordering of the
%positions has been numbered accordingly to match the encryption pattern
  number_value=0;
  position_value=0;
  filled_marker_present=0;
  unique_idx_final=unique(Idx_final);
for i=1:1:size(unique_idx_final,1)
    for j=1:1:4
        if unique_idx_final(i)==j
            number_value=number_value+2^(j-1);
        end
    end
    for k=5:1:7
        if unique_idx_final(i)==k
            position_value=position_value+2^(k-5);
        end
    end
    for l=8:1:10
        if unique_idx_final(i)==l
            filled_marker_present=filled_marker_present+1;
        end
    end
end
decoded(group_id,1)=number_value
decoded(group_id,2)=position_value
decoded(group_id,3)=filled_marker_present
if decoded(group_id,1)==4 && decoded(group_id,2)==2
    decoded(group_id,4)=1
else
    decoded(group_id,4)=0
end
header_10 = {'number' 'position' 'marker present' 'correct or not'};
decoded_data=[header_10;num2cell(decoded)];
decoded_data_name=sprintf('decoded_data_pick_num%.f.csv', pick_number);
writecell(decoded_data,fullfile(folder_name,decoded_data_name));

%plot aligned positions
fig7_name=sprintf('Final Alignment of Pick # %.f', pick_number);
fig7=figure('Name',fig7_name);
set(fig7,'visible','off');
scatter3(points_new(:,1),points_new(:,2),points_new(:,3),100,'ko','filled','MarkerFaceColor','magenta');
hold on;
scatter3(points_moved(:,1),points_moved(:,2),points_moved(:,3),1000, '*');
pbaspect([1 1 1])
hold on
k2 = convhull(double(Vcuboc(:,1)), double(Vcuboc(:,2)), double(Vcuboc(:,3)));
trisurf(k2,Vcuboc(:,1), Vcuboc(:,2), Vcuboc(:,3),'FaceColor',[0.4660 0.6740 0.1880], 'FaceAlpha',0.5)
hold on
scatter3(VO(:,1), VO(:,2), VO(:,3), 250,'ko', 'filled', 'MarkerFaceColor','blue');
scatter3(Vcuboc(:,1), Vcuboc(:,2), Vcuboc(:,3), 250,'ko', 'MarkerEdgeColor','blue');
title(fig7_name);
xlabel('x');
ylabel('y');
zlabel('z');
pbaspect([1 1 1]);
%saving
saveas(fig7,fullfile(folder_name,fig7_name),'fig');
saveas(fig7,fullfile(folder_name,fig7_name),'pdf');
header_11 = {'number' 'position' 'marker present' 'correct or not'};
decoded_data_global=[header_11;num2cell(decoded)];
decoded_data_global_name=sprintf('decoded_data_global.csv');
writecell(decoded_data_global,decoded_data_global_name);
end
header_11 = {'number' 'position' 'marker present' 'correct or not'};
decoded_data_global=[header_11;num2cell(decoded)];
decoded_data_global_name=sprintf('decoded_data_global.csv');
writecell(decoded_data_global,decoded_data_global_name);


load('b_Poly.mat');
figure,
for i = 1:size(b_Poly,1)
    plot(b_Poly{i}(:,1),b_Poly{i}(:,2),'-');hold on;
end
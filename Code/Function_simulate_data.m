function [output]=Function_simulate_data(n,xx,yy,A)       %  n refers to the number of DMU
% xx refers to the number of inputs, yy refers to the number of outputs
% A refers to the density

X=rand(n,xx)*1000;
n=size(X,1);
P=[X,zeros(n,yy)];
for k=1:yy
a=(1:xx);
random_num = a(randperm(numel(a),2));
ind1= sort(random_num);
X1=X(:,ind1(1));X2=X(:,ind1(2));
ind2=rand(1,1)*1;
ind3=1-ind2;
C=(X1.^ind2).*(X2.^ind3); 
P(:,xx+k)=C;
end
clearvars -except P and A and xx and yy and n
B=n*(1-A);
a=randperm(n);aa=a(1:B);
for k=1:B
random1=rand(1,yy);
random1=random1.*0.5
P(aa(k),xx+1:xx+yy)=P(aa(k),xx+1:xx+yy).*random1;
r2 = randi(2,1,xx);
r2 = r2;
P(aa(k),1:xx)=P(aa(k),1:xx).*r2;
end
[output]=P;

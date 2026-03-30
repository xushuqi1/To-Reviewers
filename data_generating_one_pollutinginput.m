function [output] = data_generating_one_pollutinginput (n,c,s,scale_parameter)   % 输入样本数n，c是产出间的关联系数，s是noise的分布方差,输出的最后一列是真实效率得分
X = 1 + 3 * rand(n, 4);          %   生成一个[1,4]的平均分布的投入矩阵
Y=zeros(n,1);
V=zeros(n,1);
score=zeros(n,1);                           
mu_1 = [0 0];                                %   无效程度系数模拟
Sigma_1 = [0.25, 0.25*c; 0.25*c, 0.25];      %   运行前需要check, 协方差矩阵
samples_1 = mvnrnd(mu_1, Sigma_1, n);        %   生成n个样本
samples_1 = abs(samples_1);                  %   正数
inefficiency_desirable=samples_1(:,1);
inefficiency_undesirable=samples_1(:,2);
mu_2 = [0 0];                                %   噪音系数模拟
Sigma_2 = [s, 0; 0, s];                      %   运行前需要check, 协方差矩阵
samples_2 = mvnrnd(mu_2, Sigma_2, n);        %   生成n个样本
noise_desriable=samples_2(:,1);
noise_undesriable=samples_2(:,2);
for k=1:n
    a=X(k,:);
    Y(k)=exp(scale_parameter*log(a(1))+scale_parameter*log(a(2))+scale_parameter*log(a(3))+scale_parameter*log(a(4))+ 0.5*[0.1*log(a(1))*log(a(1))-0.1*log(a(2))*log(a(1))-0.1*log(a(1))*log(a(2))+0.1*log(a(2))*log(a(2))+0.1*log(a(3))*log(a(3))-0.1*log(a(4))*log(a(3))-0.1*log(a(3))*log(a(4))+0.1*log(a(4))*log(a(4))]-inefficiency_desirable(k)+noise_desriable(k));
end                                          %   生成期望产出向量
for k=1:n
    a=X(k,:);
    V(k)=exp(4*scale_parameter*log(a(1))+ 0.5*0.1*log(a(1))*log(a(1))+inefficiency_undesirable(k)+noise_undesriable(k));
end                                          %   生成非期望产出向量
for k=1:n
score(k)=exp(-(inefficiency_desirable(k)+inefficiency_undesirable(k))/2);
end                                          %   生成效率得分
[output]=[X,Y,V,score];
end

    
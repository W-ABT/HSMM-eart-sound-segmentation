% BP神经网络
%%
%load data
clear
clc
load('PCGdata');

apha = 2;
beta = 2;
a = 1;
m = 0;
EWW1 = 0;
acc =0;

%%
% initialize first layer weights
% 设置S1，收缩期，S2，舒张期的权重
for h = 1:4
    for i = 1:3
        v(h,i) = 2*rand(1)-1;
    end
end

% initialize layers' bias
for j = 1:3
    b1(j) = 2*rand(1)-1;
    b2(j) = 2*rand(1)-1;
end

% initialize second layer weights
for i = 1:4
    for j = 1:3
        w(i,j) = 2*rand(1)-1;
    end
end

%%
%train
%comput hidden layer output and iteration

%while a == 1
%for a = 1:500000
    % 随机化处理
    B=randperm(length(PCGdata(1,:)));
    data=PCGdata(:,B(1:length(PCGdata(1,:)))); 
    % 分配训练集和测试集
    A = (1:515);
    C = (516:644);
    traindata=data(:,A(1:515)); 
    testdata = data(:,C(1:129)); 
    
    train = traindata{1,1};
    EW1 = 0;
    
    for n = 1:515
        for h = 1:4
            for i = 1:3
                x(h,i) = train(1,h) * v(h,i); 
            end
        end
        for i = 1:3
            fb(i) = x(1,i)+x(2,i)+x(3,i)+x(4,i)+b1(i);
            FB(i) = 1./(1+exp(-fb(i)));
        end
    %second layer
        for i = 1:3
            for j = 1:3
                z(i,j) = FB(i) * w(i,j);
            end
        end
        for i = 1:3
            fc(i) = z(1,i)+z(2,i)+z(3,i)+b2(i);
            FC(i) = 1./(1+exp(-fc(i)));
        end
    %comput the FC's error 
        if X{1,n}(1,5) == 0
            hot_y(1) = 0;
            hot_y(2) = 0;
            hot_y(3) = 1;
        end
    
        if X{1,n}(1,5) == 1
            hot_y(1) = 0;
            hot_y(2) = 1;
            hot_y(3) = 0;
        end
     
        if X{1,n}(1,5) == 2
            hot_y(1) = 1;
            hot_y(2) = 0;
            hot_y(3) = 0;
        end
        for i = 1:3
            y{1,n}(1,i) = hot_y(i);
        end
        for j = 1:3
            d(j) = FC(j)*(1-FC(j))*(hot_y(j)-FC(j));
        end
    %comput the FB's error
        for i = 1:3
            for j = 1:3
                er(i,j) = w(i,j)*d(j);
            end
        end
    
        for j = 1:3
            err(j) = er(j,1)+er(j,2)+er(j,3);
            e(j) = FB(j)*(1-FB(j))*err(j);
        end
    %update the second layer's weights
        for h = 1:3
            for j = 1:3
                w(h,j) =  apha*FB(h)*d(j)+w(h,j);
            end
        end
    %update the first layer's weights
        for i = 1:4
            for j = 1:3
                v(i,j) =  beta*X{1,n}(1,i)*e(j)+v(i,j);
            end
        end
    %updata the bias
        for j = 1:3
            b1(j) = beta*e(j)+b1(j);
            b2(j) = apha*d(j)+b2(j);
        end
        for i =1:3
            en(i) = y{1,n}(1,i)-FC(i);   
        end
        EW = ((en(1))^2+(en(2))^2+(en(3))^2)/2;
        EW1 = EW1 + EW;
    end
    m = m +1;
    EW1 = EW1/113;
%    if EW1 < 0.001
%        break
%    end
    EW1_fg(m) = EW1;
    
%end
t = 1:m-1;
plot(t,EW1_fg)
axis([1 m 0 1])
%%
%test    
for n = 114:150
    for h = 1:4
        for i = 1:3
            x(h,i) = X{1,n}(1,h) * v(h,i); 
        end
    end
    
    for i = 1:3
        fb(i) = x(1,i)+x(2,i)+x(3,i)+x(4,i)+b1(i);
        FB_Y(i) = 1./(1+exp(-fb(i)));
    end
    %second layer
    for i = 1:3
        for j = 1:3
            z(i,j) = FB_Y(i) * w(i,j);
        end
    end
    
    for i = 1:3
        fc(i) = z(1,i)+z(2,i)+z(3,i)+b2(i);
        FC_Y(i) = 1./(1+exp(-fc(i)));
    end
    %comput the loss
    if X{1,n}(1,5) == 0  
        hot_y(1) = 0;
        hot_y(2) = 0;
        hot_y(3) = 1;
    end
    
    if X{1,n}(1,5) == 1
        hot_y(1) = 0;
        hot_y(2) = 1;
        hot_y(3) = 0;
    end
     
    if X{1,n}(1,5) == 2
        hot_y(1) = 1;
        hot_y(2) = 0;
        hot_y(3) = 0;
    end
    for i = 1:3
        Y{1,n-113}(1,i) = hot_y(i);
    end
    Y_label{1,n-113} = FC_Y;
    

    %ew
    for i =1:3
        en(i) = Y{1,n-113}(1,i)-FC_Y(i);   
    end
    EWW(n-113) = ((en(1))^2+(en(2))^2+(en(3))^2)/2;
    EWW1 = EWW(n-113)+EWW1;
    
    
    %accuracy
    for i = 1:3
        if max(Y_label{1,n-113})-Y_label{1,n-113}(1,i) == 0
            Y_label{1,n-113}(1,i) = 1;
        else
            Y_label{1,n-113}(1,i) = 0;
        end
    end
    if Y_label{1,n-113} == Y{1,n-113}
        acc = acc+1;
    end
   
end
accuracy = acc/37;


    
    
    
    

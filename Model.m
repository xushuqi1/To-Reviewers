function L=model(P,xx,yy,zz)
n=size(P,1);
X=P(:,1:xx);Y=P(:,xx+1:xx+yy);V=P(:,xx+yy+1:xx+yy+zz);
nn1=size(X,2);nn2=size(Y,2);nn3=size(V,2);
for k=1:n   
     c=[zeros(3*n,1);-1];
     A=[zeros(nn2,n),-Y',zeros(nn2,n),Y(k,:)';
        zeros(nn3,2*n),V',V(k,:)';
        -1.*eye(n),eye(n),zeros(n,n+1);
        -1.*eye(n),zeros(n,n),eye(n),zeros(n,1);
        ];
     b=[-Y(k,:)';
         V(k,:)';
         zeros(n,1);
         zeros(n,1);
         ];
     Aeq=[X',zeros(nn1,2*n+1);];
     beq=[X(k,:)';];
     lb=[zeros(3*n+1,1)];
     ub=[inf*ones(3*n,1);0.999999999999];
     W(:,k)=cplexlp(c,A,b,Aeq,beq,lb,ub);
end
L1=W(3*n+1,:);
L=1-L1;
L=L';
end

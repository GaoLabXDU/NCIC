% X is the expression matrix , rank is the  rank of gene 
function [label]=NMF_1_type(X,Y,rank_X,rank_Y,K)
[ng,ns]=size(X);
[ng2,ns]=size(Y);
W=sparse(diag(rank_X));
W2=sparse(diag(rank_Y));
T=1000;
eps=10e-6;
res=[0];
ncs=K;
ncg=K;
% initialize the clusters using k-means
G=sparse(1:ng,kmeans(X,ncg)',ones(1,ng),ng,ncg);
F=sparse(1:ns,kmeans(X',ncs)',ones(1,ns),ns,ncs);

G2=sparse(1:ng2,kmeans(Y,ncg)',ones(1,ng2),ng2,ncg);
F2=sparse(1:ns,kmeans(Y',ncs)',ones(1,ns),ns,ncs);

for i=1:T
        % data type1
		S=(G'*W*G)\G'*W*X*F/(F'*F);
        A=X'*W*G*S;
        B=S'*G'*W*G*S;
        Aa=abs(A);
        Ba=abs(B);
        F=F.*sqrt((Aa+A+F*(Ba-B))./(Aa-A+F*(Ba+B)+eps));
        
        C=W*X*F*S';
        D=S*F'*F*S';
        Ca=abs(C);
        Da=abs(D);
        G=G.*sqrt((Ca+C+W*G*(Da-D))./(Ca-C+W*G*(Da+D)+eps));
        
        s=sum(F,2);
        iz=find(s<=0);
        ip=find(s>0);
    	F(ip,:)=F(ip,:)./repmat(s(ip),1,ncs);
    	F(iz,:)=1/ncs;
        
        s=sum(G,2);
        iz=find(s<=0);
        ip=find(s>0);
        G(ip,:)=G(ip,:)./repmat(s(ip),1,ncg);
        G(iz,:)=1/ncg;
		% data type 2
		S2=(G2'*W2*G2)\G2'*W2*Y*F2/(F2'*F2);
        A2=Y'*W2*G2*S2;
        B=S2'*G2'*W2*G2*S2;
        A2a2=abs(A2);
        Ba2=abs(B);
        F2=F2.*sqrt((A2a2+A2+F2*(Ba2-B))./(A2a2-A2+F2*(Ba2+B)+eps));
        
        C2=W2*Y*F2*S2';
        D2=S2*F2'*F2*S2';
        C2a2=abs(C2);
        D2a2=abs(D2);
        G2=G2.*sqrt((C2a2+C2+W2*G2*(D2a2-D2))./(C2a2-C2+W2*G2*(D2a2+D2)+eps));
        
        s2=sum(F2,2);
        iz2=find(s2<=0);
        ip2=find(s2>0);
    	F2(ip2,:)=F2(ip2,:)./repmat(s2(ip2),1,ncs);
    	F2(iz2,:)=1/ncs;
        
        s2=sum(G2,2);
        iz2=find(s2<=0);
        ip2=find(s2>0);
        G2(ip2,:)=G2(ip2,:)./repmat(s2(ip2),1,ncg);
        G2(iz2,:)=1/ncg;
        
		res1=trace(X'*W*X-2*X'*W*G*S*F'+F*S'*G'*W'*G*S*F')+trace(Y'*W2*Y-2*Y'*W2*G2*S2*F2'+F2*S2'*G2'*W2'*G2*S2*F2');
		res=[res res1];
    	if abs(res1-res(i))<eps
            disp(['Converge at',num2str(i)]);
			break;
        end
end
group=zeros(1,ns);
for i=1:ns
	[m index]=max(F(i,:));
	group(i)=index;
end
label=group';
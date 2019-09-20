% X is the expression matrix , rank is the  rank of gene 
function [label]=NMF_1_type(X,rank,K)
[ng,ns]=size(X);
W=sparse(diag(rank));
T=1000;
eps=10e-6;
res=[0];
ncs=K;
ncg=K;
% initialize the clusters using k-means
G=sparse(1:ng,kmeans(X,ncg)',ones(1,ng),ng,ncg);
F=sparse(1:ns,kmeans(X',ncs)',ones(1,ns),ns,ncs);

for i=1:T
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
        
        res1=trace(X'*W*X-2*X'*W*G*S*F'+F*S'*G'*W'*G*S*F');
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
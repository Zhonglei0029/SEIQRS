function [t,S,E,I,Q,R] = Program_SEIQRS(N,n,tau1,tau2,gamma,delta,MaxTime,Type)

% Sets up default parameters if necessary.
if nargin == 0
    N=100; 
    n=4;
    tau1=0.1; % 感染率
    tau2=0.05; % 潜伏期传播率
    gamma=0.1; % 康复率
    delta=0.1; % 隔离率
    MaxTime=1000;
    Type='Random';
end

% Checks all the parameters are valid
CheckGreater(N,0,'Number of individuals N');
CheckGreater(n,0,'Number of neighbours n');
CheckGreater(tau1,0,'tau1');
CheckGreater(tau2,0,'tau2');
CheckGreater(gamma,0,'gamma');
CheckGreater(delta,0,'delta');
CheckGreater(MaxTime,0,'MaxTime');

% 初始化网络
[X,Y,G,N]=Create_Network(N,n,Type);
Status=1+0*X; Status(1)=2; % 初始状态，第一个个体为感染状态
Rate=0*X; Rate(1)=gamma; Rate(find(G(:,1)))=tau1;

t=0; i=1; S=N-1; E=0; I=1; Q=0; R=0; % 初始化时间，迭代计数器和状态数量
subplot(2,1,1);
[j,k,s]=find(G);
plot([X(j) X(k)]',[Y(j) Y(k)]','-k');
hold on
Col=[0.7 0.7 0.7; 0 1 0; 1 0 0; 0 0 1; 1 1 0; 0 1 1];
for k=1:N
    H(k)=plot(X(k),Y(k),'ok','MarkerFaceColor','g');
end
set(H(1),'MarkerFaceColor','r');
hold off;
drawnow;

while (t<MaxTime & I(end)>0)
    [step,Rate,Status,e]=Iterate(Rate,Status,G,tau1,tau2,gamma,delta);
    i=i+1;
    t(i)=t(i-1)+step;
    S(i)=length(find(Status==0));
    E(i)=length(find(Status==1));
    I(i)=length(find(Status==2));
    Q(i)=length(find(Status==3));
    R(i)=length(find(Status==4));
    set(H(e),'MarkerFaceColor',Col(Status(e)+1,:));

    subplot(4,1,3);
    plot(t,S,'-g');
    ylabel 'Number of Susceptibles'
    subplot(4,1,4);
    plot(t,I,'-r');
    ylabel 'Number of Infecteds'
    xlabel 'Time'
    drawnow;
    [];
end


function [X,Y,G,N]=Create_Network(N,n,Type)

if n > (N-1)
    error('Impossible to have an average of %d contacts with a population size of %d',n,N);
end

G=sparse(1,1,0,N,N);
X=rand(N,1); Y=rand(N,1);
switch Type

    case {'Random','random'}
        contacts=0;
        while(contacts<n*N)
            i=ceil(rand(1,1)*N); j=ceil(rand(1,1)*N);
            if i~=j && G(i,j)==0
                contacts=contacts+2;
                G(i,j)=1; G(j,1)=1;
            end
        end
                
    case {'Spatial','spatial'}
        D=(X*ones(1,N) - ones(N,1)*X').^2 + (Y*ones(1,N) - ones(N,1)*Y').^2;
        Prob=exp(-D*5)-rand(N,N); Prob=triu(Prob,1)-1e100*tril(ones(N,N),0);
        [y i]=sort(reshape(Prob,N*N,1));
        p=i(end+[(1-n*N/2):1:0]);
        i=1+mod(p-1,N); j=1+floor((p-1)/N);
        G=sparse([i; j],[j; i],1,N,N);

    case {'Lattice','lattice'}
        if mod(sqrt(N),1)~=0
            warning('N=%d is not a square number, rounding to %d',N,round(sqrt(N)).^2);
            N=round(sqrt(N)).^2;
        end
        [X,Y]=meshgrid([0:(sqrt(N)-1)]/(sqrt(N)-1));
        X=reshape(X,N,1); Y=reshape(Y,N,1);
        D=(X*ones(1,N) - ones(N,1)*X').^2 + (Y*ones(1,N) - ones(N,1)*Y').^2;
        Prob=triu(D,1)+1e100*tril(ones(N,N),0);
        [y i]=sort(reshape(Prob,N*N,1));
        p=i(1:(n*N/2 - 2*sqrt(N)));
        i=1+mod(p-1,N); j=1+floor((p-1)/N);
        G=sparse([i; j],[j; i],1,N,N);

    case {'SmallWorld','Smallworld','smallworld'}
        if mod(sqrt(N),1)~=0
            warning('N=%d is not a square number, rounding to %d',N,round(sqrt(N)).^2);
            N=round(sqrt(N)).^2;
        end
        [X,Y]=meshgrid([0:(sqrt(N)-1)]/(sqrt(N)-1));
        X=reshape(X,N,1); Y=reshape(Y,N,1);
        D=(X*ones(1,N) - ones(N,1)*X').^2 + (Y*ones(1,N) - ones(N,1)*Y').^2;
        Prob=triu(D,1)+1e100*tril(ones(N,N),0);
        [y i]=sort(reshape(Prob,N*N,1));
        p=i(1:(n*N/2 - 2*sqrt(N)));
        i=1+mod(p-1,N); j=1+floor((p-1)/N);
        G=sparse([i; j],[j; i],1,N,N);
        R=10;
        for k=1:R
            i=1; j=1;
            while (i==j || G(i,j)==1)
                i=ceil(rand(1,1)*N); j=ceil(rand(1,1)*N);
            end
            G(i,j)=1; G(j,1)=1;
        end
        
    otherwise
        error('Do not recognise network type %s',Type);
end



function [step,Rate,Status,Event]=Iterate(Rate,Status,G,tau1,tau2,gamma,delta)

Sum=sum(Rate); Cum=cumsum(Rate);
Event=min(find(Cum>rand(1,1)*Sum));
Status(Event)=mod(Status(Event)+1,5);

contacts=find(G(:,Event) & Status==1);
switch Status(Event)
    case 1 % 感染者变为潜伏者
        Rate(Event)=tau2;
        Rate(contacts)=Rate(contacts)+tau1;
    case 2 % 潜伏者变为隔离者
        Rate(Event)=delta;
        Rate(contacts)=Rate(contacts)+tau2;
    case 3 % 隔离者变为康复者
        Rate(Event)=gamma;
        Rate(contacts)=Rate(contacts)+delta;
        G(Event,:)=0;
    case 0 % 康复者不再传播疾病
        Rate(Event)=0;
        Rate(contacts)=Rate(contacts)-gamma;
end
Rate=Rate.*sign(Status);
step=-log(rand(1,1))/Sum;

function []=CheckGreaterOrEqual(Parameter, Value, str)

m=find(Parameter<Value);
if length(m)>0
    error('Parameter %s(%g) (=%g) is less than %g',str,m(1),Parameter(m(1)),Value);
end

function []=CheckGreater(Parameter, Value, str)

m=find(Parameter<=Value);
if length(m)>0
    error('Parameter %s(%g) (=%g) is less than %g',str,m(1),Parameter(m(1)),Value);
end

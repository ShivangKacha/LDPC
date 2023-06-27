    % clear all
    % clc;
    % for 9*12 hmatrix
    hs=[1,0,0,0,0,1,0,1,0,1,0,0;
        1,0,0,1,1,0,0,0,0,0,1,0;
        0,1,0,0,1,0,1,0,1,0,0,0;
        0,0,1,0,0,1,0,0,0,0,1,1;
        0,0,1,0,0,0,1,1,0,0,0,1;
        0,1,0,0,1,0,0,0,1,0,1,0;
        1,0,0,1,0,0,1,0,0,1,0,0;
        0,1,0,0,0,1,0,1,0,1,0,0;
        0,0,1,1,0,0,0,0,1,0,0,1];
    
    % data=load('Hmatrix.mat'); % for hmatrix1
    % data=load('Hmatrix2.mat'); % for hmatrix2
    % hs=data.H;
    
    %cnvn = check nodes connected with variable nodes
    for i=1:size(hs,1)
        row=hs(i,:);
        indi=find(row);
        cnvn(i,:)=indi;
    end
    
    %vncn = variable nodes connected with check nodes
    for i=1:size(hs,2)
        col=hs(:,i);
        indi=find(col);
        vncn(:,i)=indi;
    end
    len=length(vncn(:,1));
    n=size(hs,2);
    u=size(hs,1);
    
    % % for decoding success probability
    step=0.02; %for 9*12 matrix
    % step=0.05; %for hmatrix 1 and 2
    answer=zeros(1,1/step);
    probalbility=0:step:1;
    temp=1;
    
    for p = 0:step:1
            
        % index=20; % for hmatrix1 and 9*12 matrix
        index=100; % for hmatrix2
        err=0;
        itr=1:index;
        erap=zeros(1,index);
        % nsiim=50;
        nsiim=1000; %for 9*12 hmatrix
        for nsim=1:nsiim
            bcap=zeros(u,n); 
            message = zeros(n,1);
            noise=rand(n,1)<p;
            send = message;
            send(noise==1) = n+1;
            b=repmat(send,1,u);
            
            
            for it=1:index
                era=0;

                %for cn to vn
                for i=1:u
                    le=length(cnvn(i,:));
                    sumi=0;
                    for k=1:le
                        sumi=sumi+b(cnvn(i,k),i);
                    end
                    
                    for j=1:le
                        su=sumi-b(cnvn(i,j),i);
                        if su>le
                            bcap(i,cnvn(i,j))=n+1;
                        else
                            bcap(i,cnvn(i,j))=mod(su,2);
                        end
                    end
                end
    
                %for vn to cn
                for j=1:n
                    le=length(vncn(:,j));
                    sumi=0;
                    for k=1:le
                        sumi=sumi+bcap(vncn(k,j),j);
                    end
                    sumi=sumi+send(j);
                    for k=1:le
                        su=sumi-bcap(vncn(k,j),j);
                        if(su==le*(n+1))
                            b(j,vncn(k,j))=n+1;
                        else
                            su=mod(su,n+1);
                            if(su==0)
                                b(j,vncn(k,j))=0;
                            else
                                b(j,vncn(k,j))=1;
                            end
                        end
                    end
                end
            end
    
            recieved = zeros(1,n);
    
            for i = 1:n
                l = 1 + len;
                recieved(i) = (sum(bcap(vncn(:,i),i)) + send(i));
                if(recieved(i)==l*(n+1))
                    recieved(i)=n+1;
                else
                    recieved(i)=mod(recieved(i),n+1)>0;
                end
            end
            recieved=transpose(recieved);
    
            if(sum(xor(recieved,message))~=0)
                err=err+1;
            end
    
        end
        
        prob = 1 - err/nsim;
        answer(temp)=prob;
        temp=temp+1;
        p % for knowing the progress of the algorithm
    end
    
    plot(probalbility,answer);
    legend('success');
    xlabel('Error probability');
    ylabel('Decode Success probability');
    
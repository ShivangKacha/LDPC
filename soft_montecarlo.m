    % clear all;
    % clc;
    % 
    hs=[1,0,0,0,0,1,0,1,0,1,0,0;
        1,0,0,1,1,0,0,0,0,0,1,0;
        0,1,0,0,1,0,1,0,1,0,0,0;
        0,0,1,0,0,1,0,0,0,0,1,1;
        0,0,1,0,0,0,1,1,0,0,0,1;
        0,1,0,0,1,0,0,0,1,0,1,0;
        1,0,0,1,0,0,1,0,0,1,0,0;
        0,1,0,0,0,1,0,1,0,1,0,0;
        0,0,1,1,0,0,0,0,1,0,0,1];
    % data=load('Hmatrix.mat');
    % data=load('Hmatrix2.mat');
    % hs=data.H;

    for i=1:size(hs,1)
        row=hs(i,:);
        indi=find(row);
        cnvn(i,:)=indi;
    end

    for i=1:size(hs,2)
        col=hs(:,i);
        indi=find(col);
        vncn(:,i)=indi;
    end

    len=length(vncn(:,1));
    n=size(hs,2);
    u=size(hs,1);
    % index=100; % for hmatrix2
    index=20; % for 9*12 hmatrix and hmatrix1

    % for monte carlo
    step=0.02; %for 9*12 matrix
    % step=0.05; %for hmatrix 1 and 2
    answer=zeros(1,1/step);
    probalbility=0:step:1;
    temp=1;

    for p = 0:step:1 % for monte carlo
        err=0;
        itr=1:index;
        erap=zeros(1,index);
        % nsiim=50; % for hmatrix1 and hmatrix2
        nsiim=1000; % for 9*12 hmatrix
        for nsim=1:nsiim
            p,nsiim
            bcap=zeros(u,n); 
            message = zeros(n,1);
            noise=rand(n,1)<p;
            send = message;
            send(noise==0) = 0.01;
            send(noise==1) = 0.5;
            b=repmat(send,1,u);


            for it=1:index
                era=0;

                % CN to VN
                for i=1:u
                    le=length(cnvn(i,:));
                    for j=1:le
                        sumi=0;
                        for k=1:le
                            if k==j
                                continue;
                            else
                                sumi=sumi*(1-b(cnvn(i,k),i))+b(cnvn(i,k),i)*(1-sumi);
                            end
                        end

                        bcap(i,cnvn(i,j))=sumi;
                    end
                end

                % VN to CN
                for j=1:n
                    le=length(vncn(:,j));
                    for i=1:le

                        p1=send(j);
                        p0=1-send(j);
                        for k=1:le
                            if i==k
                                continue;
                            else
                                p1=p1*bcap(vncn(k,j),j);
                                p0=p0*(1-bcap(vncn(k,j),j));
                            end
                        end
                        b(j,vncn(i,j))=p1/(p0+p1);
                    end
                end
            end

            recieved = zeros(1,n);

            for j=1:n
                le=length(vncn(:,j));
                p1=send(j);
                p0=1-send(j);
                for k=1:le
                    p1=p1*bcap(vncn(k,j),j);
                    p0=p0*(1-bcap(vncn(k,j),j));
                end
                if (p1/(p0+p1))==0.5
                    recieved(j)=0.5;
                elseif (p1/(p0+p1))>0.5
                    recieved(j)=1;
                else
                    recieved(j)=0;
                end
            end

            recieved=transpose(recieved);  
            if(sum(xor(recieved,message))~=0)
                err=err+1;
            end
        end
        
        prob = 1 - err/nsiim;
        answer(temp)=prob;
        temp=temp+1;
    end

    plot(probalbility,answer);
    legend('success');
    xlabel('Error probability');
    ylabel('Decode Success probability');
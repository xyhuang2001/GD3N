 clear;close all;clc;
whzzongjie_mean=[];
 whzzongjie_var=[];
for whzhigh=230   
    for whzlow=50
zongjie=[];
jilu=[];
deathtime=[];
%% ��������
   W=inline('0.5*(x/y)+0.3*(1/z)+0.2*(p/0.52)','x','y','z','p');  
   % x����ýڵ�δ�ɴ��ھ�����y����x���ֵ��z����ýڵ���ȣ�p����ýڵ�ʣ������  
   %%Etx����һ���ڵ㴫���ȥ���������ĵ�����
   
   Eelec=50e-9; %Eelec
   Efs=10e-12;%ϵ��1
   Emp=0.0013e-13;%ϵ��2
   Etx1=inline('64*50e-9+64*10e-12*d*d','d');   %Etx1 ����64bit���ݰ���Ϊ�ɴ���Ϣ��
   Etx2=inline('552*50e-9+552*10e-12*d*d','d');   %Etx2 ����552bit���ݰ���Ϊ������Ϣ��
   Etx3=inline('80*50e-9+80*10e-12*d*d','d');   %Etx2 ����80bit���ݰ���ΪIN��Ϣ��
%%Erx����һ���ڵ�Ҫ�ǽ������ݣ����ĵ�����
%%Erx����һ���ڵ�Ҫ�ǽ������ݣ����ĵ�����

    Erx=inline('k*50e-9','k'); %һԪ������k�������ݰ�������
    
%%һ��Ƚڵ�ÿһ�����������������壬�ú�������ϵͳ����

   Eperround=inline('(500-x)*k*50e-9 + 500*y','x','k','y'); 
   %x����һ��Ƚڵ�����k�������ݰ���������y����Etx(k,d)
   simulation=0;
   
while(simulation~=10)

%% ����ĳ�ʼ���������ڵ�Ĳ�������Ϣ��ʼ��
close all; clc
deathtime
simulation
my=rand(1,70000000);
death=0;
xm=400;
ym=400;
r=200; theta=0:pi/200:2*pi;
x=200+r*cos(theta); y=200+r*sin(theta);
time=0;
flag=0;
confirm=[]; 
suspect_now=[];
plot(x,y,'-')
hold on; axis equal
% ��500���ڵ㣬ͨ�Ű뾶75m
FDRM=[];
MDRM=[];
nn=0;
R=75;
E0=1.7;  %��ʼ����2J
numCH=0;   %��ͷ����
child=[];    %child�������ڱ��������ͷ�ĺ��ӱ��
energy=[];
n=500;
%SN_receive_packetnum_theory=0;%SN�հ�������
%SN_receive_packetnum_reality=0;%SN�հ�ʵ����
%����������
%throughput=0;
throughput_rate=0;
%send=0;
%node_id=[];
%tun=0;
%tu=0;

i=0;
rt=501;
x_BS=200;
y_BS=200;
plot( x_BS, y_BS, 'd','markerfacecolor',[1,0,0]); %BS����
    extra_for=[];
    judge_now=[];
    ever_judge=[];
    release_judge=[];
hold on

while(nn~=n)

        rt=rt+1;
        xd=my(1,rt)*xm;
        rt=rt+1;
        yd=my(1,rt)*ym; 
        if (sqrt((xd-200)^2+(yd-200)^2 )<r)
            i=i+1;
            %node_id=[node_id i];
            S(i).x=xd;
            S(i).y=yd;
            rt=rt+1;
            infor(i).frate = 0.8+0.2*my(1,rt); %��ÿ���ڵ��ʼ��һ������ת����
            infor(i).malicious = 0;        
            infor(i).RE = E0;    %infor(i).RE��ʾi�ڵ��ʣ����������ʼ����Ϊ2J
            infor(i).pro = 'N';    %�����Ƿ���⣬ȫ�̲���ı䣬NΪ������MΪ����
            infor(i).successive = 0;
            infor(i).type = 'MN';
            infor(i).IN = 0;
            infor(i).harsh = 0;
            infor(i).fratetemp = 0;
            infor(i).locationx = 0;
            infor(i).locationy = 0;
            infor(i).increase = 0;
            infor(i).conssuspect = 0;
            infor(i).receive=zeros(1,5000);
            infor(i).forward=zeros(1,5000);
            infor(i).q_return=zeros(10,10);
            plot(S(i).x,S(i).y,'wo');
            hold on;
            nn=nn+1;  
        end
        
end
%% ���ֽ�������ʼģ��
    k=5000;
    Receive=zeros(n,k+1);    %��ʼ����ͷ�հ�����
    Receive(:,2)=2;
    Forward=zeros(n,k+1);    %��ʼ����ͷ����������
    Forward(:,2)=2;
    CFR=zeros(n,k+1);        %��ʼ���ۼ�ת���ʾ���
    CFRtemp=zeros(n,k+1);        %��ʼ���ۼ�ת���ʾ���
    
while(death==0)     %time�������
    time=time+1
    %uu=plotnetwork(S,time,infor);
    %% ���Ȼ�ԭ
  for w=2:200
    if(time==w)
        for orig=1:n
            if(infor(orig).harsh==time-1)
                infor(orig).frate = infor(orig).fratetemp;%��ԭ
            end
        end
    end
  end
    %% ģ�ⲿ�������ŵ�������
  for w=1:200
    if(time==w)
       harshnode=0;
       harsh_node=[];
       xharsh=0;
       yharsh=0;
       while(sqrt((xharsh-200)^2+(yharsh-200)^2 )>=200)
        rt=rt+1;
        xharsh=my(1,rt)*xm;
        rt=rt+1;
        yharsh=my(1,rt)*ym; 
       end
        
       for har=1:n
          if(sqrt((S(har).x-xharsh)^2+(S(har).y-yharsh)^2)<37.5)  %harsh����
                 plot(S(har).x,S(har).y,'ro','MarkerFaceColor','r');
                 infor(har).fratetemp=infor(har).frate;  %��������ת����
                 infor(har).frate = infor(har).frate*0.8;  %�ı�ת����
                 %infor(har).frate = 0.15;  %�ı�ת����
                 infor(har).harsh=time;
                 harshnode=harshnode+1;
                 harsh_node = [harsh_node har];
          end
       end   
    end
  end
    
    
    child=[];
    ORI=[];
    numCH=0;

    numIN=0;

    numMN=0;
   %% �ɴ�
           c=zeros(n,2);
           for i=1:n  %��ʼ����һ���� infor����ڵ�����
              infor(i).clu = 0;  %clu���Դ���ýڵ�δ�ɴأ��ɴ���Ϊ1
              infor(i).dep = 0;  %dep��ʾ�ýڵ���ȣ���ʼʱ������Ϊ0
              c(i,:)=[S(i).x,S(i).y];  
              Receive(i,1)=i;
              Forward(i,1)=i;
              CFR(i,1)=i;
           end
           %���Ŀǰ���ߴ�����������ڽӾ���neighbor_matrix
           for i=1:n
              R_dist = sqrt((sum(transpose((repmat(c(i,:),n-1,1)-c([1:i-1,i+1:n],:)).^2))));
              neighbor_i = R_dist< R/2 & infor(i).clu == 0;  %ͨ�ŷ�Χ��+δ�ɴأ�����Ϊ1���������е�δ�ɴ��ھ�
              neighbor_matrix(i,:) = [neighbor_i(1:i-1),0,neighbor_i(i:n-1)];
           end
           %����������нڵ���ھ���
           neighbor_num=zeros(1,n);
           for q=1:n
               neighbor_num(1,q)=sum(neighbor_matrix(q,:));  %�洢ÿ������ھ���
           end
           %Nmax=max(max(neighbor_num));  %Nmax��NB���ֵ���ڴ˻�ȡ
           Nmax=14;
           %fprintf('Nmax%d\n',Nmax)
           %����ͳ��ÿ���ڵ���ھӱ��
           nighbor_no=zeros(n,40);
           for d=1:n
               tempn=neighbor_matrix(d,:);   %��ȡ�ض���
               neighjilu=find(tempn~=0);
               for j=1:size(neighjilu,2)
                   nighbor_no(d,j)=neighjilu(j);
               end
           end
           fid12=fopen(['20230603','.txt'],'w');
           for whz=1:n
               for whzz=1:size(nighbor_no(whz,:),2)
                   whzzz=nighbor_no(whz,whzz);
                   fprintf(fid12,'%d',whzzz);
                   fprintf(fid12,'%s',' ');
               end
               fprintf(fid12,'\n');
           end
           fclose(fid12);
           
           for i=1:n %��ʼ���ڶ�����
               infor(i).NB = sum(neighbor_matrix(i,:)~=0,2);   % infor(i).NB����i�ڵ��δ�ɴ��ھ�������ʼ���������ھ���
               infor(i).vis = 0;   %̽��״̬Ϊδ̽��
               infor(i).par = 0;  %��ʼʱ�����нڵ�ĸ�������Ϊ0
               infor(i).MNnum = 0;
               infor(i).W = 0;
               if(infor(i).type ~= 'DE')
                 infor(i).type = 'MN';  %��ʼʱ���ڵ����;�����ΪMN
               end
               %infor(i).W = 0.5*(neighbor_num(1,i)/Nmax)+0.2;  %��ʼW
           end
      
      
   MinDist_Node_BS = zeros(1,n);  %������һ������ Ϊ�˺����������
   infor(n+1).dep = 0;
   infor(n+1).MNnum = 0;
   S(n+1).x=200;
   S(n+1).y=200;
   

 %% STEP1 STEP2
 BS_NB=0;  %��¼BS�ĺ�������
   for i=1:n
      MinDist_Node_BS(i) = sqrt( (S(i).x-200)^2 + ( S(i).y-200 )^2  );  %�����������нڵ㵽BS�ľ���
      % ����ýڵ���BS����С��R�������һ����,�����в���
      if(MinDist_Node_BS(i)<R/2 & infor(i).malicious==0)
           plot([S(i).x,200],[S(i).y,200],'bo-');  %�ڵ����ͷ����
           infor(i).vis = 1;   %̽��״̬Ϊ��̽��
           infor(i).clu = 1;  %���ڳɴأ�clu��������Ϊ1
           %cluinf(1,i) = 1;   %ͬʱ���¾���
           infor(i).dep = 1;  %��һ���أ��ڵ����Ϊ1
           infor(i).par = n+1;  %������BS����501��ʾ
           infor(i).dist = MinDist_Node_BS(i);  %�����³�Ա����ͷ�ľ���
           BS_NB=BS_NB+1; %�³�Ա��+1
           child(n+1,BS_NB+1)=i;      %�����վ�ĺ���
           infor(i).W = W(infor(i).NB,Nmax,infor(i).dep,infor(i).RE);
           infor(i).RE = infor(i).RE - Etx1(infor(i).dist);%����һ�����ݣ�����һ������
           [row,col] = find(neighbor_matrix(i,:)~=0 ); % �ҵ��ýڵ���ھӣ�Ŀ���ǰ������ھӵ�NB����-1����Ϊ�ýڵ�ɴ���
           if(size(col,2)~=0)  %����ýڵ����ھӣ��������ھӵ�NB����-1
             len=size(col,2);  %�ھ���
             for j=1:len
                 infor(col(j)).NB=infor(col(j)).NB-1;
             end
           end
           % Ȼ��BS�����һ������ÿ���ڵ��Wֵ
           infor(i).dep=1;  
      end    
   end
   
   %figure(time);% ��ʾͼƬ 
  %%STEP2
    %��ȡWֵ���Ľڵ���Ϊ��һ����ͷ
    maxW=0;
    newCH=0;
    
   for i=1:n
      if((infor(i).W>maxW) & (infor(i).vis~=0) & (infor(i).type~='CH') & infor(i).NB>1)
           maxW=infor(i).W;
           newCH=i;
      end
   end
   
   if(newCH==0)
       fprintf('δ�ҵ���һ����ͷ��');
   end
   
   numCH=numCH+1;    %��ͷ����һ
   infor(numCH).sort = newCH;    %��ͷ�±��
   child(numCH,1)=newCH;   %��һ���Ǵ�ͷ���
    %�õ���һ����ͷ֮�󣬸����´�ͷ����
    zht1=plot(S(newCH).x,S(newCH).y,'ko','MarkerFaceColor','k');     % ���Ƶ�ǰ��ͷ�ڵ�
    
    infor(newCH).vis = 1;   %̽��״̬Ϊ��̽��
    infor(newCH).type = 'CH';  %���ڸýڵ��Ǵ�ͷ�ڵ㣬��������ΪC
    infor(newCH).clu = 1;  %���ڳɴأ�clu��������Ϊ1
    cluinf(1,newCH) = 1;   %ͬʱ���¾���
    infor(newCH).dep = 1;  %��ͷ���ڵ����Ϊ1
    [row,col] = find(neighbor_matrix(newCH,:)~=0 ); % �ҵ��ýڵ���ھӣ�Ŀ���ǰ������ھӵ�NB����-1����Ϊ�ýڵ�ɴ���
    if(size(col,2)~=0)  %����ýڵ����ھӣ��������ھӵ�NB����-1������δ�ɴ��ھӼ���ô�
          len=size(col,2);  %�ھ���
          count_NB=0;  % count_NB�����³�Ա������ʼ��Ϊ0
          MN2_enermax = 0; 
          IN2 = 0;
          for i=1:len
              if(infor(col(i)).clu == 0 & infor(col(i)).malicious==0) %����ھ�δ�ɴأ������ôأ�����ôص�ͬʱҪ�������ھӽڵ��NB 
                  infor(col(i)).vis = 1;   %̽��״̬Ϊ��̽��
                  plot([S(col(i)).x,S(newCH).x],[S(col(i)).y,S(newCH).y],'bo-');  %�ڵ����ͷ����
                  infor(col(i)).par=newCH;  %���ø��ף��ôش�ͷ��ΪnewCH
                  infor(col(i)).dep=infor(newCH).dep+1;  %�ڵ���ȵ��ڸ������+1
                  infor(col(i)).clu=1; %�ɴأ���־Ϊ1
                  cluinf(1,col(i)) = 1;   %ͬʱ���¾���
                  infor(col(i)).dist = sqrt( (S(col(i)).x - S(newCH).x)^2 + ( S(col(i)).y - S(newCH).y )^2);  %�����³�Ա����ͷ�ľ���
                  infor(col(i)).W=W(infor(col(i)).NB,Nmax,infor(col(i)).dep,infor(col(i)).RE);%�����³�Ա��Wֵ
                  infor(col(i)).RE = infor(col(i)).RE - Etx1(infor(col(i)).dist);  %�³�Ա����һ�����ݣ�����һ�δ�������
       
                  count_NB=count_NB+1; %�³�Ա��+1
                  child(numCH,count_NB+1)=col(i);     %!!!���´�ͷ�ĺ���
                  %���¼���Ľڵ�Ϊ��㣬�������ھӵ�NB�������Ų�����©
                  [row,col_2] = find(neighbor_matrix(col(i),:)~=0 ); 
                     if(size(col_2,2)~=0)
                         len_2=size(col_2,2);  %�ھ���
                         for ip=1:len_2
                             infor(col_2(ip)).NB = infor(col_2(ip)).NB-1; %���£�-1
                         end
                     end
              end
          end
          infor(newCH).RE = infor(newCH).RE - count_NB*Erx(64);  %��ͷ����Ҫ�������ݰ��������ܽ����ܺ�
          infor(newCH).RE = infor(newCH).RE - count_NB*Etx1(infor(newCH).dist);  %��ͷ����Ҫ�������ݰ��������ܷ����ܺ�
          %infor(newCH).W = W(infor(newCH).NB,Nmax,infor(newCH).dep,infor(newCH).RE);%�������֮���Wֵ
    end
    
    
     %figure(time);% ��ʾͼƬ 
%%STEP3


times=0;
now_CH=0;
while times<150  %�����ɴ�
   %% ��ȡWֵ���Ľڵ���Ϊ��һ����ͷnewCH
   times=times+1;
   maxW_2=0;
   for i=1:n
      if((infor(i).W>maxW_2) & (infor(i).vis~=0) & (infor(i).type~='CH') & (infor(i).NB>1))
           maxW_2=infor(i).W;
           newCH=i;  %�´�ͷ���newCH
      end
   end
   
   if(now_CH == newCH)
       break      %ֹͣ�ɴػ���
   end
   now_CH = newCH;
   
   numCH=numCH+1;
   infor(newCH).sort=numCH;
   child(numCH,1)=newCH;   %��һ���Ǵ�ͷ���
    %% �õ��´�ͷ֮�󣬸����´�ͷ���ԣ�
    plot(S(newCH).x,S(newCH).y,'ko','MarkerFaceColor','k');     % ���Ƶ�ǰ��ͷ�ڵ�
    infor(newCH).vis = 1;   %̽��״̬Ϊ��̽��
    infor(newCH).type = 'CH';  %���ڸýڵ��Ǵ�ͷ�ڵ㣬��������ΪC
    infor(newCH).clu = 1;  %���ڳɴأ�clu��������Ϊ1
    cluinf(1,newCH) = 1;   %ͬʱ���¾���
    par_CH = infor(newCH).par; %par_CH�����´�ͷ�ĸ��״�ͷ
    infor(newCH).dep = infor(par_CH).dep + 1;  %�´�ͷ���ڵ����Ϊ���״�ͷ�����+1������
    par=infor(newCH).par;
    infor(newCH).dist = sqrt((S(newCH).x - S(par).x)^2 + ( S(newCH).y - S(par).y )^2); %�����´�ͷ�븸�״�ͷ����
    MN_enermax=0; %�ò�������ʣ�����������ݴ����
    IN=0;   %IN����IN�ڵ�ı��
    
    %% ���´�ͷ��NB���γɴأ������³�Ա������Ϣ��
    [row,col] = find(neighbor_matrix(newCH,:)~=0 ); % �ҵ��ýڵ���ھӣ�Ŀ���ǰ������ھӵ�NB����-1����Ϊ�ýڵ�ɴ���
     if(size(col,2)~=0)  %����ýڵ����ھӣ������������ھӵ�NB����-1������δ�ɴ��ھӼ���ô�
          len=size(col,2);  %���ھ���
          count_NB=0;   % count_NB�����³�Ա������ʼ��Ϊ0
          for i=1:len
              
              if(infor(col(i)).clu == 0) %����ھ�δ�ɴأ������ôأ�����ôص�ͬʱҪ�������ھӽڵ��NB
               infor(col(i)).vis = 1;   %̽��״̬Ϊ��̽��   
               plot([S(col(i)).x,S(newCH).x],[S(col(i)).y,S(newCH).y],'bo-');  %�ڵ����ͷ����
               %infor(col(i)).NB=infor(col(i)).NB-1;
               count_NB=count_NB+1; %�³�Ա��+1
               child(infor(newCH).sort,count_NB+1)=col(i);   %�����´صĺ���
            
               infor(col(i)).par=newCH;  %���ø��ף��ôش�ͷ��ΪnewCH
               infor(col(i)).dep=infor(newCH).dep+1;  %�ڵ���ȵ��ڸ������+1
               infor(col(i)).clu=1; %�ɴأ���־Ϊ1
               cluinf(1,col(i)) = 1;   %ͬʱ���¾���
               infor(col(i)).dist = sqrt( (S(col(i)).x - S(newCH).x)^2 + ( S(col(i)).y - S(newCH).y )^2);  %�����³�Ա����ͷ�ľ���
               infor(col(i)).W = W(infor(col(i)).NB,Nmax,infor(col(i)).dep,infor(col(i)).RE);%�����³�Ա��Wֵ
               infor(col(i)).RE = infor(col(i)).RE - Etx1(infor(col(i)).dist);  %�³�Ա����һ�����ݣ�����һ�δ�������
                  %���¼���Ľڵ�Ϊ��㣬�������ھӵ�NB�������Ų�����©   
                  [row,col_2] = find(neighbor_matrix(col(i),:)~=0 ); 
                     if(size(col_2,2)~=0)
                         len=size(col_2,2);  %�³�Ա�ھ���
                         for j=1:len
                             infor(col_2(j)).NB = infor(col_2(j)).NB-1; %���£�-1
                         end
                     end
              end
          end
        
          %���������������
          infor(newCH).RE = infor(newCH).RE - count_NB*Erx(64) - count_NB*Etx1(infor(newCH).dist);  %��ͷ����Ҫ�������ݰ�������֮��Ҫ���״�ͷ����һ�Σ������ܽ����ܺ�
          infor(newCH).W = W(infor(newCH).NB,Nmax,infor(newCH).dep,infor(newCH).RE);%�������֮���Wֵ

          while par~=n+1   %�´���Ϣ���ݸ�SN        
              par_par = infor(par).par; %par_par����Ŀǰ�ڵ�ĸ���ͷ(үү)���������ؼ����
              dist = sqrt((S(par).x - S(par_par).x)^2 + ( S(par).y - S(par_par).y )^2);
              infor(par).RE = infor(par).RE - count_NB*Etx1(dist) - count_NB*Erx(64);%���ܡ�����һ�����ݣ�������������
              %infor(par).W=W(infor(par).NB,Nmax,infor(par).dep,infor(par).RE);%�������֮���Wֵ
              par = infor(par).par;
          end 
          
     end
    %figure(time);% ��ʾͼƬ 
end

  
%% ���濪ʼ����ģ��ѡ����ת������


% ͳ������

for i=1:n
    if(infor(i).type=='CH')
        numCH=numCH+1;   
    end
    if(infor(i).type=='MN')
        numMN=numMN+1;
    end
end


numCH=0;
numCH_num=[];   %��ͷ��ž���
numIN=0;
numIN_num=[];    %IN��ž���
numMN=0;
numMN_num=[];  %MN��ž���

for i=1:n
    copy_reci=[];
    copy_for=[];
    copy_reciMN=[];
    if(infor(i).type=='CH')
        numCH=numCH+1;
        numCH_num(1,numCH)=i;  %�����ͷ��ž���
        infor(i).sort = numCH;   %����ͷ������ţ��������
        ORI(numCH).org=i;      %����¼����¼��Ŷ�Ӧ��ͷ���
        condition(numCH,1)=i;     %��һ�д��
    end

    if(infor(i).type=='MN')
        numMN=numMN+1;
        numMN_num(1,numMN)=i;  %�����ͷ��ž���,�������
        infor(i).sort = numMN;   %��MN������ţ��������
    end
end
fprintf("CH��:%d\n",numCH)
fprintf("MN��:%d\n",numMN) 
%���ȣ�ͳһһ��child�����һ�еı�� ��һ���Ǹ��׵ĵ������ ��ʵ�ʱ��
    for i=1:numCH
        child(i,1)=infor(child(i,1)).sort;
    end
%ȥ��child�����Ŀ��У�501�и��Ƶ�numCH+1����
child(numCH+1,:) = child(n+1,:);
child(numCH+2:end,:)=[];
child(numCH+1,1) = n+1;

%% ���濪ʼ���������еĶ���ڵ㣡����ڵ㣬ʼ�ղ��� ע��ֻ��time=1 ����һ����ʱ ���д˶δ���
if(time==1)
MCH=[];   %��Ŷ���CH
mali_ratio = 0.10;   
mali_num = round(mali_ratio * n) ;  %�����ͷ����
number = 0;
x_cood = [];    %������������
y_cood = [];    %������������
while number~=mali_num        %��Բ�ھ��ȷֲ������ͷ�����㣬����ȷ�������ͷλ��
   rt=rt+1;
   xc=my(1,rt)*xm;
   rt=rt+1;
   yc=my(1,rt)*ym;
   if(sqrt((xc-200)^2+(yc-200)^2 )<r)     %�����Բ�ڣ�����
       x_cood = [x_cood xc];
       y_cood = [y_cood yc];
       number=number+1;
   end
end

for j=1:mali_num
  node_dist=zeros(1,n);     %�������нڵ㵽�����ľ���
  for i=1:n
      %s=numCH_num(1,i);     %sΪCH
      rand_x = x_cood(1,j);
      rand_y = y_cood(1,j);    %ȡ�����������
     % dist = sqrt((S(s).x-rand_x)^2+(S(s).y-rand_y)^2);    %����ô�ͷ����������
      dist = sqrt((S(i).x-rand_x)^2+(S(i).y-rand_y)^2);    %����ýڵ㣨���У�����������
      if (infor(i).pro~='M')
          node_dist(1,i)=dist;   %����ô��׻����Ƕ���ģ��ʹ���ʵ����
      else
          node_dist(1,i)=10000;   %����ô����Ѿ���ѡ���ˣ��ͰѾ�����Ϊ10000�������϶�������С��
      end
  end
  [minminmin,column]=find(node_dist==min(min(node_dist)));
  %infor(numCH_num(1,column)).pro = 'M' ;    %�ô�ͷ��Ϊ���� M
   infor(column).pro = 'M' ;    %�ô�ͷ��Ϊ���� M
  MCH = [MCH column];          %�������ڵ��ž���
end 

% ѡ�����ڵ�֮�󣬸ı���Ϣ
for i=1:mali_num
    x=MCH(1,i);
    plot(S(x).x,S(x).y,'ro','MarkerFaceColor','b');
    if(infor(x).harsh==1 && time==21)%202103�޸�Ϊ21
       infor(x).frate= infor(x).frate;
    else
       rt=rt+1;   %�ı�ת����
       infor(x).frate = 0.8*my(1,rt);  %����ڵ��ת����
       %infor(x).q_return=zeros(9,8);
    end
    infor(x).fratetemp=infor(x).frate;
end
fprintf("����ڵ㱾����ţ�\n");
sort(MCH)


end
  fid7=fopen(['20220815','.txt'],'w');
    for ii3=1:n
        fprintf(fid7,'%d',ii3);
        fprintf(fid7,'%s',' ');
        fprintf(fid7,'%f',infor(ii3).frate);
        fprintf(fid7,'\n');
    end
    fclose(fid7);

%% ����ѡȡIN
SS=0;
temp_ener=0;      %�����ݴ��������ֵ
tempIN=0;   %�����ݴ�һ��IN
for tt=1:n
    if(infor(tt).par==n+1 & infor(tt).type=='MN')
        if(infor(tt).RE>temp_ener)
            tempIN=tt;
            temp_ener=infor(tt).RE;
        end
    end
end
% ѡ�����

if(tempIN~=0)
  infor(tempIN).type='IN';
  INI(tempIN).dep = 1; 
  INI(tempIN).reach = 1;
  INI(tempIN).par = n+1;
  plot(S(tempIN).x,S(tempIN).y,'go','MarkerFaceColor','y');
  infor(n+1).IN = tempIN;
else   %���SN��Χû��ѡ��IN
  infor(n+1).IN = n+1;   
end

 
for i=1:numCH
    MN_enermax=0;
    IN=0;
    CH=numCH_num(1,i); %��ȡ��ͷ��ʵ���
    infor(CH).MNnum=0;
    fakenum=infor(CH).sort;
    for  h=1:numCH    %��child�����У���λ��ͷ����һ�У�
         if(child(h,1) == fakenum)
             SS=h;      %SS��ʾ��λ��������
             break
         end
    end
    j=2;
    while(child(SS,j)~=0)   %�ڸô���ѡ��IN
       if (infor(child(SS,j)).type=='MN' & infor(child(SS,j)).pro=='N' & infor(child(SS,j)).type~='CH')   %IN��MN��ѡȡ
         %Distance_in = sqrt((S(parentin).x-S(in).x)^2+(S(parentin).y-S(in).y)^2);
         if(MN_enermax < infor(child(SS,j)).RE)   %����ó�Ա����������Ѹó�Ա��λIN��
                   MN_enermax = infor(child(SS,j)).RE;   %�����ݴ�ֵ
                   IN=child(SS,j);   %����IN���
         end
       end
       %%���ô˴���ͳ��ÿ�����ڵ�MN��������IN��
       if(infor(child(SS,j)).type=='MN')
           infor(CH).MNnum = infor(CH).MNnum+1;
       end
       j=j+1;
       if(j>size(child,2))   %��ֹ��������ά��
             break;
       end
    end
    %%һ����ͷѡȡ��ϣ�����IN
    if(IN~=0)
        infor(IN).type ='IN';
        infor(CH).IN = IN;
        plot(S(IN).x,S(IN).y,'go','MarkerFaceColor','y');
    end
end
%�ټ�¼һ��BS�ĺ���MN��

jjj=2;
while(child(numCH+1,jjj)~=0)
  
    if(infor(child(numCH+1,jjj)).type=='MN')
      infor(n+1).MNnum=infor(n+1).MNnum+1;
    end
    jjj=jjj+1;
    if(jjj>size(child,2))   %��ֹ��������ά��
             break;
    end
end


for i=1:n
    if(infor(i).type=='IN')
        numIN=numIN+1;
        numIN_num(1,numIN)=i;  %����IN��ž���
    end
end
fprintf("IN��:%d\n",numIN)

cantreach=0;
%% ����IN
for xx=1:n
    INI(xx).dep=0;
    INI(xx).par=n+1;
    INI(xx).reach=1;
    if(infor(xx).type =='IN' & infor(xx).par~=n+1)
       parentCH = infor(xx).par;    %��ȡ��IN�ĸ���CH
       INI(xx).dep = infor(parentCH).dep+1;
       grandparentCH = infor(parentCH).par;   %��ȡ���׵ĸ���CH
       if(grandparentCH==n+1)   %���үү�ǻ�վ������ֱ��ͨ��
           parentIN = n+1;
       else
           parentIN = infor(grandparentCH).IN;
       end   
       INI(xx).par = parentIN;  %IN·�ɵĸ���IN
    elseif( infor(xx).type=='IN' & infor(xx).par==n+1 )
        INI(xx).dep = 1;
        INI(xx).par = n+1;  
    end
end


INI(n+1).par = n+1;
for nn=1:numIN
    inn=numIN_num(1,nn);
    if(INI(inn).par==0)
        pa = infor(inn).par;   %����CH
        grand = infor(pa).par;   %үүCH
        parentin=INI(inn).par;
        while(parentin == 0)
            grand=infor(grand).par;  %����������
            parentin=infor(grand).IN;
            if(parentin~=0)
                INI(inn).par = parentin;  %IN·�ɵĸ���IN
                break
            end
        end
    end
end


%%�ҵ��븸��IN����ͨ�ŵ�IN
for ff=1:numIN
    in=numIN_num(1,ff);    
    parentin = INI(in).par;
    inindis=sqrt((S(parentin).x-S(in).x)^2+(S(parentin).y-S(in).y)^2);
    INI(in).dist=inindis;
    if(inindis>75)
         INI(in).reach=0;
         cantreach=cantreach+1;
    end
end




%% ���濪ʼģ�������
 for i=1+10*(time-1):10+10*(time-1)     %һ�����֣�����20С��/202103�޸�Ϊ10��
 % fprintf("��%d��\n",i);
  if(death==0)  
  for l=1:n        %��MN����ģ�⣬�������ݰ�
      if((infor(l).type == 'MN') & (infor(l).vis~=0) & (infor(l).malicious==0))    %����ýڵ���MN����ʼģ��
        
          Receive(l,i+1)=Receive(l,i+1)+1;  %MN���հ��ɹ���+1����ʵû���հ���MN��ÿ����Ҫ�̶�+1��
          %send=send+1;%����������
          infor(l).receive(1,i+1)=1;
          %SN_receive_packetnum_theory=SN_receive_packetnum_theory+1;
%           if i>201
%             if size(judge_now,2)~=0
%               for djl=1:size(judge_now,2)
%                 if l==judge_now(djl)
%                   [infor(l).q_return,infor(l).frate]=ReinforcementLearning(infor(l).frate,infor(l).q_return,infor(l).pro);
%                 end
%               end
%             end
% %               if infor(l).pro=='M' & ismember(l,judge_now)==0
% %                   [infor(l).q_return,infor(l).frate]=ReinforcementLearning(infor(l).frate,infor(l).q_return,infor(l).pro);
% %               end
%           end
%           if infor(l).harsh==time
%               infor(l).fratetemp=infor(l).frate;  %��������ת����
%               infor(l).frate = infor(l).frate*0.8;  %�ı�ת����
%           end
          rt=rt+1;
          if(1-infor(l).frate)<my(1,rt)   %���������С�������,��ǿ��ѧϰ����������Է����������׿��Խ��յ�
                parent = infor(l).par;         %��ȡ����
                %mm=infor(l).sort;        %��ȡMN�ĵ�����ţ��Ա����ھ���
                Forward(l,i+1)=Forward(l,i+1)+1;  %MN�ķ����ɹ���+1
                infor(l).forward(1,i+1)=1;
                infor(l).RE = infor(l).RE - Etx2(infor(l).dist);  %���Ͱ�������������
                %% ������
                 while (parent~=n+1 & infor(parent).malicious==0 )
                   %yy=infor(parent).sort;      %��ȡ���׵Ĵ�ͷ������ţ��Ա����ھ���
                   Receive(parent,i+1)=Receive(parent,i+1)+1;    %��i+1�У���ΪĿǰ�ǵ�i��
                   infor(parent).receive(1,i+1)=1;
                   infor(parent).RE = infor(parent).RE - Erx(552);  %���׽��ܰ�������������
%                    if i>201
%                      if size(judge_now,2)~=0
%                         for djl=1:size(judge_now,2)
%                             if parent==judge_now(djl)
%                                 [infor(parent).q_return,infor(parent).frate]=ReinforcementLearning(infor(parent).frate,infor(parent).q_return,infor(parent).pro);
%                             end
%                         end
%                      end
% %                         if infor(parent).pro=='M' & ismember(parent,judge_now)==0
% %                              [infor(parent).q_return,infor(parent).frate]=ReinforcementLearning(infor(parent).frate,infor(parent).q_return,infor(parent).pro);
% %                         end
%                    end
%                    if infor(parent).harsh==time
%                         infor(parent).fratetemp=infor(parent).frate;  %��������ת����
%                         infor(parent).frate = infor(parent).frate*0.8;  %�ı�ת����
%                    end
                   rt=rt+1;
                   if (1-infor(parent).frate)<my(1,rt)   %������׵Ķ�����С�������,�����ת���ð�
                          Forward(parent,i+1)=Forward(parent,i+1)+1;
                          infor(parent).forward(1,i+1)=1;
                          infor(parent).RE = infor(parent).RE - Etx2(infor(parent).dist);  %���Ͱ�������������
                          parent=infor(parent).par;    %��ȡ��һ����ͷ
                   else      %����ô�ת��ʧ�ܣ���ֹͣ��������һ��ͷת��
                         break
                   end
         
                 end
          end 
       
      
      elseif(infor(l).type=='IN' & infor(l).par~=n+1)
      %% ����IN·��
              bb=l;
              %%���ȱ���ȡ�س�Ա��packetnum
              parIN=INI(bb).par;    %����IN
              parCH=infor(bb).par;   %���״���
              packetnum = infor(parCH).MNnum -1 +1;  %����������ȥMNnum�е�IN������һ�����ش���1
              infor(bb).RE = infor(bb).RE - packetnum * Etx3(INI(bb).dist);           %�����ܺ�
              while(parIN~=n+1)     %ֻҪû��SN,��һֱ������ȥ
                  infor(parIN).RE = infor(parIN).RE - packetnum * Erx(80);       %�����հ�����
                  infor(parIN).RE = infor(parIN).RE - packetnum * Etx3(INI(parIN).dist); %�����ܺ�
                  parIN=INI(parIN).par;
              end      
      
       elseif(infor(l).type=='IN' & infor(l).par==n+1)   %BS��IN
              bb=l;
              packetnum = infor(parCH).MNnum ;  %������
              infor(bb).RE = infor(bb).RE - packetnum * Etx1(INI(bb).dist);           %�����ܺ�
             
       end
  end
  
  if(size(extra_for,2)~=0)
   for extra=1:size(extra_for,2)   %������Ҫ�෢����MN
      ee=extra_for(extra);
      counte=0;
      if((infor(ee).type == 'MN') & (infor(ee).vis~=0) & (infor(ee).malicious==0))    %����ýڵ���MN����ʼģ��
          while(counte<15)
           Receive(ee,i+1)=Receive(ee,i+1)+1;  %MN���հ��ɹ���+1����ʵû���հ���MN��ÿ����Ҫ�̶�+1��
           %send=send+1;%����������
           infor(ee).receive(1,i+1)=1;
           %SN_receive_packetnum_theory=SN_receive_packetnum_theory+1;
           counte=counte+1;
%            if i>201
%              if size(judge_now,2)~=0
%                 for djl=1:size(judge_now,2)
%                     if ee==judge_now(djl)
%                         [infor(ee).q_return,infor(ee).frate]=ReinforcementLearning(infor(ee).frate,infor(ee).q_return,infor(ee).pro);
%                     end
%                 end
%              end
% %                 if infor(ee).pro=='M' & ismember(ee,judge_now)==0
% %                       [infor(ee).q_return,infor(ee).frate]=ReinforcementLearning(infor(ee).frate,infor(ee).q_return,infor(ee).pro);
% %                 end
%            end
%           if infor(ee).harsh==time
%               infor(ee).fratetemp=infor(ee).frate;  %��������ת����
%               infor(ee).frate = infor(ee).frate*0.8;  %�ı�ת����
%           end
           rt=rt+1;
           if (1-infor(ee).frate)<my(1,rt) %���������С�������,����Է����������׿��Խ��յ�
                parent = infor(ee).par;         %��ȡ����
                Forward(ee,i+1)=Forward(ee,i+1)+1;  %MN�ķ����ɹ���+1
                infor(ee).forward(1,i+1)=1;
                infor(ee).RE = infor(ee).RE - Etx2(infor(ee).dist);  %���Ͱ�������������
                %% ������
                 while (parent~=n+1 & infor(parent).malicious==0 )
                     Receive(parent,i+1)=Receive(parent,i+1)+1;    %��i+1�У���ΪĿǰ�ǵ�i��
                     infor(parent).receive(1,i+1)=1;
                     infor(parent).RE = infor(parent).RE - Erx(552);  %���׽��ܰ�������������
%                      if i>201
%                        if size(judge_now,2)~=0
%                         for djl=1:size(judge_now,2)
%                             if parent==judge_now(djl)
%                                 [infor(parent).q_return,infor(parent).frate]=ReinforcementLearning(infor(parent).frate,infor(parent).q_return,infor(parent).pro);
%                             end
%                         end
%                        end
% %                         if infor(parent).pro=='M' & ismember(parent,judge_now)==0
% %                              [infor(parent).q_return,infor(parent).frate]=ReinforcementLearning(infor(parent).frate,infor(parent).q_return,infor(parent).pro);
% %                         end
%                      end
%                      if infor(parent).harsh==time
%                         infor(parent).fratetemp=infor(parent).frate;  %��������ת����
%                         infor(parent).frate = infor(parent).frate*0.8;  %�ı�ת����
%                      end
                     rt=rt+1;
                     if (1-infor(parent).frate)<my(1,rt)  %������׵Ķ�����С�������,�����ת���ð�
                          Forward(parent,i+1)=Forward(parent,i+1)+1;
                          infor(parent).forward(1,i+1)=1;
                          infor(parent).RE = infor(parent).RE - Etx2(infor(parent).dist);  %���Ͱ�������������
                     end
                     parent=n+1;    %�������궼ǿ�ƽ���
                 end
           end 
          end
      end
   end
  end
%   %% ������Ҫ�෢����MN
  if(size(judge_now,2)~=0 & size(judge_now,1)~=0)
   for xxx=1:judge_num
     obj3=judge_now(xxx);
     if(infor(obj3).type=='MN' & (infor(obj3).vis~=0) & (infor(obj3).malicious==0) )  %�෢����ttl=1
%              obj3;
             packnum2=0;
             while(packnum2<15)
               Receive(obj3,i+1)=Receive(obj3,i+1)+1;  %MN���հ��ɹ���+1����ʵû���հ���MN��ÿ����Ҫ�̶�+1��
               %send=send+1;%����������
               infor(obj3).receive(1,i+1)=1;
               %SN_receive_packetnum_theory=SN_receive_packetnum_theory+1;
               packnum2=packnum2+1;
%                if i>201
%                  if size(judge_now,2)~=0
%                     for djl=1:size(judge_now,2)
%                         if obj3==judge_now(djl)
%                             [infor(obj3).q_return,infor(obj3).frate]=ReinforcementLearning(infor(obj3).frate,infor(obj3).q_return,infor(obj3).pro);
%                         end
%                     end
%                  end
%                     if infor(obj3).pro=='M' & ismember(obj3,judge_now)==0
%                         [infor(obj3).q_return,infor(obj3).frate]=ReinforcementLearning(infor(obj3).frate,infor(obj3).q_return,infor(obj3).pro);
%                     end
%                 end
%                if infor(obj3).harsh==time
%                     infor(obj3).fratetemp=infor(obj3).frate;  %��������ת����
%                     infor(obj3).frate = infor(obj3).frate*0.8;  %�ı�ת����
%                end
               rt=rt+1;
                if (1-infor(obj3).frate)<my(1,rt)  %���������С�������,����Է����������׿��Խ��յ�
                  Forward(obj3,i+1)=Forward(obj3,i+1)+1;  %MN�ķ����ɹ���+1
                  infor(obj3).forward(1,i+1)=1;
                  infor(obj3).RE = infor(obj3).RE - Etx2(infor(obj3).dist);  %���Ͱ�������������
                end
             end
     end
   end
  end
%һ�ֽ������Ѹ��ֶ�Ӧ�����и���һ�����һ��
 if(i<k)
  copy_reci = Receive(:,i+1);
  Receive(:,i+2) = copy_reci;
  copy_for = Forward(:,i+1);
  Forward(:,i+2) = copy_for;
 end
 
 %���㵱��cfr
    for io=1:n
            now_for = Forward(io,i+1);
            now_reci = Receive(io,i+1);
            cfr = now_for/now_reci; 
         if(i<=10) %20210527��20�޸�Ϊ10  
            if (Receive(io,i+1)-Receive(io,i)==0)  %û���հ�
               cfrtemp = CFRtemp(io,i);
            else
               cfrtemp = (Forward(io,i+1)-Forward(io,i))/(Receive(io,i+1)-Receive(io,i));
            end
         else
             symbol=1+10*(time-1);%202103�޸�ÿ10��Ϊһ����
             if (Receive(io,i+1)-Receive(io,symbol)==0)  %û���հ�
               cfrtemp = CFRtemp(io,i);
             else
               cfrtemp = (Forward(io,i+1)-Forward(io,symbol))/(Receive(io,i+1)-Receive(io,symbol));
             end
         end
            CFR(io,i+1)=cfr;
            CFRtemp(io,i+1)=cfrtemp;   %����ת���ʱ���
    end
 
 
 
 if (i==200||i==300||i==400||i==500||i==600||i==700||i==800||i==900||i==1000||i==1100||i==1200)  %%200�ֿ�ʼ��SNÿ�ֶ����ж���ڵ���
      %��ȡ��������
      %now_round=[];
      %now_round(:,1) = CFRtemp(:,i+1);         %��һ�е��ֵ���
%       now_round=[];
%       now_round(:,1) = CFR(:,i);         %��һ�е��ֵ���
%       now_round(:,2) = CFR(:,i+1);  %�ڶ��е����ۼ�
      now_round=zeros(n+1,2);
      for ioo=1:n
          cfrw=CFRtemp(ioo,i-198:i+1);
          cfrh=CFRtemp(ioo,i-198:i+1);
%           now_round(ioo,1) = exp(-std(cfrh,1));        %��һ�е��ֵ���
          if mean(cfrw(:))~=0&infor(ioo).type ~= 'DE'
               now_round(ioo,1) = std(cfrh)/mean(cfrw(:));
               now_round(ioo,2)=CFR(ioo,i+1);  %�ڶ��е����ۼ�
          elseif infor(ioo).type == 'DE'
              now_round(ioo,1)=10;
              now_round(ioo,2)=0;  %�ڶ��е����ۼ�
          else
              now_round(ioo,1)=-1;
              now_round(ioo,2)=CFR(ioo,i+1);  %�ڶ��е����ۼ�
          end
      end
  maxdists=max(now_round(:,1));
  for iwh=1:n
      cfrw=CFRtemp(iwh,i-199:i+1);
      if mean(cfrw(:))==0&infor(iwh).type~='DE'
         now_round(iwh,1) = maxdists;
      end
  end
  now_round(n+1,1)=0;
  now_round(n+1,2)=1;
     %  now_round(501,:)=[0.79 0.79];
     %  now_round(502,:)=[0.74 0.74];
     %  now_round(503,:)=[0.7 0.7];
     % fid=fopen(['forrho','.txt'],'w');%д���ļ�·��  
      fid=fopen(['202010182','.txt'],'w');%д���ļ�·��  
      [r1,c]=size(now_round);% �õ����������������   
  for i1=1:r1
   for j1=i1+1:r1
    dist= sqrt(((now_round(i1,1)-now_round(j1,1))^2+(now_round(i1,2)-now_round(j1,2))^2));
    fprintf(fid,'%d',i1);
    fprintf(fid,'%s',' ');
    fprintf(fid,'%d',j1);
    fprintf(fid,'%s',' ');
    fprintf(fid,'%f',dist);
    fprintf(fid,'\n');
%     if(i1<=500 & j1<=500)
%         fprintf(fid2,'%d',i1);
%         fprintf(fid2,'%s',' ');
%         fprintf(fid2,'%d',j1);
%         fprintf(fid2,'%s',' ');
%         fprintf(fid2,'%f',dist);
%         fprintf(fid2,'\n');
%     end
   end
  end
    fclose(fid);
    dist1=zeros(n,2);
    for i2=1:n+1
        dist1(i2,1)=i2;
        dist1(i2,2)= sqrt(((now_round(i2,1)-0)^2+(now_round(i2,2)-1)^2));
    end
   mmindist=min(dist1(:,2));
   postw=find(dist1(:,2)==mmindist);
   
   % fclose(fid2);
    %[rho_08,rho_075,rho_07]=computerho('forrho.txt');
    [condition,suspect_now,suspect_num,judge_now,judge_num,rho_everjudge]=DBSCAN_NEW('202010182.txt',condition,postw,whzhigh,whzlow,time); 
    %[condition,suspect_now,suspect_num,judge_now,judge_num,rho_everjudge,rho,delta]=DPC_SNnew_desicion_graph('202010182.txt',condition,200,whzhigh,whzlow,ever_judge); 
%     rho=rho';delta=delta';
%     figure (2);plot(rho,delta,'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');hold on
%      for www=1:500
%          if infor(www).pro=='M'
%              plot(rho(1,www),delta(1,www),'o','MarkerSize',2,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on
%          end
%      end
%    % title ('����ǿ��ѧϰ�ڵ�rho','FontSize',15.0)
%     xlabel ('\rho')
%     ylabel ('\delta')
    suspect_nowcopy=suspect_now;  %����һ��
    judge_nowcopy=judge_now;  %����һ��
    fid3=fopen(['20220809','.txt'],'w');
    for ii2=1:n
        fprintf(fid3,'%d',ii2);
        fprintf(fid3,'%s',' ');
        fprintf(fid3,'%f',now_round(ii2,1));
        fprintf(fid3,'%s',' ');
        fprintf(fid3,'%f',now_round(ii2,2));
%         fprintf(fid,'%s',' ');
%         fprintf(fid,'%f',infor(ii2).frate);
        fprintf(fid3,'\n');
    end
    fclose(fid3);
    fid4=fopen(['20220812','.txt'],'w');
    nummh=size(MCH,2);
    for ii3=1:nummh
         fprintf(fid4,'%d',MCH(ii3));
         fprintf(fid4,'\n');
    end
    fclose(fid4);
    fid11=fopen(['20230602','.txt'],'w');
    for zjy=1:suspect_num
        zjyy=suspect_now(1,zjy);
        fprintf(fid11,'%d',zjyy);
        fprintf(fid11,'%s',' ');
        fprintf(fid11,'%f',now_round(zjyy,1));
        fprintf(fid11,'%s',' ');
        fprintf(fid11,'%f',now_round(zjyy,2));
        fprintf(fid11,'\n');
    end
    fclose(fid11);

    %% judge�ڵ㣬����CH�����հ�������MN���෢�� ��IN�Ƕ��⡢������
    for kk=1:size(judge_nowcopy,2)
        if(infor(judge_nowcopy(kk)).type=='IN')
          judge_now(judge_now==judge_nowcopy(kk))=[];   %��IN������
          judge_num=judge_num-1;
        elseif(infor(judge_nowcopy(kk)).malicious~=0)
          judge_now(judge_now==judge_nowcopy(kk))=[];   %������������
          judge_num=judge_num-1;
        end
    end
    %judge�ڵ�������
    for kk2=1:size(judge_now,2)
        if(~ismember(judge_now(kk2),ever_judge))
            ever_judge=[ever_judge judge_now(kk2)];
        end 
    end
    
    %% �۲죬��Ҫjudge�Ľڵ��Ƿ�rho���ڹ���ֵ
    if(size(rho_everjudge,2)~=0)
     for  j4=1:size(rho_everjudge,2)
        if(rho_everjudge(j4)>whzhigh & ~ismember(ever_judge(j4),release_judge))
           release_judge=[release_judge ever_judge(j4)];   %����һ��������ɾ���
        end
     end
    end

    %% ��judge��release_judge�Ľڵ�ɾȥ
    judge_now
    judge_now = setdiff(judge_now,release_judge);
    judge_now
    judge_num = size(judge_now,2);
    judge_nowcopy2=judge_now;  %����һ��
   %%%%%%%���ˣ�����judge�ڵ�̶�%%%%%%
   
   %%%%%%%����ȷ������suspect�ڵ�%%%%%% 
  %ɸ��IN�ڵ�
  for nnn=1:size(suspect_nowcopy,2)       
      isub=suspect_nowcopy(nnn);
      if(infor(isub).type=='IN') %�����IN������
          suspect_now(suspect_now==isub)=[];
          suspect_num=suspect_num-1;
          fprintf('%d�Žڵ㱻�������,��IN\n',isub);
      elseif(Receive(isub,i)<50) %�����ʱ����ΪIN�����������뻷����֪
          suspect_now(suspect_now==isub)=[];
          suspect_num=suspect_num-1;
          fprintf('%d�Žڵ㱻�������\n',isub);
      end
  end
  suspect_nowcopy2=suspect_now;  %�ٸ���һ��
  %�ų������ڵ������
%   for iii=1:size(suspect_nowcopy2,2)
%       subj=suspect_nowcopy2(iii);
%       if now_round(subj,1)/now_round(subj,2)<0.8
%            suspect_now(suspect_now==subj)=[];   %�������
%            suspect_num=suspect_num-1;
%            if(infor(subj).malicious==0)
%              fprintf('%d�Žڵ㱻������ɣ������ǻ�������\n',subj);
%            end
%       end
%   end
%           
  %����suspect �ھӼ�⣬������ھ�ת���ʺܽӽ���˵����harsh������0.7+�Ķ���ڵ㣬���Խ�����ɣ�����ܽӽ�����©�죬������MDR
%     for iii=1:size(suspect_nowcopy2,2)
%         FR_ne=[];
%         subj=suspect_nowcopy2(iii); 
%         subj_ne=nighbor_no(subj,:);   %subj_ne��Ÿö���ڵ���ھӱ�ž���
%         for j=1:size(subj_ne,2)
%             if(subj_ne(j)~=0)
%               FR_ne=[FR_ne now_round(subj_ne(j),2)]; 
%             else
%                 break
%             end
%         end
%         %��FR_ne�������ֵ����׼��
%         FRmean=mean(FR_ne);
%         FRstd=std(FR_ne);
%         subj_FR=now_round(subj,2);   %��ȷ�϶���ڵ�ת���ʣ��ж϶���
%         if(subj_FR-FRmean>0||FRmean-subj_FR<FRstd)
%             suspect_now(suspect_now==subj)=[];   %�������
%             suspect_num=suspect_num-1;
%             if(infor(subj).malicious==0)
%               fprintf('%d�Žڵ㱻������ɣ������ǻ�������\n',subj);
%             end
%         end
%     end
%     
  %�����ھ�ͶƱ�ķ�ʽ������ӻ����е������ڵ�͸�ת���ʵĶ���ڵ�
   for iii=1:size(suspect_nowcopy2,2)
        FR_ne=[];
        subj=suspect_nowcopy2(iii); 
        subj_ne=nighbor_no(subj,:);   %subj_ne��Ÿö���ڵ���ھӱ�ž���
        for j=1:size(subj_ne,2)
            if(subj_ne(j)~=0)
              FR_ne=[FR_ne (sin(5/8*pi*now_round(subj_ne(j),1)+0.5*pi)+1)/now_round(subj_ne(j),2)]; 
            else
                break
            end
        end
        
        %��FR_ne�������ֵ����׼��
%         FRmean=mean(FR_ne);
%         FRstd=std(FR_ne);
        subj_FR=(sin(5/8*pi*now_round(subj,1)+0.5*pi)+1)/now_round(subj,2);
        appr0val=0;
        opp0se=0;
        for vote=1:size(FR_ne,2)%���������ھӣ��ھ�CFR<0.8����֧�֣���֮�򷴶�
            if abs(FR_ne(vote)-subj_FR)<=0.27
                appr0val=appr0val+1;
            else
                opp0se=opp0se+1;
            end
        end
        appr0val_rate=appr0val/(appr0val+opp0se);%����֧����
        %֧���ʴ�����ֵ����ýڵ����Χ�������ơ���ֵ���㣺��������ӻ����ص���/�������
        %�����ӻ������������������ʱ����ֵ=0.391
        if(appr0val_rate>0.391)
            suspect_now(suspect_now==subj)=[];   %�������
            suspect_num=suspect_num-1;
            if(infor(subj).malicious==0)
              fprintf('%d�Žڵ㱻������ɣ������ǻ�������\n',subj);
            end
        end
    end
  %%%%%���ˣ�����suspectȫ��ȷ��%%%%%%
%    www=zeros(size(suspect_nowcopy2,2),2);
%    sumM=zeros(size(suspect_nowcopy2,2),1);
%    nm=size(suspect_nowcopy2,2);
%    neighbor_matrixs=zeros(nm,nm);
%    neighbor_all=[];
% %   R_dist1=zeros(size(suspect_nowcopy2,2),size(suspect_nowcopy2,2)-1);
%    for iii=1:size(suspect_nowcopy2,2)
%        subj=suspect_nowcopy2(iii);
%        www(iii,:)=[S(subj).x,S(subj).y];
%    end
%    for iiii=1:size(suspect_nowcopy2,2)
%        R_dist1 = sqrt((sum(transpose((repmat(www(iiii,:),nm-1,1)-www([1:iiii-1,iiii+1:nm],:)).^2))));
%        neighbor_iiii = R_dist< R/2;
%        neighbor_matrixs(iiii,:) = [neighbor_iiii(1:iiii-1),0,neighbor_iiii(iiii:nm-1)];
% %        sumM(iiii)=sum(sum(R_dist1));
%    end
%     neighbor_nums=zeros(nm,2);
%    for qq=1:nm
%        neighbor_nums(1,qq)=sum(neighbor_matrix(qq,:));
%        subjj=suspect_nowcopy2(qq); 
%        if neighbor_nums(1,qq)>18
%            suspect_now(suspect_now==subjj)=[];   %�������
%            suspect_num=suspect_num-1;
%            if(infor(subjj).malicious==0)
%              fprintf('%d�Žڵ㱻������ɣ������ǻ�������\n',subjj);
%            end
%        end
%    end
%     www=zeros(size(suspect_nowcopy2,2),2);
%    sumM=zeros(size(suspect_nowcopy2,2),1);
%    nm=size(suspect_nowcopy2,2);
%    neighbor_matrixs=zeros(nm,nm);
%    neighbor_all=[];
%    DC=zeros(nm,2);
% %   R_dist1=zeros(size(suspect_nowcopy2,2),size(suspect_nowcopy2,2)-1);
%    for iii=1:size(suspect_nowcopy2,2)
%        subj=suspect_nowcopy2(iii);
%        www(iii,:)=[S(subj).x,S(subj).y];
%    end
%    for iiii=1:size(suspect_nowcopy2,2)
%        R_dist1 = sum(transpose((repmat(www(iiii,:),nm-1,1)-www([1:iiii-1,iiii+1:nm],:)).^2));
%        neighbor_iiii = R_dist1< R/2;
%        neighbor_matrixs(iiii,:) = [neighbor_iiii(1:iiii-1),0,neighbor_iiii(iiii:nm-1)];
%        DC(iiii,1)=iiii;
%        DC(iiii,2)=sum(R_dist1);
%        
% %        sumM(iiii)=sum(sum(R_dist1));
%    end
%     neighbor_nums=zeros(nm,2);
%    for qq=1:nm
%        neighbor_nums(1,qq)=sum(neighbor_matrix(qq,:));
%        subjj=suspect_nowcopy2(qq); 
%        if neighbor_nums(1,qq)>18
%            suspect_now(suspect_now==subjj)=[];   %�������
%            suspect_num=suspect_num-1;
%            if(infor(subjj).malicious==0)
%              fprintf('%d�Žڵ㱻������ɣ������ǻ�������\n',subjj);
%            end
%        end
%    end
%     for qq=1:nm
%         neighbor_nums(qq,2)=sum(neighbor_matrix(qq,:));
%     end
%     for hz=1:5%ȡ������ֵ
%         xxhz=min(min(DC));
%         [x,y]=find(DC==xxhz);%������ֵ��λ��
%         for ziii=1:size(x,1)%�ҵ���������ھ������Ľڵ�
%           zii=x(ziii);
%           subjj=suspect_nowcopy2(zii);
%           suspect_now(suspect_now==subjj)=[];   %�������
%           suspect_num=suspect_num-1;
%           if(infor(subjj).malicious==0)
%             fprintf('%d�Žڵ㱻������ɣ������ǻ�������\n',subjj);
%           end
%           neighbor_all=[neighbor_all,nighbor_no(subjj,:)];%�ھ�֮��
%           DC(zii,2)=max(max(DC));
%         end
%     end
%         
%         
%      for zi=1:nm
%          subjjj=suspect_nowcopy2(zi);
%         if ismember(subjjj, neighbor_all)>0
%            suspect_now(suspect_now==subjjj)=[];   %�������
%            suspect_num=suspect_num-1;
%            if(infor(subjjj).malicious==0)
%              fprintf('%d�Žڵ㱻������ɣ������ǻ�������\n',subjjj);
%            end
%         end
%      end   
  
  
  %% judge suspect�Ѿ�ɸѡһ�飬��¼������Ϣ
  %�ѷ�judge�ڵ��increase��0
  for nn2=1:n       
      if(ismember(nn2,judge_now) == 0 & infor(nn2).malicious == 0) %�������judge�ڵ��У�����0
          infor(nn2).increase=0;
      end
  end
  %�Ѵ�ʱjudge�ڵ�increase��1,conssuspect+1
 if(size(judge_now,2)~=0 & size(judge_now,1)~=0)
  for kk=1:size(judge_now,2)  
      if(infor(judge_now(kk)).malicious==0)
         infor(judge_now(kk)).increase=1;
         infor(judge_now(kk)).conssuspect=infor(judge_now(kk)).conssuspect+1;
      end
  end
 end
  %��¼��Ҫ����CH�෢���Ľڵ�
    extra_for=[];
    for lll=1:numCH
        cc=numCH_num(lll);
      if(infor(cc).increase==1)  %��judgeCH�ĺ��ӣ������Ƕ෢�����ڴ˼�¼
           fff=infor(cc).sort;
           for ppp=1:numCH
               if(child(ppp,1)==fff)
                  break
               end
           end
           lie=2;
           while(child(ppp,lie)~=0)   %��¼��Ҫ���ⷢ����MN
               extra_for=[extra_for child(ppp,lie)];
               lie=lie+1;
                 if(lie>size(child,2))   %��ֹ��������ά��
                   break
                 end
           end
      end
    end
  %����suspect��ڵ�
  for nn=1:n       %�����нڵ�������ѷǻ��ɽڵ��successive��0
      if(ismember(nn,suspect_now) == 0 & infor(nn).malicious == 0) %������ڻ��ɽڵ��У�����0
          infor(nn).successive=0;
      end
  end
  
  for yy=1:suspect_num     %�Ѹ��ֻ��ɽڵ��עΪ1
      zz=suspect_now(1,yy);
      if(infor(zz).malicious == 0)      
         infor(zz).successive = infor(zz).successive+1;
      end
  end
  
   for mm=1:n      %�۲��Ƿ���ĳ�ڵ��������ֱ�����
      if(infor(mm).successive == 1 & infor(mm).malicious==0 & infor(mm).vis~=0)    %����ýڵ��������ֱ�����(���һ�δ��ȷ��)����ȷ�϶���
          infor(mm).malicious = i;    %��¼���ڵڼ��ֱ��������ģ�
          infor(mm).type = 'DE';
          infor(mm).locationx = S(mm).x;
          infor(mm).locationy = S(mm).y;
          S(mm).x=500;
          S(mm).y=500;
          confirm=[confirm mm];    %�Ѹ�ȷ�ϵ�ӵ�confirm������
       elseif(infor(mm).conssuspect == 25 & infor(mm).malicious==0 & infor(mm).vis~=0)
           infor(mm).malicious = i;    %��¼���ڵڼ��ֱ��������ģ�
           infor(mm).type = 'DS';
           infor(mm).locationx = S(mm).x;
           infor(mm).locationy = S(mm).y;
           S(mm).x=500;
           S(mm).y=500;
           confirm=[confirm mm];    %�Ѹ�ȷ�ϵ�ӵ�confirm������
      end
      
      
   end
  miss=setdiff(MCH,confirm);    %���ɽڵ��л��ɴ�ģ�Ҳ����©��
  fault=setdiff(confirm,MCH) ;    %���ɽڵ�û�ҳ����ģ�Ҳ�������
  misscopy=miss; %����һ��  
   if size(fault,2)~=0
   fprintf("���ڵ㣺")
   fault
   FDR = size(fault,2)/(n-mali_num);   %����FDR
  else
   FDR = 0;   %����FDR
  end
  if size(miss,2)~=0
   fprintf("©��ڵ㣺")
   miss
   MDR = size(miss,2)/mali_num;   %����MDR
  else
   MDR = 0;   %����FDR
  end
  FDRM=[FDRM,FDR];
  MDRM=[MDRM,MDR];
  end
  
 %% ÿ�ּ���Ƿ���һ���ڵ�����������У�death=1
 for vv=1:n
     % ��¼���нڵ�����
     %energy(vv,i)=infor(vv).RE;
     if(infor(vv).RE<=0)
         vv
         death=1;
         fprintf('$$$$$$$$$$$$$$')
         deathtime=[deathtime i];
         if(i>225)   %�����˼��
             flag=1;
         end
         %break
     end
 end
 
 if(death==1)
     break
 end
 
 end
 

%plotmali(S,time,infor)    
    
end
 %jilu=[jilu i];
if(flag==1)  %�˴η������
  simulation=simulation+1;
  confirm
  miss=setdiff(MCH,confirm);    %���ɽڵ��л��ɴ�ģ�Ҳ����©��
  fault=setdiff(confirm,MCH) ;    %���ɽڵ�û�ҳ����ģ�Ҳ�������
  misscopy=miss; %����һ��
%   fid10=fopen(['FDR and MDR','.txt'],'w');
%   for www=1:29
%       fprintf(fid7,'%f',FDRM(www));
%       fprintf(fid7,'%s',' ');
%       fprintf(fid7,'%f',MDRM(www));
%       fprintf(fid7,'\n');
%   end
%   fclose(fid10);
  %% ��©��ڵ�ɸѡ
%   for i=1:size(misscopy,2)
%         FRm_ne=[];
%         subjm=misscopy(i); 
%         subjm_ne=nighbor_no(subjm,:);   %subj_ne��Ÿö���ڵ���ھӱ�ž���
%         for j=1:size(subjm_ne,2)
%             if(subjm_ne(j)~=0)
%               FRm_ne=[FRm_ne now_round(subjm_ne(j),2)]; 
%             else
%                 break
%             end
%         end
%         %��FRm_ne�������ֵ����׼��
%         FRmm=mean(FRm_ne);
%         FRstdm=std(FRm_ne);
%         subjm_FR=now_round(subjm,2);   
%         if(abs(subjm_FR-FRmm)<2*FRstdm)
%             miss(miss==subjm)=[];   %������MDR
%             fprintf('%d�Žڵ㲻��MDR\n',subjm);
%         end
%    end
  
  
  
  if size(fault,2)~=0
   fprintf("���ڵ㣺")
   fault
   FDR = size(fault,2)/(n-mali_num);   %����FDR
  else
   FDR = 0;   %����FDR
  end
  if size(miss,2)~=0
   fprintf("©��ڵ㣺")
   miss
   MDR = size(miss,2)/mali_num;   %����MDR
  else
   MDR = 0;   %����FDR
  end
  Accuracy=(n-size(fault,2)-size(miss,2))/n;%��������ȷ��
  %����throughput
%   for i=1:n
%       if (infor(i).par==501)
%           SN_receive_packetnum_reality=SN_receive_packetnum_reality+Forward(i);
%       end
%   end
  %throughput=SN_receive_packetnum_reality;

  receive_sum=sum(Receive);
  forward_sum=sum(Forward);
  lastRE=1;
 lastFO=1;
  while(lastRE~=k)
      if (receive_sum(lastRE)==0)
          receive_final=lastRE-1;
          break
      else
          lastRE=lastRE+1;
      end
  end
  while(lastFO~=k)
      if (forward_sum(lastFO)==0)
          forward_final=lastFO-1;
          break
      else
          lastFO=lastFO+1;
      end
  end
  throughput_rate=forward_sum(forward_final)/receive_sum(receive_final);
% child_id=child(end,:);
% child_id=child_id(child_id~=0);
% for dd=2:size(child_id,2)
%     forward_id=Receive(child_id(dd),:);
%     forward_id=forward_id(forward_id~=0);
%     tu=tu+forward_id(end);
% end
% for ddd=1:500
%     if infor(ddd).type=='MN'
%         receive_id=Receive(ddd,:);
%         receive_id=receive_id(receive_id~=0);
%         tun=tun+receive_id(end);
%     end
% end
% throughput_rate=tu/tun;
          
  FDRmatix(1,simulation)=FDR;
  MDRmatix(1,simulation)=MDR;
 fprintf("FDR:%f\n",FDR);
 fprintf("MDR:%f\n",MDR);
 zongjie_simu(simulation).receive=Receive;
 zongjie_simu(simulation).forward=Forward;
 zongjie_simu(simulation).infor=infor;
 zongjie_simu(simulation).CFR=CFR;
 zongjie_simu(simulation).FDR=FDR;
 zongjie_simu(simulation).MDR=MDR;
 zongjie_simu(simulation).fault_num=size(fault,2);
 zongjie_simu(simulation).CHnum=numCH;
 zongjie_simu(simulation).miss_num=size(miss,2);
 zongjie_simu(simulation).accuracy=Accuracy;
 zongjie_simu(simulation).throughput_rate=throughput_rate;
   zongjie(simulation,1)=FDR;%��һ��FDR
   zongjie(simulation,2)=MDR; %�ڶ���MDR
%   zongjie(simulation,3)=size(fault,2);%���������ڵ����
%   zongjie(simulation,4)=numCH;%�����д�ͷ����
%   zongjie(simulation,5)=size(miss,2);%������©��ڵ����
  zongjie(simulation,3)=Accuracy;%�����м�������ȷ��
%   %�����еڰ�������������
%   %zongjie(simulation,7)=throughput;
   zongjie(simulation,4)=throughput_rate;
  
end
  
end
end
 
 zongjie
 %fprintf("�����ֵ����һ��fdr �ڶ���mdr��\n");
% whzzongjie_mean(whzlow-81,:)=100*mean(zongjie,1);
whzzongjie_mean(whzlow-49,whzhigh-235:whzhigh-232)=100*mean(zongjie,1);
%whzzongjie(whzlow-81,whzhigh-177:whzhigh-176)=zongjie_simu;
wxc(whzlow-49).zongjie=zongjie_simu;
 %fprintf("��������һ��fdr �ڶ���mdr��\n");
 % whzzongjie_var(whzlow-81,:)=var(zongjie);
 whzzongjie_var(whzlow-49,whzhigh-235:whzhigh-232)=var(zongjie);
 %whzzongjie_var(whzlow-81,whzhigh-171:whzhigh-165)=var(zongjie);
    end
end

  fprintf("�����ֵ����һ��FDR���ڶ���MDR����������ȷ�ʣ������У�throughput_rate��\n");
  whzzongjie_mean
  fprintf("��������һ��FDR���ڶ���MDR����������ȷ�ʣ������У�throughput_rate��\n");
  whzzongjie_var

for i=1:n
figure (2);plot(i,infor(i).frate,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');hold on
if infor(i).pro=='M'
plot(i,infor(i).frate,'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on
end
end
plot(16,infor(16).frate,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(227,infor(227).frate,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(255,infor(255).frate,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(345,infor(345).frate,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(14,infor(14).frate,'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on
plot(165,infor(165).frate,'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on
plot(185,infor(185).frate,'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on
plot(410,infor(410).frate,'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on
xlabel('Nodes ID');
ylabel('Average forwarding rate (%)');
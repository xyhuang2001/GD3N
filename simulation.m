 clear;close all;clc;
whzzongjie_mean=[];
 whzzongjie_var=[];
for whzhigh=230   
    for whzlow=50
zongjie=[];
jilu=[];
deathtime=[];
%% 函数定义
   W=inline('0.5*(x/y)+0.3*(1/z)+0.2*(p/0.52)','x','y','z','p');  
   % x代表该节点未成簇邻居数，y代表x最大值，z代表该节点深度，p代表该节点剩余能量  
   %%Etx代表一个节点传输出去数据所消耗的能量
   
   Eelec=50e-9; %Eelec
   Efs=10e-12;%系数1
   Emp=0.0013e-13;%系数2
   Etx1=inline('64*50e-9+64*10e-12*d*d','d');   %Etx1 代表64bit数据包，为成簇信息包
   Etx2=inline('552*50e-9+552*10e-12*d*d','d');   %Etx2 代表552bit数据包，为数据信息包
   Etx3=inline('80*50e-9+80*10e-12*d*d','d');   %Etx2 代表80bit数据包，为IN信息包
%%Erx代表一个节点要是接受数据，消耗的能量
%%Erx代表一个节点要是接受数据，消耗的能量

    Erx=inline('k*50e-9','k'); %一元函数，k代表数据包比特数
    
%%一深度节点每一轮消耗能量函数定义，该函数决定系统寿命

   Eperround=inline('(500-x)*k*50e-9 + 500*y','x','k','y'); 
   %x代表一深度节点数，k代表数据包比特数，y代表Etx(k,d)
   simulation=0;
   
while(simulation~=10)

%% 网络的初始化，包括节点的布局与信息初始化
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
% 内500个节点，通信半径75m
FDRM=[];
MDRM=[];
nn=0;
R=75;
E0=1.7;  %初始能量2J
numCH=0;   %簇头数量
child=[];    %child矩阵用于保存各个簇头的孩子编号
energy=[];
n=500;
%SN_receive_packetnum_theory=0;%SN收包理论数
%SN_receive_packetnum_reality=0;%SN收包实际数
%网络吞吐量
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
plot( x_BS, y_BS, 'd','markerfacecolor',[1,0,0]); %BS绘制
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
            infor(i).frate = 0.8+0.2*my(1,rt); %给每个节点初始化一个正常转发率
            infor(i).malicious = 0;        
            infor(i).RE = E0;    %infor(i).RE表示i节点的剩余能量，初始设置为2J
            infor(i).pro = 'N';    %代表是否恶意，全程不会改变，N为正常，M为恶意
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
%% 布局结束，开始模拟
    k=5000;
    Receive=zeros(n,k+1);    %初始化簇头收包数矩
    Receive(:,2)=2;
    Forward=zeros(n,k+1);    %初始化簇头发包数矩阵
    Forward(:,2)=2;
    CFR=zeros(n,k+1);        %初始化累计转发率矩阵
    CFRtemp=zeros(n,k+1);        %初始化累计转发率矩阵
    
while(death==0)     %time代表大轮
    time=time+1
    %uu=plotnetwork(S,time,infor);
    %% 首先还原
  for w=2:200
    if(time==w)
        for orig=1:n
            if(infor(orig).harsh==time-1)
                infor(orig).frate = infor(orig).fratetemp;%还原
            end
        end
    end
  end
    %% 模拟部分区域信道质量差
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
          if(sqrt((S(har).x-xharsh)^2+(S(har).y-yharsh)^2)<37.5)  %harsh区域
                 plot(S(har).x,S(har).y,'ro','MarkerFaceColor','r');
                 infor(har).fratetemp=infor(har).frate;  %保存正常转发率
                 infor(har).frate = infor(har).frate*0.8;  %改变转发率
                 %infor(har).frate = 0.15;  %改变转发率
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
   %% 成簇
           c=zeros(n,2);
           for i=1:n  %初始化第一部分 infor代表节点属性
              infor(i).clu = 0;  %clu属性代表该节点未成簇，成簇则为1
              infor(i).dep = 0;  %dep表示该节点深度，初始时均设置为0
              c(i,:)=[S(i).x,S(i).y];  
              Receive(i,1)=i;
              Forward(i,1)=i;
              CFR(i,1)=i;
           end
           %获得目前无线传感器网络的邻接矩阵neighbor_matrix
           for i=1:n
              R_dist = sqrt((sum(transpose((repmat(c(i,:),n-1,1)-c([1:i-1,i+1:n],:)).^2))));
              neighbor_i = R_dist< R/2 & infor(i).clu == 0;  %通信范围内+未成簇，设置为1，即论文中的未成簇邻居
              neighbor_matrix(i,:) = [neighbor_i(1:i-1),0,neighbor_i(i:n-1)];
           end
           %下面计算所有节点的邻居数
           neighbor_num=zeros(1,n);
           for q=1:n
               neighbor_num(1,q)=sum(neighbor_matrix(q,:));  %存储每个点的邻居数
           end
           %Nmax=max(max(neighbor_num));  %Nmax是NB最大值，在此获取
           Nmax=14;
           %fprintf('Nmax%d\n',Nmax)
           %下面统计每个节点的邻居编号
           nighbor_no=zeros(n,40);
           for d=1:n
               tempn=neighbor_matrix(d,:);   %提取特定行
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
           
           for i=1:n %初始化第二部分
               infor(i).NB = sum(neighbor_matrix(i,:)~=0,2);   % infor(i).NB代表i节点的未成簇邻居数，初始等于所有邻居数
               infor(i).vis = 0;   %探测状态为未探测
               infor(i).par = 0;  %初始时，所有节点的父亲设置为0
               infor(i).MNnum = 0;
               infor(i).W = 0;
               if(infor(i).type ~= 'DE')
                 infor(i).type = 'MN';  %初始时，节点类型均设置为MN
               end
               %infor(i).W = 0.5*(neighbor_num(1,i)/Nmax)+0.2;  %初始W
           end
      
      
   MinDist_Node_BS = zeros(1,n);  %先设置一个矩阵 为了后续存放数据
   infor(n+1).dep = 0;
   infor(n+1).MNnum = 0;
   S(n+1).x=200;
   S(n+1).y=200;
   

 %% STEP1 STEP2
 BS_NB=0;  %记录BS的孩子总数
   for i=1:n
      MinDist_Node_BS(i) = sqrt( (S(i).x-200)^2 + ( S(i).y-200 )^2  );  %计算其余所有节点到BS的距离
      % 如果该节点与BS距离小于R，则划入第一个簇,并进行操作
      if(MinDist_Node_BS(i)<R/2 & infor(i).malicious==0)
           plot([S(i).x,200],[S(i).y,200],'bo-');  %节点与簇头连线
           infor(i).vis = 1;   %探测状态为已探测
           infor(i).clu = 1;  %由于成簇，clu属性设置为1
           %cluinf(1,i) = 1;   %同时更新矩阵
           infor(i).dep = 1;  %第一个簇，节点深度为1
           infor(i).par = n+1;  %父亲是BS，用501表示
           infor(i).dist = MinDist_Node_BS(i);  %计算新成员到簇头的距离
           BS_NB=BS_NB+1; %新成员数+1
           child(n+1,BS_NB+1)=i;      %保存基站的孩子
           infor(i).W = W(infor(i).NB,Nmax,infor(i).dep,infor(i).RE);
           infor(i).RE = infor(i).RE - Etx1(infor(i).dist);%传输一次数据，消耗一次能量
           [row,col] = find(neighbor_matrix(i,:)~=0 ); % 找到该节点的邻居，目的是把他的邻居的NB属性-1，因为该节点成簇了
           if(size(col,2)~=0)  %如果该节点有邻居，把他的邻居的NB属性-1
             len=size(col,2);  %邻居数
             for j=1:len
                 infor(col(j)).NB=infor(col(j)).NB-1;
             end
           end
           % 然后BS计算第一个簇内每个节点的W值
           infor(i).dep=1;  
      end    
   end
   
   %figure(time);% 显示图片 
  %%STEP2
    %获取W值最大的节点作为第一个簇头
    maxW=0;
    newCH=0;
    
   for i=1:n
      if((infor(i).W>maxW) & (infor(i).vis~=0) & (infor(i).type~='CH') & infor(i).NB>1)
           maxW=infor(i).W;
           newCH=i;
      end
   end
   
   if(newCH==0)
       fprintf('未找到第一个簇头！');
   end
   
   numCH=numCH+1;    %簇头数加一
   infor(numCH).sort = newCH;    %簇头新编号
   child(numCH,1)=newCH;   %第一列是簇头编号
    %得到第一个簇头之后，更新新簇头属性
    zht1=plot(S(newCH).x,S(newCH).y,'ko','MarkerFaceColor','k');     % 绘制当前簇头节点
    
    infor(newCH).vis = 1;   %探测状态为已探测
    infor(newCH).type = 'CH';  %由于该节点是簇头节点，类型设置为C
    infor(newCH).clu = 1;  %由于成簇，clu属性设置为1
    cluinf(1,newCH) = 1;   %同时更新矩阵
    infor(newCH).dep = 1;  %簇头，节点深度为1
    [row,col] = find(neighbor_matrix(newCH,:)~=0 ); % 找到该节点的邻居，目的是把他的邻居的NB属性-1，因为该节点成簇了
    if(size(col,2)~=0)  %如果该节点有邻居，把他的邻居的NB属性-1，并且未成簇邻居加入该簇
          len=size(col,2);  %邻居数
          count_NB=0;  % count_NB代表新成员数，初始化为0
          MN2_enermax = 0; 
          IN2 = 0;
          for i=1:len
              if(infor(col(i)).clu == 0 & infor(col(i)).malicious==0) %如果邻居未成簇，则加入该簇，加入该簇的同时要更新其邻居节点的NB 
                  infor(col(i)).vis = 1;   %探测状态为已探测
                  plot([S(col(i)).x,S(newCH).x],[S(col(i)).y,S(newCH).y],'bo-');  %节点与簇头连线
                  infor(col(i)).par=newCH;  %设置父亲（该簇簇头）为newCH
                  infor(col(i)).dep=infor(newCH).dep+1;  %节点深度等于父亲深度+1
                  infor(col(i)).clu=1; %成簇，标志为1
                  cluinf(1,col(i)) = 1;   %同时更新矩阵
                  infor(col(i)).dist = sqrt( (S(col(i)).x - S(newCH).x)^2 + ( S(col(i)).y - S(newCH).y )^2);  %计算新成员到簇头的距离
                  infor(col(i)).W=W(infor(col(i)).NB,Nmax,infor(col(i)).dep,infor(col(i)).RE);%计算新成员的W值
                  infor(col(i)).RE = infor(col(i)).RE - Etx1(infor(col(i)).dist);  %新成员传输一次数据，消耗一次传输能量
       
                  count_NB=count_NB+1; %新成员数+1
                  child(numCH,count_NB+1)=col(i);     %!!!存新簇头的孩子
                  %以新加入的节点为起点，更新其邻居的NB，这样才不会遗漏
                  [row,col_2] = find(neighbor_matrix(col(i),:)~=0 ); 
                     if(size(col_2,2)~=0)
                         len_2=size(col_2,2);  %邻居数
                         for ip=1:len_2
                             infor(col_2(ip)).NB = infor(col_2(ip)).NB-1; %更新，-1
                         end
                     end
              end
          end
          infor(newCH).RE = infor(newCH).RE - count_NB*Erx(64);  %簇头由于要接受数据包，计算总接受能耗
          infor(newCH).RE = infor(newCH).RE - count_NB*Etx1(infor(newCH).dist);  %簇头由于要发送数据包，计算总发送能耗
          %infor(newCH).W = W(infor(newCH).NB,Nmax,infor(newCH).dep,infor(newCH).RE);%计算更新之后的W值
    end
    
    
     %figure(time);% 显示图片 
%%STEP3


times=0;
now_CH=0;
while times<150  %持续成簇
   %% 获取W值最大的节点作为新一个簇头newCH
   times=times+1;
   maxW_2=0;
   for i=1:n
      if((infor(i).W>maxW_2) & (infor(i).vis~=0) & (infor(i).type~='CH') & (infor(i).NB>1))
           maxW_2=infor(i).W;
           newCH=i;  %新簇头编号newCH
      end
   end
   
   if(now_CH == newCH)
       break      %停止成簇环节
   end
   now_CH = newCH;
   
   numCH=numCH+1;
   infor(newCH).sort=numCH;
   child(numCH,1)=newCH;   %第一列是簇头编号
    %% 得到新簇头之后，更新新簇头属性：
    plot(S(newCH).x,S(newCH).y,'ko','MarkerFaceColor','k');     % 绘制当前簇头节点
    infor(newCH).vis = 1;   %探测状态为已探测
    infor(newCH).type = 'CH';  %由于该节点是簇头节点，类型设置为C
    infor(newCH).clu = 1;  %由于成簇，clu属性设置为1
    cluinf(1,newCH) = 1;   %同时更新矩阵
    par_CH = infor(newCH).par; %par_CH代表新簇头的父亲簇头
    infor(newCH).dep = infor(par_CH).dep + 1;  %新簇头，节点深度为父亲簇头的深度+1！！！
    par=infor(newCH).par;
    infor(newCH).dist = sqrt((S(newCH).x - S(par).x)^2 + ( S(newCH).y - S(par).y )^2); %计算新簇头与父亲簇头距离
    MN_enermax=0; %该参数代表剩余能量，是暂存变量
    IN=0;   %IN代表IN节点的编号
    
    %% 找新簇头的NB，形成簇，簇内新成员更新信息：
    [row,col] = find(neighbor_matrix(newCH,:)~=0 ); % 找到该节点的邻居，目的是把他的邻居的NB属性-1，因为该节点成簇了
     if(size(col,2)~=0)  %如果该节点有邻居，把他的所有邻居的NB属性-1，并且未成簇邻居加入该簇
          len=size(col,2);  %总邻居数
          count_NB=0;   % count_NB代表新成员数，初始化为0
          for i=1:len
              
              if(infor(col(i)).clu == 0) %如果邻居未成簇，则加入该簇，加入该簇的同时要更新其邻居节点的NB
               infor(col(i)).vis = 1;   %探测状态为已探测   
               plot([S(col(i)).x,S(newCH).x],[S(col(i)).y,S(newCH).y],'bo-');  %节点与簇头连线
               %infor(col(i)).NB=infor(col(i)).NB-1;
               count_NB=count_NB+1; %新成员数+1
               child(infor(newCH).sort,count_NB+1)=col(i);   %保存新簇的孩子
            
               infor(col(i)).par=newCH;  %设置父亲（该簇簇头）为newCH
               infor(col(i)).dep=infor(newCH).dep+1;  %节点深度等于父亲深度+1
               infor(col(i)).clu=1; %成簇，标志为1
               cluinf(1,col(i)) = 1;   %同时更新矩阵
               infor(col(i)).dist = sqrt( (S(col(i)).x - S(newCH).x)^2 + ( S(col(i)).y - S(newCH).y )^2);  %计算新成员到簇头的距离
               infor(col(i)).W = W(infor(col(i)).NB,Nmax,infor(col(i)).dep,infor(col(i)).RE);%计算新成员的W值
               infor(col(i)).RE = infor(col(i)).RE - Etx1(infor(col(i)).dist);  %新成员传输一次数据，消耗一次传输能量
                  %以新加入的节点为起点，更新其邻居的NB，这样才不会遗漏   
                  [row,col_2] = find(neighbor_matrix(col(i),:)~=0 ); 
                     if(size(col_2,2)~=0)
                         len=size(col_2,2);  %新成员邻居数
                         for j=1:len
                             infor(col_2(j)).NB = infor(col_2(j)).NB-1; %更新，-1
                         end
                     end
              end
          end
        
          %下面进行能量清算
          infor(newCH).RE = infor(newCH).RE - count_NB*Erx(64) - count_NB*Etx1(infor(newCH).dist);  %簇头由于要接受数据包，并且之后还要向父亲簇头传输一次，计算总接受能耗
          infor(newCH).W = W(infor(newCH).NB,Nmax,infor(newCH).dep,infor(newCH).RE);%计算更新之后的W值

          while par~=n+1   %新簇信息传递给SN        
              par_par = infor(par).par; %par_par代表目前节点的父簇头(爷爷)，方便计算簇间距离
              dist = sqrt((S(par).x - S(par_par).x)^2 + ( S(par).y - S(par_par).y )^2);
              infor(par).RE = infor(par).RE - count_NB*Etx1(dist) - count_NB*Erx(64);%接受、传输一次数据，消耗两种能量
              %infor(par).W=W(infor(par).NB,Nmax,infor(par).dep,infor(par).RE);%计算更新之后的W值
              par = infor(par).par;
          end 
          
     end
    %figure(time);% 显示图片 
end

  
%% 下面开始进行模拟选择性转发攻击


% 统计数量

for i=1:n
    if(infor(i).type=='CH')
        numCH=numCH+1;   
    end
    if(infor(i).type=='MN')
        numMN=numMN+1;
    end
end


numCH=0;
numCH_num=[];   %簇头编号矩阵
numIN=0;
numIN_num=[];    %IN编号矩阵
numMN=0;
numMN_num=[];  %MN编号矩阵

for i=1:n
    copy_reci=[];
    copy_for=[];
    copy_reciMN=[];
    if(infor(i).type=='CH')
        numCH=numCH+1;
        numCH_num(1,numCH)=i;  %存入簇头编号矩阵
        infor(i).sort = numCH;   %给簇头单独编号，方便后面
        ORI(numCH).org=i;      %反记录，记录编号对应簇头编号
        condition(numCH,1)=i;     %第一列存放
    end

    if(infor(i).type=='MN')
        numMN=numMN+1;
        numMN_num(1,numMN)=i;  %存入簇头编号矩阵,方便后面
        infor(i).sort = numMN;   %给MN单独编号，方便后面
    end
end
fprintf("CH数:%d\n",numCH)
fprintf("MN数:%d\n",numMN) 
%首先，统一一下child矩阵第一列的编号 第一列是父亲的单独编号 非实际编号
    for i=1:numCH
        child(i,1)=infor(child(i,1)).sort;
    end
%去掉child矩阵后的空行，501行复制到numCH+1行上
child(numCH+1,:) = child(n+1,:);
child(numCH+2:end,:)=[];
child(numCH+1,1) = n+1;

%% 下面开始决定网络中的恶意节点！恶意节点，始终不变 注意只有time=1 即第一大轮时 运行此段代码
if(time==1)
MCH=[];   %存放恶意CH
mali_ratio = 0.10;   
mali_num = round(mali_ratio * n) ;  %恶意簇头数量
number = 0;
x_cood = [];    %存放随机横坐标
y_cood = [];    %存放随机纵坐标
while number~=mali_num        %向圆内均匀分布恶意簇头数个点，用于确定恶意簇头位置
   rt=rt+1;
   xc=my(1,rt)*xm;
   rt=rt+1;
   yc=my(1,rt)*ym;
   if(sqrt((xc-200)^2+(yc-200)^2 )<r)     %如果在圆内，则保留
       x_cood = [x_cood xc];
       y_cood = [y_cood yc];
       number=number+1;
   end
end

for j=1:mali_num
  node_dist=zeros(1,n);     %保存所有节点到随机点的距离
  for i=1:n
      %s=numCH_num(1,i);     %s为CH
      rand_x = x_cood(1,j);
      rand_y = y_cood(1,j);    %取出随机点坐标
     % dist = sqrt((S(s).x-rand_x)^2+(S(s).y-rand_y)^2);    %计算该簇头到随机点距离
      dist = sqrt((S(i).x-rand_x)^2+(S(i).y-rand_y)^2);    %计算该节点（所有）到随机点距离
      if (infor(i).pro~='M')
          node_dist(1,i)=dist;   %如果该簇首还不是恶意的，就存真实距离
      else
          node_dist(1,i)=10000;   %如果该簇首已经被选过了，就把距离设为10000，让他肯定不会最小了
      end
  end
  [minminmin,column]=find(node_dist==min(min(node_dist)));
  %infor(numCH_num(1,column)).pro = 'M' ;    %该簇头设为恶意 M
   infor(column).pro = 'M' ;    %该簇头设为恶意 M
  MCH = [MCH column];          %存入恶意节点编号矩阵
end 

% 选完恶意节点之后，改变信息
for i=1:mali_num
    x=MCH(1,i);
    plot(S(x).x,S(x).y,'ro','MarkerFaceColor','b');
    if(infor(x).harsh==1 && time==21)%202103修改为21
       infor(x).frate= infor(x).frate;
    else
       rt=rt+1;   %改变转发率
       infor(x).frate = 0.8*my(1,rt);  %恶意节点的转发率
       %infor(x).q_return=zeros(9,8);
    end
    infor(x).fratetemp=infor(x).frate;
end
fprintf("恶意节点本来编号：\n");
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

%% 簇内选取IN
SS=0;
temp_ener=0;      %用于暂存最大能量值
tempIN=0;   %用于暂存一跳IN
for tt=1:n
    if(infor(tt).par==n+1 & infor(tt).type=='MN')
        if(infor(tt).RE>temp_ener)
            tempIN=tt;
            temp_ener=infor(tt).RE;
        end
    end
end
% 选择完毕

if(tempIN~=0)
  infor(tempIN).type='IN';
  INI(tempIN).dep = 1; 
  INI(tempIN).reach = 1;
  INI(tempIN).par = n+1;
  plot(S(tempIN).x,S(tempIN).y,'go','MarkerFaceColor','y');
  infor(n+1).IN = tempIN;
else   %如果SN周围没能选出IN
  infor(n+1).IN = n+1;   
end

 
for i=1:numCH
    MN_enermax=0;
    IN=0;
    CH=numCH_num(1,i); %获取簇头真实编号
    infor(CH).MNnum=0;
    fakenum=infor(CH).sort;
    for  h=1:numCH    %在child矩阵中，定位簇头的那一行！
         if(child(h,1) == fakenum)
             SS=h;      %SS表示定位到的行数
             break
         end
    end
    j=2;
    while(child(SS,j)~=0)   %在该簇内选举IN
       if (infor(child(SS,j)).type=='MN' & infor(child(SS,j)).pro=='N' & infor(child(SS,j)).type~='CH')   %IN从MN中选取
         %Distance_in = sqrt((S(parentin).x-S(in).x)^2+(S(parentin).y-S(in).y)^2);
         if(MN_enermax < infor(child(SS,j)).RE)   %如果该成员能量更大，则把该成员定位IN点
                   MN_enermax = infor(child(SS,j)).RE;   %更新暂存值
                   IN=child(SS,j);   %更新IN编号
         end
       end
       %%利用此处，统计每个簇内的MN数（不含IN）
       if(infor(child(SS,j)).type=='MN')
           infor(CH).MNnum = infor(CH).MNnum+1;
       end
       j=j+1;
       if(j>size(child,2))   %防止超出矩阵维度
             break;
       end
    end
    %%一个簇头选取完毕，绘制IN
    if(IN~=0)
        infor(IN).type ='IN';
        infor(CH).IN = IN;
        plot(S(IN).x,S(IN).y,'go','MarkerFaceColor','y');
    end
end
%再记录一下BS的孩子MN数

jjj=2;
while(child(numCH+1,jjj)~=0)
  
    if(infor(child(numCH+1,jjj)).type=='MN')
      infor(n+1).MNnum=infor(n+1).MNnum+1;
    end
    jjj=jjj+1;
    if(jjj>size(child,2))   %防止超出矩阵维度
             break;
    end
end


for i=1:n
    if(infor(i).type=='IN')
        numIN=numIN+1;
        numIN_num(1,numIN)=i;  %存入IN编号矩阵
    end
end
fprintf("IN数:%d\n",numIN)

cantreach=0;
%% 处理IN
for xx=1:n
    INI(xx).dep=0;
    INI(xx).par=n+1;
    INI(xx).reach=1;
    if(infor(xx).type =='IN' & infor(xx).par~=n+1)
       parentCH = infor(xx).par;    %获取该IN的父亲CH
       INI(xx).dep = infor(parentCH).dep+1;
       grandparentCH = infor(parentCH).par;   %获取父亲的父亲CH
       if(grandparentCH==n+1)   %如果爷爷是基站，可以直接通信
           parentIN = n+1;
       else
           parentIN = infor(grandparentCH).IN;
       end   
       INI(xx).par = parentIN;  %IN路由的父亲IN
    elseif( infor(xx).type=='IN' & infor(xx).par==n+1 )
        INI(xx).dep = 1;
        INI(xx).par = n+1;  
    end
end


INI(n+1).par = n+1;
for nn=1:numIN
    inn=numIN_num(1,nn);
    if(INI(inn).par==0)
        pa = infor(inn).par;   %父亲CH
        grand = infor(pa).par;   %爷爷CH
        parentin=INI(inn).par;
        while(parentin == 0)
            grand=infor(grand).par;  %继续向上找
            parentin=infor(grand).IN;
            if(parentin~=0)
                INI(inn).par = parentin;  %IN路由的父亲IN
                break
            end
        end
    end
end


%%找到与父亲IN不可通信的IN
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




%% 下面开始模拟产生包
 for i=1+10*(time-1):10+10*(time-1)     %一个大轮，包含20小轮/202103修改为10轮
 % fprintf("第%d轮\n",i);
  if(death==0)  
  for l=1:n        %对MN进行模拟，产生数据包
      if((infor(l).type == 'MN') & (infor(l).vis~=0) & (infor(l).malicious==0))    %如果该节点是MN，则开始模拟
        
          Receive(l,i+1)=Receive(l,i+1)+1;  %MN的收包成功数+1，其实没有收包，MN在每轮是要固定+1的
          %send=send+1;%计算吞吐量
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
%               infor(l).fratetemp=infor(l).frate;  %保存正常转发率
%               infor(l).frate = infor(l).frate*0.8;  %改变转发率
%           end
          rt=rt+1;
          if(1-infor(l).frate)<my(1,rt)   %如果丢包率小于随机数,且强化学习发包，则可以发包，即父亲可以接收到
                parent = infor(l).par;         %获取父亲
                %mm=infor(l).sort;        %获取MN的单独编号，以便存放在矩阵
                Forward(l,i+1)=Forward(l,i+1)+1;  %MN的发包成功数+1
                infor(l).forward(1,i+1)=1;
                infor(l).RE = infor(l).RE - Etx2(infor(l).dist);  %发送包，消耗能量！
                %% 处理父亲
                 while (parent~=n+1 & infor(parent).malicious==0 )
                   %yy=infor(parent).sort;      %获取父亲的簇头单独编号，以便存放在矩阵
                   Receive(parent,i+1)=Receive(parent,i+1)+1;    %第i+1列，因为目前是第i轮
                   infor(parent).receive(1,i+1)=1;
                   infor(parent).RE = infor(parent).RE - Erx(552);  %父亲接受包，消耗能量！
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
%                         infor(parent).fratetemp=infor(parent).frate;  %保存正常转发率
%                         infor(parent).frate = infor(parent).frate*0.8;  %改变转发率
%                    end
                   rt=rt+1;
                   if (1-infor(parent).frate)<my(1,rt)   %如果父亲的丢包率小于随机数,则可以转发该包
                          Forward(parent,i+1)=Forward(parent,i+1)+1;
                          infor(parent).forward(1,i+1)=1;
                          infor(parent).RE = infor(parent).RE - Etx2(infor(parent).dist);  %发送包，消耗能量！
                          parent=infor(parent).par;    %获取下一跳簇头
                   else      %如果该次转发失败，则停止继续向下一簇头转发
                         break
                   end
         
                 end
          end 
       
      
      elseif(infor(l).type=='IN' & infor(l).par~=n+1)
      %% 处理IN路由
              bb=l;
              %%首先本获取簇成员数packetnum
              parIN=INI(bb).par;    %父亲IN
              parCH=infor(bb).par;   %父亲簇首
              packetnum = infor(parCH).MNnum -1 +1;  %发包数，减去MNnum中的IN，增加一个本簇簇首1
              infor(bb).RE = infor(bb).RE - packetnum * Etx3(INI(bb).dist);           %发包能耗
              while(parIN~=n+1)     %只要没到SN,就一直传递下去
                  infor(parIN).RE = infor(parIN).RE - packetnum * Erx(80);       %父亲收包耗能
                  infor(parIN).RE = infor(parIN).RE - packetnum * Etx3(INI(parIN).dist); %发包能耗
                  parIN=INI(parIN).par;
              end      
      
       elseif(infor(l).type=='IN' & infor(l).par==n+1)   %BS的IN
              bb=l;
              packetnum = infor(parCH).MNnum ;  %发包数
              infor(bb).RE = infor(bb).RE - packetnum * Etx1(INI(bb).dist);           %发包能耗
             
       end
  end
  
  if(size(extra_for,2)~=0)
   for extra=1:size(extra_for,2)   %处理需要多发包的MN
      ee=extra_for(extra);
      counte=0;
      if((infor(ee).type == 'MN') & (infor(ee).vis~=0) & (infor(ee).malicious==0))    %如果该节点是MN，则开始模拟
          while(counte<15)
           Receive(ee,i+1)=Receive(ee,i+1)+1;  %MN的收包成功数+1，其实没有收包，MN在每轮是要固定+1的
           %send=send+1;%计算吞吐量
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
%               infor(ee).fratetemp=infor(ee).frate;  %保存正常转发率
%               infor(ee).frate = infor(ee).frate*0.8;  %改变转发率
%           end
           rt=rt+1;
           if (1-infor(ee).frate)<my(1,rt) %如果丢包率小于随机数,则可以发包，即父亲可以接收到
                parent = infor(ee).par;         %获取父亲
                Forward(ee,i+1)=Forward(ee,i+1)+1;  %MN的发包成功数+1
                infor(ee).forward(1,i+1)=1;
                infor(ee).RE = infor(ee).RE - Etx2(infor(ee).dist);  %发送包，消耗能量！
                %% 处理父亲
                 while (parent~=n+1 & infor(parent).malicious==0 )
                     Receive(parent,i+1)=Receive(parent,i+1)+1;    %第i+1列，因为目前是第i轮
                     infor(parent).receive(1,i+1)=1;
                     infor(parent).RE = infor(parent).RE - Erx(552);  %父亲接受包，消耗能量！
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
%                         infor(parent).fratetemp=infor(parent).frate;  %保存正常转发率
%                         infor(parent).frate = infor(parent).frate*0.8;  %改变转发率
%                      end
                     rt=rt+1;
                     if (1-infor(parent).frate)<my(1,rt)  %如果父亲的丢包率小于随机数,则可以转发该包
                          Forward(parent,i+1)=Forward(parent,i+1)+1;
                          infor(parent).forward(1,i+1)=1;
                          infor(parent).RE = infor(parent).RE - Etx2(infor(parent).dist);  %发送包，消耗能量！
                     end
                     parent=n+1;    %发不发完都强制结束
                 end
           end 
          end
      end
   end
  end
%   %% 处理需要多发包的MN
  if(size(judge_now,2)~=0 & size(judge_now,1)~=0)
   for xxx=1:judge_num
     obj3=judge_now(xxx);
     if(infor(obj3).type=='MN' & (infor(obj3).vis~=0) & (infor(obj3).malicious==0) )  %多发包，ttl=1
%              obj3;
             packnum2=0;
             while(packnum2<15)
               Receive(obj3,i+1)=Receive(obj3,i+1)+1;  %MN的收包成功数+1，其实没有收包，MN在每轮是要固定+1的
               %send=send+1;%计算吞吐量
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
%                     infor(obj3).fratetemp=infor(obj3).frate;  %保存正常转发率
%                     infor(obj3).frate = infor(obj3).frate*0.8;  %改变转发率
%                end
               rt=rt+1;
                if (1-infor(obj3).frate)<my(1,rt)  %如果丢包率小于随机数,则可以发包，即父亲可以接收到
                  Forward(obj3,i+1)=Forward(obj3,i+1)+1;  %MN的发包成功数+1
                  infor(obj3).forward(1,i+1)=1;
                  infor(obj3).RE = infor(obj3).RE - Etx2(infor(obj3).dist);  %发送包，消耗能量！
                end
             end
     end
   end
  end
%一轮结束，把该轮对应的整列复制一遍给下一列
 if(i<k)
  copy_reci = Receive(:,i+1);
  Receive(:,i+2) = copy_reci;
  copy_for = Forward(:,i+1);
  Forward(:,i+2) = copy_for;
 end
 
 %计算当轮cfr
    for io=1:n
            now_for = Forward(io,i+1);
            now_reci = Receive(io,i+1);
            cfr = now_for/now_reci; 
         if(i<=10) %20210527由20修改为10  
            if (Receive(io,i+1)-Receive(io,i)==0)  %没有收包
               cfrtemp = CFRtemp(io,i);
            else
               cfrtemp = (Forward(io,i+1)-Forward(io,i))/(Receive(io,i+1)-Receive(io,i));
            end
         else
             symbol=1+10*(time-1);%202103修改每10轮为一大轮
             if (Receive(io,i+1)-Receive(io,symbol)==0)  %没有收包
               cfrtemp = CFRtemp(io,i);
             else
               cfrtemp = (Forward(io,i+1)-Forward(io,symbol))/(Receive(io,i+1)-Receive(io,symbol));
             end
         end
            CFR(io,i+1)=cfr;
            CFRtemp(io,i+1)=cfrtemp;   %当轮转发率保存
    end
 
 
 
 if (i==200||i==300||i==400||i==500||i==600||i==700||i==800||i==900||i==1000||i==1100||i==1200)  %%200轮开始，SN每轮都进行恶意节点检测
      %获取簇首数据
      %now_round=[];
      %now_round(:,1) = CFRtemp(:,i+1);         %第一列当轮单轮
%       now_round=[];
%       now_round(:,1) = CFR(:,i);         %第一列当轮单轮
%       now_round(:,2) = CFR(:,i+1);  %第二列当轮累计
      now_round=zeros(n+1,2);
      for ioo=1:n
          cfrw=CFRtemp(ioo,i-198:i+1);
          cfrh=CFRtemp(ioo,i-198:i+1);
%           now_round(ioo,1) = exp(-std(cfrh,1));        %第一列当轮单轮
          if mean(cfrw(:))~=0&infor(ioo).type ~= 'DE'
               now_round(ioo,1) = std(cfrh)/mean(cfrw(:));
               now_round(ioo,2)=CFR(ioo,i+1);  %第二列当轮累计
          elseif infor(ioo).type == 'DE'
              now_round(ioo,1)=10;
              now_round(ioo,2)=0;  %第二列当轮累计
          else
              now_round(ioo,1)=-1;
              now_round(ioo,2)=CFR(ioo,i+1);  %第二列当轮累计
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
     % fid=fopen(['forrho','.txt'],'w');%写入文件路径  
      fid=fopen(['202010182','.txt'],'w');%写入文件路径  
      [r1,c]=size(now_round);% 得到矩阵的行数和列数   
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
%    % title ('加入强化学习节点rho','FontSize',15.0)
%     xlabel ('\rho')
%     ylabel ('\delta')
    suspect_nowcopy=suspect_now;  %复制一份
    judge_nowcopy=judge_now;  %复制一份
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

    %% judge节点，若是CH，多收包；若是MN，多发包 ；IN非恶意、不考虑
    for kk=1:size(judge_nowcopy,2)
        if(infor(judge_nowcopy(kk)).type=='IN')
          judge_now(judge_now==judge_nowcopy(kk))=[];   %若IN，不算
          judge_num=judge_num-1;
        elseif(infor(judge_nowcopy(kk)).malicious~=0)
          judge_now(judge_now==judge_nowcopy(kk))=[];   %若死亡，不算
          judge_num=judge_num-1;
        end
    end
    %judge节点存入矩阵
    for kk2=1:size(judge_now,2)
        if(~ismember(judge_now(kk2),ever_judge))
            ever_judge=[ever_judge judge_now(kk2)];
        end 
    end
    
    %% 观察，需要judge的节点是否rho高于过阈值
    if(size(rho_everjudge,2)~=0)
     for  j4=1:size(rho_everjudge,2)
        if(rho_everjudge(j4)>whzhigh & ~ismember(ever_judge(j4),release_judge))
           release_judge=[release_judge ever_judge(j4)];   %保存一个解除嫌疑矩阵
        end
     end
    end

    %% 把judge中release_judge的节点删去
    judge_now
    judge_now = setdiff(judge_now,release_judge);
    judge_now
    judge_num = size(judge_now,2);
    judge_nowcopy2=judge_now;  %复制一份
   %%%%%%%至此，当轮judge节点固定%%%%%%
   
   %%%%%%%下面确定当轮suspect节点%%%%%% 
  %筛除IN节点
  for nnn=1:size(suspect_nowcopy,2)       
      isub=suspect_nowcopy(nnn);
      if(infor(isub).type=='IN') %如果是IN，不算
          suspect_now(suspect_now==isub)=[];
          suspect_num=suspect_num-1;
          fprintf('%d号节点被解除嫌疑,是IN\n',isub);
      elseif(Receive(isub,i)<50) %如果长时间作为IN，几乎不参与环境感知
          suspect_now(suspect_now==isub)=[];
          suspect_num=suspect_num-1;
          fprintf('%d号节点被解除嫌疑\n',isub);
      end
  end
  suspect_nowcopy2=suspect_now;  %再复制一份
  %排除正常节点的嫌疑
%   for iii=1:size(suspect_nowcopy2,2)
%       subj=suspect_nowcopy2(iii);
%       if now_round(subj,1)/now_round(subj,2)<0.8
%            suspect_now(suspect_now==subj)=[];   %解除嫌疑
%            suspect_num=suspect_num-1;
%            if(infor(subj).malicious==0)
%              fprintf('%d号节点被解除怀疑，可能是环境因素\n',subj);
%            end
%       end
%   end
%           
  %处理suspect 邻居监测，如果与邻居转发率很接近，说明是harsh或者是0.7+的恶意节点，可以解除嫌疑；如果很接近但被漏检，不算作MDR
%     for iii=1:size(suspect_nowcopy2,2)
%         FR_ne=[];
%         subj=suspect_nowcopy2(iii); 
%         subj_ne=nighbor_no(subj,:);   %subj_ne存放该恶意节点的邻居编号矩阵
%         for j=1:size(subj_ne,2)
%             if(subj_ne(j)~=0)
%               FR_ne=[FR_ne now_round(subj_ne(j),2)]; 
%             else
%                 break
%             end
%         end
%         %对FR_ne矩阵求均值、标准差
%         FRmean=mean(FR_ne);
%         FRstd=std(FR_ne);
%         subj_FR=now_round(subj,2);   %待确认恶意节点转发率（判断对象）
%         if(subj_FR-FRmean>0||FRmean-subj_FR<FRstd)
%             suspect_now(suspect_now==subj)=[];   %解除嫌疑
%             suspect_num=suspect_num-1;
%             if(infor(subj).malicious==0)
%               fprintf('%d号节点被解除怀疑，可能是环境因素\n',subj);
%             end
%         end
%     end
%     
  %采用邻居投票的方式区别恶劣环境中的正常节点和高转发率的恶意节点
   for iii=1:size(suspect_nowcopy2,2)
        FR_ne=[];
        subj=suspect_nowcopy2(iii); 
        subj_ne=nighbor_no(subj,:);   %subj_ne存放该恶意节点的邻居编号矩阵
        for j=1:size(subj_ne,2)
            if(subj_ne(j)~=0)
              FR_ne=[FR_ne (sin(5/8*pi*now_round(subj_ne(j),1)+0.5*pi)+1)/now_round(subj_ne(j),2)]; 
            else
                break
            end
        end
        
        %对FR_ne矩阵求均值、标准差
%         FRmean=mean(FR_ne);
%         FRstd=std(FR_ne);
        subj_FR=(sin(5/8*pi*now_round(subj,1)+0.5*pi)+1)/now_round(subj,2);
        appr0val=0;
        opp0se=0;
        for vote=1:size(FR_ne,2)%遍历所有邻居，邻居CFR<0.8，则支持，反之则反对
            if abs(FR_ne(vote)-subj_FR)<=0.27
                appr0val=appr0val+1;
            else
                opp0se=opp0se+1;
            end
        end
        appr0val_rate=appr0val/(appr0val+opp0se);%计算支持率
        %支持率大于阈值，则该节点和周围环境相似。阈值计算：邻域与恶劣环境重叠面/邻域面积
        %当恶劣环境面积和邻域面积相等时，阈值=0.391
        if(appr0val_rate>0.391)
            suspect_now(suspect_now==subj)=[];   %解除嫌疑
            suspect_num=suspect_num-1;
            if(infor(subj).malicious==0)
              fprintf('%d号节点被解除怀疑，可能是环境因素\n',subj);
            end
        end
    end
  %%%%%至此，本轮suspect全部确定%%%%%%
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
%            suspect_now(suspect_now==subjj)=[];   %解除嫌疑
%            suspect_num=suspect_num-1;
%            if(infor(subjj).malicious==0)
%              fprintf('%d号节点被解除怀疑，可能是环境因素\n',subjj);
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
%            suspect_now(suspect_now==subjj)=[];   %解除嫌疑
%            suspect_num=suspect_num-1;
%            if(infor(subjj).malicious==0)
%              fprintf('%d号节点被解除怀疑，可能是环境因素\n',subjj);
%            end
%        end
%    end
%     for qq=1:nm
%         neighbor_nums(qq,2)=sum(neighbor_matrix(qq,:));
%     end
%     for hz=1:5%取五次最大值
%         xxhz=min(min(DC));
%         [x,y]=find(DC==xxhz);%求出最大值的位置
%         for ziii=1:size(x,1)%找到具有最大邻居数量的节点
%           zii=x(ziii);
%           subjj=suspect_nowcopy2(zii);
%           suspect_now(suspect_now==subjj)=[];   %解除嫌疑
%           suspect_num=suspect_num-1;
%           if(infor(subjj).malicious==0)
%             fprintf('%d号节点被解除怀疑，可能是环境因素\n',subjj);
%           end
%           neighbor_all=[neighbor_all,nighbor_no(subjj,:)];%邻居之和
%           DC(zii,2)=max(max(DC));
%         end
%     end
%         
%         
%      for zi=1:nm
%          subjjj=suspect_nowcopy2(zi);
%         if ismember(subjjj, neighbor_all)>0
%            suspect_now(suspect_now==subjjj)=[];   %解除嫌疑
%            suspect_num=suspect_num-1;
%            if(infor(subjjj).malicious==0)
%              fprintf('%d号节点被解除怀疑，可能是环境因素\n',subjjj);
%            end
%         end
%      end   
  
  
  %% judge suspect已经筛选一遍，记录最终信息
  %把非judge节点的increase置0
  for nn2=1:n       
      if(ismember(nn2,judge_now) == 0 & infor(nn2).malicious == 0) %如果不在judge节点中，则置0
          infor(nn2).increase=0;
      end
  end
  %把此时judge节点increase置1,conssuspect+1
 if(size(judge_now,2)~=0 & size(judge_now,1)~=0)
  for kk=1:size(judge_now,2)  
      if(infor(judge_now(kk)).malicious==0)
         infor(judge_now(kk)).increase=1;
         infor(judge_now(kk)).conssuspect=infor(judge_now(kk)).conssuspect+1;
      end
  end
 end
  %记录需要向父亲CH多发包的节点
    extra_for=[];
    for lll=1:numCH
        cc=numCH_num(lll);
      if(infor(cc).increase==1)  %对judgeCH的孩子，让他们多发包，在此记录
           fff=infor(cc).sort;
           for ppp=1:numCH
               if(child(ppp,1)==fff)
                  break
               end
           end
           lie=2;
           while(child(ppp,lie)~=0)   %记录需要额外发包的MN
               extra_for=[extra_for child(ppp,lie)];
               lie=lie+1;
                 if(lie>size(child,2))   %防止超出矩阵维度
                   break
                 end
           end
      end
    end
  %处理suspect类节点
  for nn=1:n       %对所有节点遍历，把非怀疑节点的successive置0
      if(ismember(nn,suspect_now) == 0 & infor(nn).malicious == 0) %如果不在怀疑节点中，则置0
          infor(nn).successive=0;
      end
  end
  
  for yy=1:suspect_num     %把该轮怀疑节点标注为1
      zz=suspect_now(1,yy);
      if(infor(zz).malicious == 0)      
         infor(zz).successive = infor(zz).successive+1;
      end
  end
  
   for mm=1:n      %观察是否有某节点连续五轮被怀疑
      if(infor(mm).successive == 1 & infor(mm).malicious==0 & infor(mm).vis~=0)    %如果该节点连续五轮被怀疑(并且还未被确认)，则确认恶意
          infor(mm).malicious = i;    %记录是在第几轮被检测出来的！
          infor(mm).type = 'DE';
          infor(mm).locationx = S(mm).x;
          infor(mm).locationy = S(mm).y;
          S(mm).x=500;
          S(mm).y=500;
          confirm=[confirm mm];    %把该确认点加到confirm矩阵里
       elseif(infor(mm).conssuspect == 25 & infor(mm).malicious==0 & infor(mm).vis~=0)
           infor(mm).malicious = i;    %记录是在第几轮被检测出来的！
           infor(mm).type = 'DS';
           infor(mm).locationx = S(mm).x;
           infor(mm).locationy = S(mm).y;
           S(mm).x=500;
           S(mm).y=500;
           confirm=[confirm mm];    %把该确认点加到confirm矩阵里
      end
      
      
   end
  miss=setdiff(MCH,confirm);    %怀疑节点中怀疑错的！也就是漏检
  fault=setdiff(confirm,MCH) ;    %怀疑节点没找出来的！也就是误检
  misscopy=miss; %复制一份  
   if size(fault,2)~=0
   fprintf("误检节点：")
   fault
   FDR = size(fault,2)/(n-mali_num);   %计算FDR
  else
   FDR = 0;   %计算FDR
  end
  if size(miss,2)~=0
   fprintf("漏检节点：")
   miss
   MDR = size(miss,2)/mali_num;   %计算MDR
  else
   MDR = 0;   %计算FDR
  end
  FDRM=[FDRM,FDR];
  MDRM=[MDRM,MDR];
  end
  
 %% 每轮检查是否有一个节点死亡，如果有，death=1
 for vv=1:n
     % 记录所有节点能量
     %energy(vv,i)=infor(vv).RE;
     if(infor(vv).RE<=0)
         vv
         death=1;
         fprintf('$$$$$$$$$$$$$$')
         deathtime=[deathtime i];
         if(i>225)   %启动了检测
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
if(flag==1)  %此次仿真可用
  simulation=simulation+1;
  confirm
  miss=setdiff(MCH,confirm);    %怀疑节点中怀疑错的！也就是漏检
  fault=setdiff(confirm,MCH) ;    %怀疑节点没找出来的！也就是误检
  misscopy=miss; %复制一份
%   fid10=fopen(['FDR and MDR','.txt'],'w');
%   for www=1:29
%       fprintf(fid7,'%f',FDRM(www));
%       fprintf(fid7,'%s',' ');
%       fprintf(fid7,'%f',MDRM(www));
%       fprintf(fid7,'\n');
%   end
%   fclose(fid10);
  %% 对漏检节点筛选
%   for i=1:size(misscopy,2)
%         FRm_ne=[];
%         subjm=misscopy(i); 
%         subjm_ne=nighbor_no(subjm,:);   %subj_ne存放该恶意节点的邻居编号矩阵
%         for j=1:size(subjm_ne,2)
%             if(subjm_ne(j)~=0)
%               FRm_ne=[FRm_ne now_round(subjm_ne(j),2)]; 
%             else
%                 break
%             end
%         end
%         %对FRm_ne矩阵求均值、标准差
%         FRmm=mean(FRm_ne);
%         FRstdm=std(FRm_ne);
%         subjm_FR=now_round(subjm,2);   
%         if(abs(subjm_FR-FRmm)<2*FRstdm)
%             miss(miss==subjm)=[];   %不算作MDR
%             fprintf('%d号节点不算MDR\n',subjm);
%         end
%    end
  
  
  
  if size(fault,2)~=0
   fprintf("误检节点：")
   fault
   FDR = size(fault,2)/(n-mali_num);   %计算FDR
  else
   FDR = 0;   %计算FDR
  end
  if size(miss,2)~=0
   fprintf("漏检节点：")
   miss
   MDR = size(miss,2)/mali_num;   %计算MDR
  else
   MDR = 0;   %计算FDR
  end
  Accuracy=(n-size(fault,2)-size(miss,2))/n;%计算检测正确率
  %计算throughput
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
   zongjie(simulation,1)=FDR;%第一列FDR
   zongjie(simulation,2)=MDR; %第二列MDR
%   zongjie(simulation,3)=size(fault,2);%第三列误检节点个数
%   zongjie(simulation,4)=numCH;%第四列簇头个数
%   zongjie(simulation,5)=size(miss,2);%第五列漏检节点个数
  zongjie(simulation,3)=Accuracy;%第六列计算检测正确率
%   %第七列第八列网络吞吐量
%   %zongjie(simulation,7)=throughput;
   zongjie(simulation,4)=throughput_rate;
  
end
  
end
end
 
 zongjie
 %fprintf("结果均值：第一列fdr 第二列mdr：\n");
% whzzongjie_mean(whzlow-81,:)=100*mean(zongjie,1);
whzzongjie_mean(whzlow-49,whzhigh-235:whzhigh-232)=100*mean(zongjie,1);
%whzzongjie(whzlow-81,whzhigh-177:whzhigh-176)=zongjie_simu;
wxc(whzlow-49).zongjie=zongjie_simu;
 %fprintf("结果方差：第一列fdr 第二列mdr：\n");
 % whzzongjie_var(whzlow-81,:)=var(zongjie);
 whzzongjie_var(whzlow-49,whzhigh-235:whzhigh-232)=var(zongjie);
 %whzzongjie_var(whzlow-81,whzhigh-171:whzhigh-165)=var(zongjie);
    end
end

  fprintf("结果均值：第一列FDR，第二列MDR，第三列正确率，第四列：throughput_rate：\n");
  whzzongjie_mean
  fprintf("结果方差：第一列FDR，第二列MDR，第三列正确率，第四列：throughput_rate：\n");
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
function [condition,m,n,q,p,r] = DBSCAN_NEW(filename,condition2,clusterhead,a,b,times)
xx=load(filename,'s');
ND=max(xx(:,2));
NL=max(xx(:,1));
if (NL>ND)
  ND=NL; %% 确保 DN 取为第一二列最大值中的较大者，并将其作为数据点总数
end
N=size(xx,1); %% xx 第一个维度的长度，相当于文件的行数（即距离的总个数）
dist=zeros(ND,ND);%节点距离矩阵
time=zeros(ND,1);%节点深度
judge=[];   %可疑节点矩阵
judge_num=0;   %可疑节点数量
suspect_MCH=[];  %怀疑节点矩阵，负责存放怀疑节点编号
suspect_num=0;    %由算法产生的怀疑节点数量值
whz=size(clusterhead,1);
% Neigh=zeros(ND,1);
% ep=0.10;
%% 利用 xx 为 dist 数组赋值，注意输入只存了 0.5*DN(DN-1) 个值，这里将其补成了满矩阵
%% 这里不考虑对角线元素
for i=1:N
    ii=xx(i,1);
    jj=xx(i,2);
    dist(ii,jj)=xx(i,3);
    dist(jj,ii)=xx(i,3);
end
DIST=zeros(whz,2);
if whz>1
   for i=1:whz
       DIST(i,1)=clusterhead(i);
   end
   for ii=1:(whz-1)
       for j=ii+1:whz
           DIST(ii,2)=DIST(ii,2)+dist(clusterhead(ii),clusterhead(j));
           DIST(j,2)=DIST(j,2)+dist(clusterhead(ii),clusterhead(j));
       end
   end
   DIST1=sortrows(DIST,2);
   depth=DIST1(1,1);
else
   depth=clusterhead;
end
percent=2.0;
% for i=1:ND-1
%     for j=i+1:ND
%         if dist(i,j)<=ep
%             Neigh(i,1)=Neigh(i,1)+1;
%             Neigh(j,1)=Neigh(j,1)+1;
%         end
%     end
% end
 
position=round(N*percent/100);   %% round 是一个四舍五入函数
sda=sort(xx(:,3));   %% 对所有距离值作升序排列
dc=sda(position);
rho=zeros(ND,2);
% for i=1:ND
%   rho(i,1)=i;
%   rho(i,2)=0;
% end
% for i=1:ND-1
%   for j=i+1:ND
%      rho(i,2)=rho(i,2)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
%      rho(j,2)=rho(j,2)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
%   end
% end
% ordrho=sortrows(rho,-2);
% depth=ordrho(1,1);
%% 确定 epsilon
% if times<41
%     epsilon=0.20;
% elseif times<81
%     epsilon=0.17;
% elseif times<121
%     epsilon=0.14;
% elseif times<161
%     epsilon=0.14;
% elseif times<201  
%     epsilon=0.12;
% else
    epsilon=0.13;
% end

%% 开始计算邻居节点
hz1=0;
for i=1:15
    hz=size(depth,2);
    for ii=hz1+1:hz
        clusterhead1=depth(ii);
        for iii=1:ND
            if dist(clusterhead1,iii)<epsilon&&time(iii)==0
                time(iii)=i;
                depth=[depth,iii];
            end
        end
    end
    hz1=hz;
    epsilon=epsilon*0.8;
end
for i1=1:ND
%     if time(i1)>4&&time(i1)<8
%        suspect_num = suspect_num+1;
%        suspect_MCH(1,suspect_num) = i1;
%     end
    if time(i1)==0||time(i1)>15
        suspect_num = suspect_num+1;
        suspect_MCH(1,suspect_num) = i1;
    end
end
condition=condition2;
m=suspect_MCH;
n=suspect_num;
q=judge;
p=judge_num;
r=[];
fid8=fopen(['20220820','.txt'],'w');
for i=1:ND
    fprintf(fid8,'%d',i);
    fprintf(fid8,'%s',' ');
    fprintf(fid8,'%f',time(i));
    fprintf(fid8,'\n');
end
fclose(fid8);
end














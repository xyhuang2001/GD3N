 xx=load('20220809.txt','s');
 yy=load('20220812.txt','s');
 zz=load('20220820.txt','s');
 ww=load('20230602.txt','s');
 max(zz)
 find(zz(:,2)>10);
 mate=zeros(500,1);
 pro=zeros(500,1);
 for ii4=1:75
     for ii5=1:500
     if yy(ii4)==ii5
        mate(ii5)=100;
     end
     end
 end
 for i=1:size(ww(:,1),1)
     a=ww(i,1);
     pro(a,1)=100;
 end
 for ii6=1:500
     if mate(ii6)==0
         plot(xx(ii6,2),xx(ii6,3),'o','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','g');hold on
     end
 end
  for ii7=1:500
     if mate(ii7)==100
         plot(xx(ii7,2),xx(ii7,3),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on
     end
  end
 xlabel('Std/Mean');
 ylabel('CFR');
 title('Suspicious Data Point Distribution');
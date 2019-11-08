%空间二次曲面拟合算法
num_points = length(x);
x=x*100/65535;
y=y*100/65535;
%z=z*100/65535;


 x_max=max(x);
 y_max=max(y);
 x_min=min(x);
 y_min=min(y);
 %z_max=max(z);
 %z_min=min(z);
 
 x=x-(x_max+x_min)/2;
 y=y-(y_max+y_min)/2;
 %z=z-(z_max+z_min)/2;
 
 
figure;
plot(x,y,'*');
axis equal
 

 x_avr = sum(x)/num_points;
 y_avr = sum(y)/num_points;
 %z_avr = sum(z)/num_points;
 
 xx_avr = sum(x.*x)/num_points;
 yy_avr = sum(y.*y)/num_points;
 %zz_avr = sum(z.*z)/num_points;
 xy_avr = sum(x.*y)/num_points;
 %xz_avr = sum(x.*z)/num_points;
 %yz_avr = sum(y.*z)/num_points;
 
 xxx_avr = sum(x.*x.*x)/num_points;
 xxy_avr = sum(x.*x.*y)/num_points;
  %xxz_avr = sum(x.*x.*z)/num_points;
 xyy_avr = sum(x.*y.*y)/num_points;
  %xzz_avr = sum(x.*z.*z)/num_points;
 yyy_avr = sum(y.*y.*y)/num_points;
  %yyz_avr = sum(y.*y.*z)/num_points;
  %yzz_avr = sum(y.*z.*z)/num_points;
  %zzz_avr = sum(z.*z.*z)/num_points;
 
  %xyz_avr = sum(z.*z.*z)/num_points;
 
 yyyy_avr = sum(y.*y.*y.*y)/num_points;
  %zzzz_avr = sum(z.*z.*z.*z)/num_points;
 xxyy_avr = sum(x.*x.*y.*y)/num_points;
  %xxzz_avr = sum(x.*x.*z.*z)/num_points;
  %yyzz_avr = sum(y.*y.*z.*z)/num_points;
 
 
 xxxx_avr = sum(x.*x.*x.*x)/num_points;
 
 xxxy_avr = sum(x.*x.*x.*y)/num_points;
  %xxxz_avr = sum(x.*x.*x.*z)/num_points; 
  %xxyz_avr = sum(x.*x.*y.*z)/num_points;
 xyyy_avr = sum(x.*y.*y.*y)/num_points;
  %xyyz_avr = sum(x.*y.*y.*z)/num_points;
  %xyzz_avr = sum(x.*y.*z.*z)/num_points;
 % xzzz_avr = sum(x.*z.*z.*z)/num_points;
 
 yyxy_avr = sum(y.*y.*x.*y)/num_points;
 % yyxz_avr = sum(y.*y.*x.*z)/num_points;
  %yyyz_avr = sum(y.*y.*y.*z)/num_points;
  %yzzz_avr = sum(y.*z.*z.*z)/num_points;
 
 % zzxy_avr = sum(z.*z.*x.*y)/num_points;
 % zzxz_avr = sum(z.*z.*x.*z)/num_points;
 % zzyz_avr = sum(z.*z.*y.*z)/num_points;
%----------------带约束条件的最小二乘。。结果是天文数字。。-------------
%syms A B C D E F M
%[A,B,C,D,E,F]=solve([(A*xxxx_avr+B*xxyy_avr+C*xxxy_avr+D*xxx_avr+E*xxy_avr+F*xx_avr)],[(A*xxyy_avr+B*yyyy_avr+C*xyyy_avr+D*xyy_avr+E*yyy_avr+F*yy_avr)],[(A*xxxy_avr+B*xyyy_avr+C*xxyy_avr+D*xxy_avr+E*xyy_avr+F*xy_avr)],[(A*xxx_avr+B*xyy_avr+C*xxy_avr+D*xx_avr+E*xy_avr+F*x_avr)],[(A*xxy_avr+B*yyy_avr+C*xyy_avr+D*xy_avr+E*yy_avr+F*y_avr)],[(A*xx_avr+B*yy_avr+C*xy_avr+D*x_avr+E*y_avr+F)])
%-----------解算一般椭圆方程（带角度旋转）首位系数为1------------------------------------
% N1=F1*[x';y'];
%  x_2=N1(1,:);
%  y_2=N1(2,:);
%   figure;
%  plot(x_2,y_2,'*');
%-----------解算一般椭圆方程（带角度旋转）首位系数为a------------------------------------ 
% C0=[xxxx_avr xxyy_avr xxzz_avr xxxy_avr xxxz_avr xxyz_avr xxx_avr xxy_avr xxz_avr;
%     xxyy_avr yyyy_avr yyzz_avr xyyy_avr xyyz_avr yyyz_avr xyy_avr yyy_avr yyz_avr;
%     xxzz_avr yyzz_avr zzzz_avr xyzz_avr xzzz_avr yzzz_avr xzz_avr yzz_avr zzz_avr;
%     xxxy_avr xyyy_avr xyzz_avr xxyy_avr xxyz_avr xyyz_avr xxy_avr xyy_avr xyz_avr;
%     xxxz_avr xyyz_avr xzzz_avr xxyz_avr xxzz_avr xyzz_avr xxz_avr xyz_avr xzz_avr;
%     xyyz_avr yyyz_avr yzzz_avr xyyz_avr xyzz_avr yyzz_avr xyz_avr yyz_avr yzz_avr;
%     xxx_avr xyy_avr xzz_avr xxy_avr xxz_avr xyz_avr xx_avr xy_avr xz_avr;
%     xxy_avr yyy_avr yzz_avr xyy_avr xyz_avr yyz_avr xy_avr yy_avr yz_avr;
%     xxz_avr yyz_avr zzz_avr xyz_avr xzz_avr yzz_avr xz_avr yz_avr zz_avr 
%     ];
% B0=[20000*xx_avr;
%     20000*yy_avr;
%     20000*zz_avr;
%     20000*xy_avr;
%     20000*xz_avr;
%     20000*yz_avr;
%     20000*x_avr;
%     20000*y_avr;
%     20000*z_avr;
%     ];
   C0=[xxxx_avr xxyy_avr xxxy_avr xxx_avr xxy_avr 
       xxyy_avr yyyy_avr xyyy_avr xyy_avr yyy_avr 
       xxxy_avr xyyy_avr xxyy_avr xxy_avr xyy_avr 
       xxx_avr  xyy_avr  xxy_avr  xx_avr  xy_avr  
       xxy_avr  yyy_avr  xyy_avr  xy_avr  yy_avr  
       
       ];
   B0=[-xx_avr;
       -yy_avr;
       -xy_avr;
       -x_avr;
       -y_avr;
       ];
A0=inv(C0)*B0;
%resoult = solution_equations_n_yuan(A0,b);
 
% x00 = -resoult(3)/2;                 
% y00 = -resoult(4)/(2*resoult(1));     
% z00 = -resoult(5)/(2*resoult(2));     
 
% AA = sqrt(x00*x00 + resoult(1)*y00*y00 + resoult(2)*z00*z00 - resoult(6)); 
% BB = AA/sqrt(resoult(1));                                                    
% CC = AA/sqrt(resoult(2));                                                

%*****************************************************软磁矩阵计算
% A3=[resoult2(1)    resoult2(4)/2  resoult2(5)/2;
%     resoult2(4)/2  resoult2(2)    resoult2(6)/2;
%     resoult2(5)/2  resoult2(6)/2  0
%     ];
% fprintf('拟合结果\n');
% fprintf('A = %f, B = %f, C = %f\n',resoult);
% F=inv(Q*D2);
% N1=F*[x'; y'; z'];
% x1=N1(1,:);
% y1=N1(2,:);
% z1=N1(3,:);

fun = @(a)a(1).*x.*x+a(2).*y.*y+a(3).*x.*y+a(4).*x+a(5).*y+a(6);
a0 = [1 1 1 1 1 1];
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
sym1 = lsqnonlin(fun,a0,[],[],options)
L1=[0.827  0;
    0      1];
L2=[0.9869    -0.1611;
    -0.1611    0.9869];
L3=inv(L1*L2);
 N1=inv(L1*L2)*[x'; y'];
 x2=N1(1,:);
 y2=N1(2,:);
 figure;
plot(x2,y2,'*');
axis equal
 



%sym2 = lsqcurvefit(fun,a0,x,y)
    



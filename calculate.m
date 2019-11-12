%set the format of data 
num_points = length(x);
x=x*100/65535;
y=y*100/65535;

%calculate the max&min value of x y
x_max=max(x);
y_max=max(y);
x_min=min(x);
y_min=min(y);

%compensate the hard iron distortion
x=x-(x_max+x_min)/2;
y=y-(y_max+y_min)/2;

%plot image   
figure;
plot(x,y,'*');
axis equal
 
%calculate the matrix parameter
x_avr = sum(x)/num_points;
y_avr = sum(y)/num_points;
xx_avr = sum(x.*x)/num_points;
yy_avr = sum(y.*y)/num_points;
xy_avr = sum(x.*y)/num_points;
xxx_avr = sum(x.*x.*x)/num_points;
xxy_avr = sum(x.*x.*y)/num_points;
xyy_avr = sum(x.*y.*y)/num_points;
yyy_avr = sum(y.*y.*y)/num_points;
yyyy_avr = sum(y.*y.*y.*y)/num_points;
xxyy_avr = sum(x.*x.*y.*y)/num_points;
xxxx_avr = sum(x.*x.*x.*x)/num_points;
xxxy_avr = sum(x.*x.*x.*y)/num_points;
xyyy_avr = sum(x.*y.*y.*y)/num_points;
yyxy_avr = sum(y.*y.*x.*y)/num_points;

%parameter matrixs
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
%calculate the parameter of ellipse 
A0=inv(C0)*B0;

%calculate the distortion matrix 

a1=0.5*asin(-0.5*A0(3)/(A0(1)*A0(2))^0.5);
Ky=(A0(1)/A0(2))^0.5;

L1=[1, 0
    0,Ky];
L2=[cos(a1),sin(a1)
    sin(a1),cos(a1)];
L3=inv(L1*L2);
N1=inv(L1*L2)*[x'; y'];
x2=N1(1,:);
y2=N1(2,:);
figure;
plot(x2,y2,'*');
axis equal
 



%sym2 = lsqcurvefit(fun,a0,x,y)
    



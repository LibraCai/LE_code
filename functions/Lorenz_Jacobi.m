function Jac=Lorenz_Jacobi(X)

SIGMA = 10;R = 28;BETA = 8/3; %请务必检查Lorenz_Jacobi.m 以及Lorenz.m 这两个文件的三个参数的设置是否一致
x=X(1); y=X(2); z=X(3);

%f = [SIGMA*(y-x);
%    -x*z+R*x-y  ;
%    x*y-BETA*z  ];

Jac = [-SIGMA SIGMA     0
        R - z    -1    -x
            y     x -BETA];

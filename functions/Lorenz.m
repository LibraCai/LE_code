function f=Lorenz(X)

SIGMA = 10;R = 56;BETA = 8/3;
x=X(1); y=X(2); z=X(3);

f = [SIGMA*(y-x);
    -x*z+R*x-y  ;
    x*y-BETA*z  ];

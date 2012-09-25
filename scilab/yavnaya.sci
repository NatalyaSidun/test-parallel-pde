funcprot(0);
global u x t start finish acc K L N T a; 
//Правая часть дифференциального уравнения.
function y=f(x,t)
    y=sin(x*t)
endfunction
//Начальное условие
function y=fi(x)
    y=exp(0.15*x)
endfunction
//Условие на левой границе
function y=myu(t)
    y=1

endfunction
//Условие на правой границе
function y=nyu(x)
    y=2.117
endfunction
function [u,x,t]=parabol(N,K,L,T,a)
    h=L/N;
    delta=T/K;
    for i=1:N+1
        x(i)=(i-1)*h;
        u(i,1)=fi(x(i));
    end

    for j=1:K+1
        t(j)=(j-1)*delta;
        u(1,j)=myu(t(j));
        u(N,j)=nyu(t(j));
    end
    gam=a^2*delta/h^2;
    
   
    
    for j=1:K
        for i=2:N
            u(i,j+1)=gam*u(i-1,j)+(1-2*gam)*u(i,j)+gam*u(i+1,j)+delta*...
            f(x(i),t(j));
        end
    end
    

endfunction

function [u,x,t]=parabol_parallel(u)
    h=L/N;
    delta=T/K;
    
    for i=1:N+1
        x(i)=(i-1)*h;
        u(i)=fi(x(i));
    end
    
    for j=1:K+1
         for i=1:N+1
            t(j)=(j-1)*delta;
            u(i*N + j)=myu(t(j));
            u(i*N+ N +j)=nyu(t(j));
         end
     end
     gam=a^2*delta/h^2;
     for j=1:K
        for i=2:N
            u(i*N+(j+1))=gam*u((i-1)* N + j)+(1-2*gam)*u(i*N+j)+gam*u((i+1)*N + j)+delta*f(x(i),t(j));
        end
    end
endfunction
Nstart = 100;
Nfinish = 1000;
Nh = 100;
Tstart = 100;
Tfinish = 1000;
Th = 100;

L = 1000;
K = 1000;

a = 1;
u = zeros(K * N);
write(%io(2),'Linear');
//for N = Nstart:Nh:Nfinish
 //   for T = Tstart:Th:Tfinish
  //      mprintf('N = %d\n' , N);
   //     mprintf('T = %d\n' , T);
    //    tic();
     //       [u,x,t]=parabol(N,K,L,T,a);
     //   mprintf('%.3f \n', toc());
    //end
//end
//write(%io(2),'Parallel');
//for N = Nstart:Nh:Nfinish
 //   for T = Tstart:Th:Tfinish
   //     mprintf('N = %d\n' , N);
    //    mprintf('T = %d\n' , T);

    N = 100;
    T = 100;
      tic();
           [u,x,t]=parallel_run(u, "parabol_parallel");
      mprintf('%.3f \n', toc());
    //end
//end

//[u,x,t]=parabol(N,K,L,T,a);

//mprintf('%.3f \n', finish);

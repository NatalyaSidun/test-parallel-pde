funcprot(0);
global Nstart Nfinish Nh  Tstart  Tfinish Th linear duration;

function  MatrixVector(N,T,)
        A = rand(N,N);
        B = rand(T);
        C = A * B;
       
    
endfunction

function main()
    Nstart = 100;
    Nfinish = 1000;
    Nh = 100;
    Tstart = 100;
    Tfinish = 1000;
    Th = 100;
        write(%io(2),'Linear');
        for N = Nstart:Nh:Nfinish
                for T = Tstart:Th:Tfinish
                        if T==N
                            mprintf('N = %d\n' , N);
                            mprintf('T = %d\n' , T);
                            tic();
                              MatrixVector(N, T);
                            mprintf('%.3f\n' , toc());
                                 
                            end
                 end
        end
        write(%io(2),'Parallel');
         for N = Nstart:Nh:Nfinish
                for T = Tstart:Th:Tfinish
                        if T==N
                            mprintf('N = %d\n' , N);
                            mprintf('T = %d\n' , T);
                            tic();  
                            parallel_run(N,T, "MatrixVector");  
                            mprintf('%.3f\n' , toc());
                        end
                 end
         end
endfunction

main();
function r=f(n)
   if(n == 0)
     r= 1;
   else
     r= n-m(f(n-1));
   end;
endfunction

function r=m(n)
   if(n == 0)
     r= 0;
   else
     r= n - f(m(n-1));
   end;
endfunction;

n_max=40;
tic();
t0=getdate();
for i=1:n_max
   r(i)= m(i);
end;
//etime(getdate(), t0); // output the wallclock time for the explicit loop computation
mprintf('%.3f\n' ,toc());
tic();
//t0=getdate();
r= parallel_run(1:n_max,"m");
// output the wallclock time for the parallel_run computation,
// it should be lower that the previous on multicore architectures.
//etime(getdate(), t0);
mprintf('%.3f\n' ,toc());
start = tic();
realtimeinit(1);
realtime(0);
realtime(10);
fin = toc();
mprintf('%.3f\n' , fin);
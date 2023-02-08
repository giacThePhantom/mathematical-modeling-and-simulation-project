[V1,M1,N1,t1,Mtot1,Ntot1] = ml_4_rtc();
[V2,M2,N2,t2,Mtot2,Ntot2] = ml_4_gill();

csvwrite("rtc.csv", cat(1, V1, M1, N1).')
csvwrite("gill.csv", cat(1, V2, M2, N2).')

figure
   plot(t1, V1, "g")
   hold on;
   plot(t2, V2, "b")


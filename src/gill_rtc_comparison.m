n_channels = ["params/1_1.csv", "params/10_10.csv", "params/50_50.csv", "params/100_100.csv"];

for i=1:length(n_channels)

   tau_mat = [ ; ; ; ]; 
   [V1,M1,N1,t1,tau_mat] = ml_4_rtc(n_channels(i), tau_mat, 0);

   tau_arr = [];
   eps_arr = [];
   [V2,M2,N2,t2,tau_arr, eps_arr] = ml_4_gill(n_channels(i), tau_arr, eps_arr);

   csvwrite(join([n_channels(i) "rtc"], "_"), cat(1, V1, M1, N1).')
   csvwrite(join([n_channels(i) "gill"], "_"), cat(1, V2, M2, N2).')

end


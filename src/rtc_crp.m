function statistics = compute_rtc_statistics(f1, f2)

    params1 = readtable(f1);
    params2 = readtable(f2);

    params1_arr = table2array(params1);
    params2_arr = table2array(params2);

    diff_params = params1_arr - params2_arr;
    norm = sqrt(sum(diff_params.^2));
    N = 5;
    spike_diff = 0;
    for i=1:N
        tau_mat = [ ; ; ; ];

        [V1,M1,N1,t1,tau_mat] = ml_4_rtc(f1, tau_mat, 0);
        [V2,M2,N2,t2,tau_mat] = ml_4_rtc(f2, tau_mat, 0);

        times1 = length(spike_detection(V1, 0));
        times2 = length(spike_detection(V2, 0));
        spike_diff = spike_diff + (times1 - times2);
    end

    result = 1/(norm*N)*spike_diff;
    fprintf("%f",result);
end

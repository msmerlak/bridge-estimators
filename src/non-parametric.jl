function Rt(X::RegularObservations; window_length = 5)
    [   mean(X.values[t - window_length::t])/mean(X.values[t - 1 - window_length:t-1])     
        for t in 1:X.n
    ]
end
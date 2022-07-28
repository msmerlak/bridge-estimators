import BridgeEstimators.RegularObservations

function Rt(X::RegularObservations; window_length = 5)
    [   mean(X.values[i - window_length:i, j])/mean(X.values[i - 1 - window_length:i-1, j])     
        for i in window_length + 2:size(X.values, 1), j in 1:size(X.values, 2)
    ]
end

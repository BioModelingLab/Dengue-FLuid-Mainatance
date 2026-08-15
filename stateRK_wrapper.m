function dx = stateRK_wrapper(t, x, params)
    dx = stateRK(x, params); % Call the original stateRK function
end

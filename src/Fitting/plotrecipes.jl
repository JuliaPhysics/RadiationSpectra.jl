@recipe function f(fit::FitFunction; npoints = 501, use_initial_parameters = false)
    x = collect(range(fit.fitrange[1], stop=fit.fitrange[2], length=npoints))
    par = use_initial_parameters ? fit.initial_parameters : fit.parameters
    y = fit.model(x, par)
    linecolor --> :red
    x,y
end


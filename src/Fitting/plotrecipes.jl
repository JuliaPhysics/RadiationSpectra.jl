@recipe function f(ff::FitFunction; npoints = 501, use_initial_parameters = false)
    x = collect(range(ff.fitrange[1], stop=ff.fitrange[2], length=npoints))
    par = use_initial_parameters ? ff.initial_parameters : ff.parameters
    y = ff.model(x, par)
    linecolor --> :red
    x,y
end


function lam = lambda(v,d,visc,delta)

    Re = abs(v)*d/visc;
    
    if Re < 2320
        lam = 64/Re;
    elseif Re < 10*d/delta
        lam = 0.3164 / Re^0.25;
    elseif Re < 500*d/delta
        lam = 0.11 * (68/Re + delta/d)^0.25;
    else
        lam = 0.11 * (delta/d)^0.25;
    end

end
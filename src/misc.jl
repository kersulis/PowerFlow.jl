import Base.LinAlg.scale!

"""
Scale all fields of `c` which have physical units by default by `s`. Convenient way to
switch between physical units and pu.
"""
function scale!(c::MatpowerCases.Case, s::Float64)
    fields = [:bus; :gen; :branch]
    for f in fields
        scale!(c.(f), s)
    end
end

function scale!(b::MatpowerCases.Bus, s::Float64)
    fields = [:PD; :QD]
    for f in fields
        scale!(b.(f), s)
    end
end

function scale!(g::MatpowerCases.Gen, s::Float64)
    fields = [:PG; :QG; :QMAX; :QMIN; :MBASE; :PMAX; :PMIN;
                :PC1; :PC2; :QC1MIN; :QC1MAX; :QC2MIN; :QC2MAX;
                :RAMP_AGC; :RAMP_10; :RAMP_30; :RAMP_Q]
    for f in fields
        scale!(g.(f), s)
    end
end

function scale!(b::MatpowerCases.Branch, s::Float64)
    fields = [:RATE_A; :RATE_B; :RATE_C; :PF; :QF; :PT; :QT]
    for f in fields
        scale!(b.(f), s)
    end
end

"""
Renumber elements of v so all instances of old[i] are replaced with new[i].
"""
function renumber!(v, old, new)
    for i in 1:length(v)
        v[i] = new[findfirst(old, v[i])]
    end
    nothing
end

function renumber!(c::MatpowerCases.Case, old, new)
    renumber!(c.bus.BUS_I, old, new)
    renumber!(c.gen.GEN_BUS, old, new)
    renumber!(c.branch.F_BUS, old, new)
    renumber!(c.branch.T_BUS, old, new)
    nothing
end

"""
Aggregate generation at each bus.
"""
function aggregate_gen(g::MatpowerCases.Gen, b::MatpowerCases.Bus)
    idx = b.BUS_I
    P = Vector{Float64}(length(idx))
    Q = Vector{Float64}(length(idx))
    for i in idx
        P[i] = sum(g.PG[findin(g.GEN_BUS, i)])
        Q[i] = sum(g.QG[findin(g.GEN_BUS, i)])
    end
    return P, Q
end

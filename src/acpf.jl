using MatpowerCases
using JLD

"""
An AC power flow routine based on Newton's method. Accommodates voltage and reactive power limits.

IN

* `V` [pu], Vector of voltage magnitudes
* `T` [rad], Vector of voltage angles
* `Ps` [pu], Vector of scheduled active power injections
* `Qs` [pu], Vector of scheduled reactive power injections
* `Pc` [pu], Vector of calculated active power injections (initialize to zeros)
* `Qc` [pu], Vector of calculated reactive power injections (initialize to zeros)
* `ty` [-], Vector of bus types (1: PQ, 2: PV, 3: slack)
* `Y` [pu], Admittance matrix
* `J`, Jacobian matrix (initialize to zeros)
* `Vmax` [pu], Vector of maximum voltage magnitudes for PQ buses
* `Qmin` [pu], Vector of minimum reactive power injections
* `Qmax` [pu], Vector of maximum reactive power injections
* `tol` [-], Newton iteration tolerance, 1e-5 by default
* `maxiter` [-], Maximum allowable iterations
* `silent` [-], Set to true to suppress messages

OUT: Modifies the following input arguments.

* `V` [pu], Vector of converged voltage magnitudes
* `T` [pu], Vector of converged voltage angles
* `Pc` [pu], Vector of calculated active power injections
* `Qc` [pu], Vector of calculated reactive power injections
* `Qs` [pu], May change if reactive power limit is encountered
* `ty` [-], If a PQ bus reaches its maximum voltage, the algorithm converts it
to a PV bus. If a PV bus encounters a reactive power limit, the algorithm converts it to a PQ bus.
"""
function acpf!{I<:Integer, F<:AbstractFloat}(
    V::Vector{Float64},
    T::Vector{Float64},
    Ps::Vector{Float64},
    Qs::Vector{Float64},
    Pc::Vector{Float64},
    Qc::Vector{Float64},
    ty::Vector{I},
    Y::Union{Array{Complex{Float64},2},SparseMatrixCSC{Complex{Float64},Int64}},
    J::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int64}};
    Vmax::Vector{F} = Float32[],
    Qmin::Vector{F} = Float32[],
    Qmax::Vector{F} = Float32[],
    tol::Float64 = 1e-5,
    maxiter::Int64 = 8,
    silent::Bool = false
    )

    nb = length(ty)

    # determine whether PQ nodes will switch to PV or vice versa
    Vlim = !isempty(Vmax)
    Qlim = !isempty(Qmin) || !isempty(Qmax)

    # initialize mismatch, update vectors
    mis = Vector{Float64}(2nb)
    update = Vector{Float64}(2nb)

    # Newton iteration
    converged = false
    iter = 0
    while iter < maxiter
        injection!(Pc, Qc, V, T, Y)
        mismatch!(mis, Ps, Pc, Qs, Qc, ty)

        # are we done?
        if maxabs(mis) < tol
            # see if any buses need to change type
            tycheck = deepcopy(ty)
            Vlim && pq2pv!(ty, V, Qc, Vmax; silent = silent)
            Qlim && pv2pq!(ty, V, Qs, Qc, Qmin, Qmax; silent = silent)

            # if no type changes, terminate
            if ty == tycheck
                converged = true
                !silent && info("Converged, $(iter) iterations")
                break
            else
                # recompute injection, mismatch
                continue
            end
        end

        # if not done, update and continue
        jacobian!(J, ty, V, T, Pc, Qc, Y)
        update[:] = lufact!(-J)\mis
        T[:] += update[1:nb]
        V[:] += update[nb+1:2nb]
        iter += 1
    end
    if !converged
        warn("Powerflow not converged")
        # save("dump.jld", "V", V, "T", T, "ty", ty, "Pc", Pc, "Qc", Qc, "mis", mis, "Ps", Ps, "Qs", Qs, "Qmin", Qmin, "Qmax", Qmax, "Y", Y, "J", J)
    end
    nothing
end

"""
Update real and reactive power injections.

IN:

* `Pc`, vector of active power injections [pu]
* `Qc`, vector of reactive power injections [pu]
* `V`, vector of bus voltages [pu]
* `T`, vector of bus angles [rad]
* `Y`, admittance matrix (sparse or dense)

OUT: modifies input arguments `Pc` and `Qc`.
"""
function injection!(Pc, Qc, V, T, Y)
    Pc[:] = 0.0
    Qc[:] = 0.0
    @inbounds for i in 1:length(V)
        for k in find(Y[i,:])
            Tik = T[i] - T[k]
            Gik, Bik = real(Y[i,k]), imag(Y[i,k])
            Pc[i] += V[i]*V[k]*(Gik*cos(Tik) + Bik*sin(Tik))
            Qc[i] += V[i]*V[k]*(Gik*sin(Tik) - Bik*cos(Tik))
        end
    end
    nothing
end

"""
Update power flow Jacobian.

IN:

* `J`, Jacobian (sparse or dense)
* `ty`, vector of bus types [1 for PQ, 2 for PV, 3 for reference]
* `V`, vector of bus voltages [pu]
* `T`, vector of bus angles [rad]
* `Pc`, vector of calculated active power injections [pu]
* `Qc`, vector of calculated reactive power injections [pu]
* `Y`, admittance matrix (sparse or dense)

OUT: modifies input argument `J`.
"""
function jacobian!(J, ty, V, T, Pc, Qc, Y)
    nb = length(V)
    dPda = sub(J, 1:nb, 1:nb)
    dPdV = sub(J, 1:nb, (nb + 1):(2nb))
    dQda = sub(J, (nb + 1):(2nb), 1:nb)
    dQdV = sub(J, (nb + 1):(2nb), (nb + 1):(2nb))

    @inbounds for i = 1:nb
        Vi = V[i]
        for k = find(Y[i,:])
            if i == k
                dPda[i,i] = -imag(Y[i,i])*Vi^2 - Qc[i]
                dPdV[i,i] =  real(Y[i,i])*Vi   + Pc[i]/Vi
                dQda[i,i] = -real(Y[i,i])*Vi^2 + Pc[i]
                dQdV[i,i] = -imag(Y[i,i])*Vi   + Qc[i]/Vi
            else
                Vk = V[k]
                Tik = T[i] - T[k]
                Gik, Bik = real(Y[i,k]), imag(Y[i,k])
                dPda[i,k] =  Vi*Vk*(Gik*sin(Tik) - Bik*cos(Tik))
                dPdV[i,k] =  Vi*(Gik*cos(Tik)    + Bik*sin(Tik))
                dQda[i,k] = -Vi*Vk*(Gik*cos(Tik) + Bik*sin(Tik))
                dQdV[i,k] =  Vi*(Gik*sin(Tik)    - Bik*cos(Tik))
            end
        end

        if ty[i] != 1 # only PQ buses have reactive mismatch
            dPdV[:,i] = 0.0
            dQda[i,:] = 0.0

            dQdV[i,:] = 0.0
            dQdV[:,i] = 0.0
            dQdV[i,i] = 1.0
        end
        if ty[i] == 3 # slack bus has no mismatch eqs
            dPda[i,:] = 0.0
            dPda[:,i] = 0.0
            dPda[i,i] = 1.0

            dPdV[i,:] = 0.0
            dQda[:,i] = 0.0
        end
    end
    nothing
end

function mismatch!(mis, Ps, Pc, Qs, Qc, ty)
    nb = length(Ps)
    Pm = sub(mis, 1:nb)
    Qm = sub(mis, (nb + 1):2nb)
    @inbounds for i = 1:nb
        Pm[i] = (ty[i] != 3 ? Pc[i] - Ps[i] : 0)
        Qm[i] = (ty[i] == 1 ? Qc[i] - Qs[i] : 0)
    end
    nothing
end

"""
Loop through PQ buses with specified voltage limits.
If a bus voltage exceeds its limit, switch it to a PV node.

IN:

* `ty`, vector of bus types [1 for PQ, 2 for PV, 3 for reference]
* `V`, vector of bus voltages [pu]
* `Qc`, vector of calculated reactive power injections [pu]
* `Vmax`, vector of maximum bus voltages [pu]

OUT: Modifies input arguments `ty` and `V`.
"""
function pq2pv!(ty, V, Qc, Vmax; silent = false)
    @inbounds for i in findin(ty, 1)
        if V[i] > Vmax[i]
            V[i] = Vmax[i]
            ty[i] = 2 # switch from PQ to PV bus
            !silent && info("Bus $i has hit voltage limit when Q = $(Qc[i])")
        end
    end
    nothing
end

"""
Loop through PV buses with specified reactive power limits.
If a bus reactive injection exceeds its limit, switch it to a PQ node.

Modiifies input arguments `ty` and `Qs`.

IN:

* `ty`, vector of bus types [1 for PQ, 2 for PV, 3 for reference]
* `V`, vector of bus voltages [pu]
* `Qs`, vector of specified reactive power injections [pu]
* `Qc`, vector of calculated reactive power injections [pu]
* `Qmin`, vector of minimum reactive power injections [pu]
* `Qmax`, vector of maximum reactive power injections [pu]

OUT: This function modifies input arguments `ty` and `Qs`.
"""
function pv2pq!(ty, V, Qs, Qc, Qmin, Qmax; silent = false)
    @inbounds for i in findin(ty, 2)
        if Qc[i] > Qmax[i]
            Qs[i] = Qmax[i]
            ty[i] = 1
            !silent && info("Bus $i has hit Qmax limit when V = $(V[i])")
        elseif Qc[i] < Qmin[i]
            Qs[i] = Qmin[i]
            ty[i] = 1
            !silent && info("Bus $i has hit Qmin limit when V = $(V[i])")
        end
    end
    nothing
end

"""
    Y = getY(f, t, r, x, b [, tap, ysh, lineOut])
This function builds admittance matrices.

INPUTS:

*  `f`: nline Vector of from buses (range 1 to nbus)
*  `t`: nline Vector of to buses (range 1 to nbus)
*  `r`: nline Vector of series resistances (pu)
*  `x`: nline Vector of series reactances (pu)
*  `b`: nline Vector of line charging shunt susceptances (pu)
*  `tap`: nline vector of complex tap ratios (pu)
*  `ysh`: nbus vector of complex bus shunt admittances (pu)
*  `lineOut`: [nline x 1, logical] which lines are out of service

for each branch, compute the elements of the branch admittance matrix where
```
  | If |   | Yff  Yft |   | Vf |
  |    | = |          | * |    |
  | It |   | Ytf  Ytt |   | Vt |
```
*Assumes the grid is connected (i.e. no island nodes).*
"""
function getY(
    from::Vector{Int64},
    to::Vector{Int64},
    r::Vector{Float64},
    x::Vector{Float64},
    b::Vector{Float64},
    tap::Vector{Complex{Float64}} = fill(Complex(1.), length(from)),
    ysh::Vector{Complex{Float64}} = fill(0.0im, length(unique([from;to])));
    lineOut::Vector{Bool} = fill(false, length(from)),
    reindex::Bool = true
    )

    if reindex
        old = sort(unique([from;to]))
        new = collect(1:length(old))
        f = copy(from)
        t = copy(to)
        renumber!(f, old, new)
        renumber!(t, old, new)
    else
        f = from
        t = to
    end

    nbus,nline = length(ysh),length(lineOut)
    inService = find(!lineOut)
    tap = tap[inService]
    fis = f[inService]
    tis = t[inService]

    Ys = (1./(r + x*im))[inService]
    Bc = b[inService]
    Ytt = Ys + Bc*im/2
    Yff = Ytt./(abs2(tap))
    Yft = -Ys./conj(tap)
    Ytf = -Ys./tap

    # connection matrix for line and from buses
    Cf = sparse(1:nline,f,ones(nline),nline,nbus)
    # connection matrix for line and to buses
    Ct = sparse(1:nline,t,ones(nline),nline,nbus)

    # build Yf and Yt such that Yf*V is the vector of complex branch currents injected.
    i = [inService;inService]
    Yf = sparse(i,[fis;tis],[Yff;Yft],nline,nbus)
    Yt = sparse(i,[fis;tis],[Ytf;Ytt],nline,nbus)

    # build Ybus
    Y = Cf'*Yf + Ct'*Yt + sparse(1:nbus,1:nbus,ysh,nbus,nbus)
    return full(Y) #, Yf, Yt
end

function acpf!(c::MatpowerCases.Case; silent = false)
    # convert to per unit
    scale!(c, 1/c.baseMVA)
    b, g = c.bus, c.gen
    # convert to rad
    b.VA[:] = deg2rad(b.VA)

    # convert to internal numbering
    old = copy(b.BUS_I)
    new = collect(1:length(old))
    old != new && renumber!(c, old, new)

    nb = length(b.PD)
    PG, QG = aggregate_gen(g, b)

    Ps = PG - b.PD
    Qs = QG - b.QD

    Pc = zeros(nb)
    Qc = zeros(nb)

    Y = getY(c; reindex = false)
    J = zeros(2nb, 2nb)

    Qmin = zeros(nb)
    Qmin[g.GEN_BUS] = g.QMIN
    Qmax = zeros(nb)
    Qmax[g.GEN_BUS] = g.QMAX

    acpf!(b.VM, b.VA, Ps, Qs, Pc, Qc,
        b.BUS_TYPE, Y, J,
        Vmax = b.VMAX, Qmin = Qmin, Qmax = Qmax, silent = silent)

    # update generation at slack bus
    if length(findin(b.BUS_TYPE, 3)) > 1
        warn("More than one slack bus!")
    else
        bidx = findfirst(b.BUS_TYPE, 3)
        gidx = findfirst(g.GEN_BUS, bidx)
        g.PG[gidx] = Pc[bidx]
        g.QG[gidx] = Qc[bidx]
    end

    # update reactive generation at PV buses
    bidx = findin(b.BUS_TYPE, 2)
    for i in bidx
        gidx = findfirst(g.GEN_BUS, i)
        g.QG[gidx] = Qc[i]
    end

    # convert back to physical units
    scale!(c, c.baseMVA)
    # convert back to degrees
    b.VA[:] = rad2deg(b.VA)

    # convert back to old bus labels
    old != new && renumber!(c, new, old)
    nothing
end

function getY(c::MatpowerCases.Case; reindex = true)
    b = c.branch
    s = c.bus
    ysh = (s.GS + s.BS*1im )/c.baseMVA
    tap = Vector{Complex{Float64}}(b.TAP + 1.0)
    getY(b.F_BUS, b.T_BUS, b.BR_R, b.BR_X, b.BR_B, tap, ysh; lineOut = !b.BR_STATUS, reindex = reindex)
end

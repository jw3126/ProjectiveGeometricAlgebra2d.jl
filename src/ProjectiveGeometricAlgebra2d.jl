module ProjectiveGeometricAlgebra2d
#https://bivector.net/2DPGA.pdf
using StaticArrays
using ArgCheck
using Random

struct PGA2d{T}
    e::T
    e0::T
    e1::T
    e2::T
    e12::T
    e02::T
    e01::T
    e012::T
end

function Base.show(io::IO, o::PGA2d)
    print(io, "PGA2d")
    print(io, "(")
    pieces = String[]
    for pname in propertynames(o)
        val = getproperty(o, pname)
        if val != 0
            push!(pieces, "$pname=$val")
        end
    end
    print(io, join(pieces, ", "))
    print(io, ")")
end

function _serialize(x::PGA2d)
    @SVector[
        x.e,x.e0,x.e1,x.e2,x.e12,x.e02,x.e01,x.e012,
    ]
end

function Base.isapprox(x::PGA2d, y::PGA2d; kw...)
    isapprox(_serialize(x), _serialize(y); kw...)
end

@inline function PGA2d(;e=0, e0=0, e1=0, e2=0, e12=0, e02=0, e01=0, e012=0)
    e, e0, e1, e2, e12, e02, e01, e012 = promote(e, e0, e1, e2, e12, e02, e01, e012)
    PGA2d(
        e   ,
        e0  ,
        e1  ,
        e2  ,
        e12 ,
        e02 ,
        e01 ,
        e012,
    )
end

function wedge(x,y)
    PGA2d(
        e = x.e*y.e,

        e0 = x.e*y.e0 + x.e0*y.e,
        e1 = x.e*y.e1 + x.e1*y.e,
        e2 = x.e*y.e2 + x.e2*y.e,

        e12 = x.e*y.e12 + x.e1*y.e2 - x.e2*y.e1 + x.e12*y.e,
        e02 = x.e*y.e02 + x.e0*y.e2 - x.e2*y.e0 + x.e02*y.e,
        e01 = x.e*y.e01 + x.e0*y.e1 - x.e1*y.e0 + x.e01*y.e,

        e012 = x.e*y.e012 +
        x.e0*y.e12 - x.e1*y.e02 + x.e2*y.e01 +
        x.e12*y.e0 - x.e02*y.e1 + x.e01*y.e2 +
        x.e012 * y.e,
    )
end

function geomul(x,y)
    PGA2d(
        e = x.e*y.e + x.e1*y.e1 + x.e2*y.e2 - x.e12*y.e12,

        e0 = x.e*y.e0 + x.e0*y.e - x.e1*y.e01 - x.e2*y.e02 - x.e12*y.e012 + x.e02*y.e2 + x.e01*y.e1 - x.e012*y.e12,
        e1 = x.e*y.e1 + x.e1*y.e - x.e2*y.e12 + x.e12*y.e2,
        e2 = x.e*y.e2 + x.e2*y.e + x.e1*y.e12 - x.e12*y.e1,

        e12 = x.e*y.e12 + x.e1*y.e2 - x.e2*y.e1 + x.e12*y.e,
        e02 = x.e*y.e02 + x.e0*y.e2 - x.e1*y.e012 - x.e2*y.e0 - x.e12*y.e01 + x.e02*y.e + x.e01*y.e12 - x.e012*y.e1,
        e01 = x.e*y.e01 + x.e0*y.e1 - x.e1*y.e0 + x.e2*y.e012 + x.e12*y.e02 - x.e02*y.e12 + x.e01*y.e + x.e012*y.e2,
        e012 = x.e*y.e012 + x.e0*y.e12 - x.e1*y.e02 + x.e2*y.e01 + x.e12*y.e0 - x.e02*y.e1 + x.e01*y.e2 + x.e012*y.e,
    )
end

@inline function dual(x::PGA2d)
    PGA2d(
        e   =  x.e012,
        e0  =  x.e12,
        e1  = -x.e02,
        e2  =  x.e01,
        e12 =  x.e0,
        e02 = -x.e1,
        e01 =  x.e2,
        e012= x.e,
    )
end

function reverse(x::PGA2d)
    PGA2d(
        e    = e    ,

        e0   = e0   ,
        e1   = e1   ,
        e2   = e2   ,

        e12  = -e12 ,
        e02  = -e02 ,
        e01  = -e01 ,

        e012 = -e012,
    )

end

function pga_point(v)
    @argcheck length(v) == 2
    x,y = v
    PGA2d(e02=-x, e01=y, e12=1)
end

function pga_line(line)
    # ax + by + c = 0
    PGA2d(e1=line.a, e2=line.b, e0=line.c)
end

function euclidean_line_abc(pga)
    (a=pga.e1, b=pga.e2, c=pga.e0)
end

function euclidean_point(pga)
    w = pga.e12
    x = - pga.e02 / w
    y = pga.e01 / w
    @SVector[x,y]
end

function join_points(pt1, pt2)
    l_pga = anti_wedge(pga_point(pt1), pga_point(pt2))
    euclidean_line_abc(l_pga)
end

function meet_lines(l1, l2)
    pt_pga = wedge(pga_line(l1), pga_line(l2))
    return euclidean_point(pt_pga)
end

function anti_wedge_using_dual(x1,x2)
    dual(wedge(dual(x1), dual(x2)))
end

function randn_pga(rng=Random.GLOBAL_RNG)
    blades = ntuple(_->randn(rng), Val(8))
    PGA2d(blades...)
end

function rand_pga(args...)
    blades = ntuple(_->rand(args...), Val(8))
    PGA2d(blades...)
end

function anti_wedge(x1,x2)
    # TODO more efficient implementation
    anti_wedge_using_dual(x1,x2)
end

function pseudo_one(T=Int)
    PGA2d(e012=one(T))
end

end#module

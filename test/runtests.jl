using ProjectiveGeometricAlgebra2d
using Test
using ProjectiveGeometricAlgebra2d: PGA2d
using ProjectiveGeometricAlgebra2d: wedge, anti_wedge, geomul
using ProjectiveGeometricAlgebra2d: dual, reverse
using ProjectiveGeometricAlgebra2d: euclidean_point, euclidean_line_abc
using ProjectiveGeometricAlgebra2d: pga_point, pga_line
using ProjectiveGeometricAlgebra2d: meet_lines, join_points
using ProjectiveGeometricAlgebra2d: randn_pga

function test_associative(op, x, y, z)
    lhs = op(op(x,y), z)
    rhs = op(x, op(y, z))
    @test typeof(lhs)  == typeof(rhs)
    if !(lhs ≈ rhs)
        @show op
        for pname in propertynames(lhs)
            val_lhs = getproperty(lhs, pname)
            val_rhs = getproperty(rhs, pname)
            if !(val_lhs ≈ val_rhs)
                msg = """
                lhs.$(pname) = $(val_lhs)
                rhs.$(pname) = $(val_rhs)
                """
                println(msg)
            end
        end
    end
    @test lhs ≈ rhs
end

function test_associative(op)

    x = randn_pga()
    y = randn_pga()
    z = randn_pga()
    test_associative(op,x,y,z)
end

@testset "conversion" begin
    @testset "point" begin
        pt_euclid = randn(2)
        pt_pga = pga_point(pt_euclid)
        @test euclidean_point(pt_pga) == pt_euclid
    end
    @testset "line" begin
        l_euclid = NamedTuple{(:a,:b,:c)}(randn(3))
        l_pga = pga_line(l_euclid)
        @test euclidean_line_abc(l_pga) == l_euclid
    end
end

@testset "join and meet" begin

    @test join_points([0,1], [1,0]) === (a=1,b=1,c=-1)
    @test join_points([0,0], [1,0]) === (a=0,b=1,c=0) # TODO correct sign?
    @test join_points([0,1], [0,0]) === (a=1,b=0,c=0) # TODO correct sign?

    l1 = (a=1, b=0, c=-2)
    l2 = (a=0, b=1, c=-1)
    @test meet_lines(l1, l2) == [2,1]
end

@testset "mwe geomul" begin
    @test geomul(PGA2d(e1=1), PGA2d(e1=1)) == PGA2d(e=1)

end

@testset "mwe associativity geomul" begin
    x = PGA2d(e2=1)
    y = PGA2d(e12=1, e012=1)
    z = PGA2d( e02=1)
    xy = geomul(x,y)
    @test xy == PGA2d(e1=-1, e01=1)
    xy_z = geomul(xy, z)
    @test xy_z == PGA2d(e012=1)

    yz = geomul(y,z)
    @test yz == PGA2d(e01=1)
    x_yz = geomul(x, yz)
    @test x_yz == PGA2d(e012=1)
    @test xy_z == x_yz
end

@testset "mwe associtivity geomul" begin
    x = PGA2d(e01=1)
    y = PGA2d(e1=1)
    z = PGA2d(e1=1)

    xy = geomul(x,y)
    @test xy == PGA2d(e0=1)
    xy_z = geomul(xy, z)
    @test xy_z == PGA2d(e01=1)
    yz = geomul(y,z)
    @test yz == PGA2d(e=1)
    x_yz = geomul(x, yz)
    @test x_yz == PGA2d(e01=1)
    @test x_yz == xy_z
end

@testset "mwe associativity wedge" begin
    x = PGA2d(e=1, e1=1, e12=1)
    y = PGA2d(e=1, e0=1, e1=1, e2=1, e01=1)
    z = PGA2d(e0=1, e1=1, e02=1, e01=1)
    @test wedge(wedge(x,y), z) == wedge(x, wedge(y,z))
end

@testset "associativity" begin

    for _ in 1:100
        test_associative(wedge)
        test_associative(anti_wedge)
        test_associative(geomul)
    end
end

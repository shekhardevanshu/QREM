using Printf
using DelimitedFiles, Statistics
using FLoops
using LinearAlgebra

function binaryRead(fname)
    data = Array{Float64}(undef, filesize(fname)÷8)
    read!(fname, data)
    return data
end

function main(args)
    
    # barr = sort(vcat(collect(0.0:1:5), collect(10:5:20), collect(30:20:100),  collect(150:50:256), [120, 180, 500], collect(7:2:14), collect(0.1:0.2:1))) # 10, 12
    # barr = sort(vcat([0.1], collect(1:1:5), collect(10:5:20), collect(30:20:100), collect(35:20:100), collect(150:50:256), [120, 180], collect(500:200:1500), collect(7:2:14))) # 14
    barr = sort(vcat([0.1], collect(1:1:5), collect(10:5:20), collect(30:20:100), collect(25:10:50), collect(150:50:256), [120, 180, 350, 450], collect(500:200:1500), collect(7:2:14))) # 16
    # barr = sort(vcat(collect(0.1:0.2:10), collect(10:1:30), collect(35:5:100))) # 12

    bname = map(x -> @sprintf("%0.1f", x), barr)
    
    l = parse(Int, args[1])
    N = 2^l
    nconv = parse(Int, args[2])
    s = parse(Int, args[3])
    # a = parse(Float64, args[4])
    # c = parse(Float64, args[5])
    # a = @sprintf("%0.2f", a)
    # c = @sprintf("%0.1f", c)
    e = "0.0"

    ipr_s = Array{Float64}(undef, s*nconv)
    ipr_g = Matrix{Float64}(undef, length(barr), length(ipr_s))
    ipr_pgi = Vector{Float64}(undef, nconv)

    for (j, g) in enumerate(bname)

        @floop for m in 1:s
            
            eVec = binaryRead("./dataBinQREM_L=$(l)/eigvec_L=$(l)_b=$(g)_itr_$(m)_e=$e.dat")
	    # eVec = binaryRead("./dataBinQREM_L=$l/eigvec_L=$(l)_b=$(g)_a=$(a)_g=$(c)_itr_$(m)_e=$e.dat")
            
            for n in 1:nconv
                ipr_pgi[n] = sum((eVec[(n-1)*N+1:n*N]).^4)
            end

            ipr_s[(m-1)*nconv+1:nconv*m] = ipr_pgi

        end

	print("completed g = $g ipr = $(mean(ipr_s))\n")
        flush(stdout)

        ipr_g[j,:] = ipr_s
        
    end

    # writedlm("ipr_qrem_pe,L=$l,nev=$nconv,a=$a,g=$c,e=$e.txt", ipr_g)
    writedlm("ipr_qrem_pe,L=$l,nev=$nconv,e=$e.txt", [ipr_g barr])

end

main(ARGS)

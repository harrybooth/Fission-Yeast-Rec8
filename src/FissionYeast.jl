mutable struct MitoticCell # Haploid
    rm :: Float64
    rp :: Float64
    gamete_fate :: Any
end

mutable struct MGam # Haploid
    rm :: Float64
    rp :: Float64
end

mutable struct PGam # Haploid
    rm :: Float64
    rp :: Float64
end

mutable struct Zygote # Diploid
    rm :: Float64
    rp :: Float64
end

mutable struct Spore # Haploid
    rm :: Float64
    rp :: Float64
end

mutable struct PopulationParameters
    max_population_size :: Int
    max_population_pairs :: Int
    n_mitosis_rounds :: Int
    p_sister_mating :: Float64
    mutation_distribution :: Any
    mutation_cap :: Any
end

######################

function mitosis(FY::Union{MitoticCell,MGam,PGam}, noise_distribution :: Distribution, cap::Float64)
    
    rm = max(0.,min(FY.rm + rand(noise_distribution),cap))
    rp = max(0.,min(FY.rp + rand(noise_distribution),cap)) 

    return [MitoticCell(rm,rp,:P),MitoticCell(rm,rp,:M)] # should sisters be genetically identical. with average inheritance and sister mating effectivly asexual
end

function gametogenesis(FY::MitoticCell)

    if FY.gamete_fate == :P
        return PGam(FY.rm,FY.rp)
    else
        return MGam(FY.rm,FY.rp) 
    end
end

function mating_success(G1::Union{MGam,PGam},G2::Union{MGam,PGam},Pm,p)

    if typeof(G1) == MGam
        prob = Pm(G1.rm,G2.rp,p)
    else
        prob = Pm(G2.rm,G1.rp,p)
    end

    prob
end

function survive_mitosis(G1::Union{MGam,PGam},q,p)

    if typeof(G1) == MGam
        prob = q(G1.rm,p)
    else
        prob = q(G1.rp,p)
    end

    prob
end

function fusion(G1::Union{MGam,PGam},G2::Union{MGam,PGam})
    Zygote((G1.rm + G2.rm)/2,(G1.rp + G2.rp)/2)
end

function sporulate(Z::Zygote)
    return [Spore(Z.rm,Z.rp) for _ in 1:4] # 4 spores
end

function germinate(S::Spore)
    return MitoticCell(S.rm,S.rp,nothing)
end

#######################

function create_gamete_mate_pairs(gamete_population::Vector{Vector{Any}},PP)

    N_gamete_pairs = length(gamete_population)

    sister_mating_id = rand(Uniform(0,1),N_gamete_pairs) .< PP.p_sister_mating

    with_sisters = gamete_population[sister_mating_id]
    with_others = gamete_population[.! sister_mating_id]

    if length(with_others) > 1
        with_others_paired = vcat([[with_others[end][1],with_others[1][2]]], [[with_others[i][1],with_others[i+1][2]] for i in 1:length(with_others)-1])

        mating_pairs = vcat(with_sisters,with_others_paired)
    else
        mating_pairs = with_sisters
    end

    return mating_pairs

end

#######################

function fission_yeast_lifecycle(population,Pm,q,p,PP::PopulationParameters)

    for i in 1:PP.n_mitosis_rounds
        population_v = reduce(vcat,population)
        population_new = [mitosis(i,PP.mutation_distribution,PP.mutation_cap) for i in population_v]
        population = sample(population_new,PP.max_population_pairs)
    end

    gamete_population = [gametogenesis.(i) for i in population]

    if PP.p_sister_mating != 0
        mating_pairs = create_gamete_mate_pairs(gamete_population,PP)
    else
        mating_pairs = gamete_population
    end

    returning_daughters = []
    returning_spores = []

    for gametes in mating_pairs

        if rand() < p[:f]

            zygote = fusion(gametes[1],gametes[2])
            spores = sporulate(zygote)

            for spore in spores
                if rand() < mating_success(gametes[1],gametes[2],Pm,p)
                    push!(returning_spores,germinate(spore))
                end
            end

        else
            for g in gametes
                for daughter in mitosis(g,PP.mutation_distribution,PP.mutation_cap) 
                    if rand() < survive_mitosis(g,q,p)
                        push!(returning_daughters,daughter)
                    end
                end
            end
        end

    end

    return sample(vcat(returning_daughters,returning_spores),PP.max_population_size)
end

########################

function evolution(initial_spore,N_gen,Pm,q,p,PP,record_population_state = false)

    population = [germinate(initial_spore) for i in 1:PP.max_population_size]

    all_rp = []
    all_rm = []

    rp = []
    rm = []

    gen = 0

    while (gen < N_gen)

        population = fission_yeast_lifecycle(population,Pm,q,p,PP)

        push!(all_rp,mean(map(cell->cell.rp,population)))
        push!(all_rm,mean(map(cell->cell.rm,population)))

        if record_population_state
            push!(rp,map(cell->cell.rp,population))
            push!(rm,map(cell->cell.rm,population))
        end
        
        gen +=1
    end

    return all_rp,all_rm,rp,rm
end

########################

Pm(rm,rp,p) = p[:κm]*((rm + rp)^p[:hm] / (p[:am]^p[:hm] + (rm + rp)^p[:hm])) + (1-p[:κm])*p[:Δm]
Pm(r,p) = p[:κm]*(r^p[:hm] / (p[:am]^p[:hm] + r^p[:hm])) + (1-p[:κm])*p[:Δm]

q(ri,p) = p[:κq]*(1 - (ri^p[:hq] / (p[:aq]^p[:hq] + ri^p[:hq]))) + (1-p[:κq])*p[:Δq]

#######################

J(p) = -((1-p.f)*p.κq*p.hq*p.aq^p.hq)/(4*p.f*p.κm*p.hm*p.am^p.hm)

deriv_fitness_p(rp,rm,p) = J(p)*rp^(p.hq-1)*(p.am^p.hm + (rp+rm)^p.hm)^2 + (p.aq^p.hq + rp^p.hq)^2
deriv_fitness_m(rp,rm,p) = J(p)*rm^(p.hq-1)*(p.am^p.hm + (rp+rm)^p.hm)^2 + (p.aq^p.hq + rm^p.hq)^2

dvect(x,p) = Point2f([deriv_fitness_m(x[1],x[2],p),deriv_fitness_p(x[1],x[2],p)])

dvect_grad(G,u,p) = -1 .* [deriv_fitness_m(u[1],u[2],p),deriv_fitness_p(u[1],u[2],p)]
# -------------------------------------------------- VARIABLES --------------------------------------------------

# -------------------------------------------------- CONSTANTS --------------------------------------------------
const BASE_COMPLEMENT = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')


# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# returns the reversed complement codon set
function create_complement_reversed_codons(data::CodonGraphData)
    temp_codons = reverse_codon_set(complement_codon_set(data.codon_set))
    @debug "Original codon set: $(data.codon_set)
    -> reversed complement codon set: $temp_codons"
    return temp_codons
end

# returns the complement codon set
function complement_codon_set(codons::Vector{String})
    complement_codons = Vector{String}()
    for codon in codons
        # add the reversed complement codon to the reversed_codons set
        push!(complement_codons, complement_codon(codon))
    end
    @debug "Original codon set: $codons -> complement codon set: $complement_codons"
    return complement_codons
end


# returns the reversed codon set
function reverse_codon_set(codons::Vector{String})
    reversed_codons = Vector{String}()
    for codon in codons
        # add the reversed complement codon to the reversed_codons set
        push!(reversed_codons, reverse_codon(codon))
    end
    @debug "Original codon set: $codons -> reversed codon set: $reversed_codons"
    return reversed_codons
end


# returns the complement base
function complement_base(base::Char)
    @assert base in keys(BASE_COMPLEMENT)
    "Base is invalid. Only A, C, G, T are allowed."
    @debug "Original base: $base, -> complement base: $(BASE_COMPLEMENT[base])"
    return BASE_COMPLEMENT[base]
end


# returns the complement codon
function complement_codon(codon::AbstractString)
    @assert all(c in keys(BASE_COMPLEMENT) for c in codon)
    "Codon contains invalid characters. Only A, C, G, T are allowed."

    if length(codon) == 3
        @debug "Original codon: $codon, -> complement codon: $(String([complement_base(c) for c in codon]))"
        return String([complement_base(c) for c in codon])
    end
end


# returns the reversed codon
function reverse_codon(codon::AbstractString)
    @assert all(c in keys(BASE_COMPLEMENT) for c in codon)
    "Codon contains invalid characters. Only A, C, G, T are allowed."

    if length(codon) == 3
        reversed = reverse(codon)
        @debug "Original codon: $codon, -> reversed codon: $reversed"
        return reversed
    end
end
test_codon1 = "ACG"
test_codon2 = "CGT"

# returns the anticodon of a given Codon
function complement(seq::AbstractString)
    if length(seq) == 3
        comp = map(c -> Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')[c], seq)
        return reverse(String(comp))
    end
end

# @debug "Anticodon für ", test_codon1, " ist ", complement(test_codon1))
# @debug "Anticodon für ", test_codon2, " ist ", complement(test_codon2))

code_set_true = ["ACG", "CGT",]
code_set_false = ["ACC", "CGT",]

function is_self_complementary(set::AbstractArray)
    for codon in set
        if complement(codon) in set
            @debug "Das Codon ", codon, " hat ein komplementäres Gegenstück in der Menge: ", complement(codon)
        else
            @debug "Das Codon ", codon, " hat kein komplementäres Gegenstück in der Menge. Es fehlt: ",
            complement(codon)
            return false
        end
    end
    return true
end

is_self_complementary(code_set_false)
is_self_complementary(code_set_true)
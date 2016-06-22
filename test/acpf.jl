using PowerFlow, MAT, MatpowerCases
m = extract_case(matread("../data/acpf_case9.mat")["results"])
j = extract_case("case9")
acpf!(j)

for f in fieldnames(Bus)
    if !all(isnan(m.bus.(f)))
        println(f, "\t\t", maxabs(j.bus.(f) - m.bus.(f)))
    end
end

for f in fieldnames(Branch)
    if !all(isnan(m.branch.(f)))
        println(f, "\t\t", maxabs(j.branch.(f) - m.branch.(f)))
    end
end

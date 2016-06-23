using PowerFlow, MAT, MatpowerCases
m = extract_case(matread("../data/acpf_case9.mat")["results"])
j = extract_case("case9")
acpf!(j)

for f in fieldnames(Bus)
    if !all(isnan(j.bus.(f)))
        println(f, "\t\t", maxabs(j.bus.(f) - m.bus.(f)))
    end
end

for f in fieldnames(Branch)
    if !all(isnan(j.branch.(f)))
        println(f, "\t\t", maxabs(j.branch.(f) - m.branch.(f)))
    end
end

for f in fieldnames(Gen)
    if !all(isnan(j.gen.(f)))
        println(f, "\t\t", maxabs(j.gen.(f) - m.gen.(f)))
    end
end

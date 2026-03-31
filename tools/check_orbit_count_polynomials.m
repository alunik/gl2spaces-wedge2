load "src/paper_data.m";

R<q> := PolynomialRing(Integers());
K := FieldOfFractions(R);

forward KernelOrderNoRadicalsPolynomial;
forward KernelOrderPolynomial;
forward StabilizerOrderPolynomial;
forward OrbitSizePolynomial;

GLOrderPolynomial := function(n)
    out := R!1;
    for i in [0..n - 1] do
        out *:= q^n - q^i;
    end for;
    return out;
end function;

GaussianTwoSubspacesPolynomial := function(m)
    num := (q^m - 1) * (q^m - q);
    den := (q^2 - 1) * (q^2 - q);
    quo, rem := Quotrem(num, den);
    assert rem eq 0;
    return quo;
end function;

SupportKeyPolynomial := function(block)
    if block`kind eq "J" then
        return "R:" cat block`label;
    elif block`kind eq "P" then
        return Sprintf("P:%o", block`degree);
    end if;
    error "SupportKeyPolynomial only applies to regular blocks";
end function;

RegularClustersPolynomial := function(blocks)
    regular := [block : block in blocks | block`kind ne "S"];
    if #regular eq 0 then
        return [];
    end if;

    clusters := [];
    current := [regular[1]];
    current_key := SupportKeyPolynomial(regular[1]);
    for i in [2..#regular] do
        key := SupportKeyPolynomial(regular[i]);
        if key eq current_key then
            Append(~current, regular[i]);
        else
            Append(~clusters, current);
            current := [regular[i]];
            current_key := key;
        end if;
    end for;
    Append(~clusters, current);
    return clusters;
end function;

QuotientOrderPolynomial := function(blocks)
    qtype := FiniteQuotientStringFromBlocks(blocks);
    if qtype eq "GL_2(q)" then
        return GLOrderPolynomial(2);
    elif qtype eq "B(q)" then
        return q * (q - 1)^2;
    elif qtype eq "T_s(q)" then
        return (q - 1)^2;
    elif qtype eq "N_s(q)" then
        return 2 * (q - 1)^2;
    elif qtype eq "N_ns(q)" then
        return 2 * (q^2 - 1);
    elif qtype eq "~C_2(q)" then
        return 2 * (q - 1);
    elif qtype eq "~C_3(q)" then
        return 3 * (q - 1);
    elif qtype eq "~S_3(q)" then
        return 6 * (q - 1);
    end if;
    error "Unsupported quotient type";
end function;

Sp4OrderPolynomial := function()
    return q^4 * (q^2 - 1) * (q^4 - 1);
end function;

KernelOrderNoRadicalsPolynomial := function(blocks)
    if &and[block`kind eq "S" : block in blocks] then
        if #blocks eq 1 then
            d := blocks[1]`d;
            return q^(2 * d - 2) * (q - 1);
        elif #blocks eq 2 and blocks[1]`d eq 2 and blocks[2]`d eq 2 then
            return q^6 * GLOrderPolynomial(2);
        end if;
        error "Unsupported pure singular core";
    end if;

    nontrivial_s := [block : block in blocks | block`kind eq "S" and block`d gt 1];
    regular_core := [block : block in blocks | block`kind ne "S"];

    if #nontrivial_s eq 0 then
        order := R!1;
        for cluster in RegularClustersPolynomial(blocks) do
            if cluster[1]`kind eq "P" then
                m := cluster[1]`degree;
                order *:= q^m * (q^(2 * m) - 1);
            else
                ds := [block`d : block in cluster];
                if #ds eq 1 then
                    d := ds[1];
                    order *:= q^(3 * d - 2) * (q^2 - 1);
                elif ds eq [1, 1] then
                    order *:= Sp4OrderPolynomial();
                elif ds eq [2, 1] then
                    order *:= q^9 * (q^2 - 1)^2;
                else
                    error "Unsupported regular cluster";
                end if;
            end if;
        end for;
        return order;
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 2 and #regular_core gt 0 then
        return q^(TotalDimension(regular_core) + 2) * (q - 1) *
            KernelOrderNoRadicalsPolynomial(regular_core);
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 3 and #regular_core eq 1 and
        regular_core[1]`kind eq "J" and regular_core[1]`d eq 1 then
        return q^7 * (q - 1) * (q^2 - 1);
    end if;

    error "Unsupported non-radical family";
end function;

KernelOrderPolynomial := function(blocks)
    r, core := SplitOffS1Blocks(blocks);
    if r gt 0 then
        return q^(r * TotalDimension(core)) * GLOrderPolynomial(r) *
            KernelOrderNoRadicalsPolynomial(core);
    end if;
    return KernelOrderNoRadicalsPolynomial(blocks);
end function;

StabilizerOrderPolynomial := function(blocks)
    return KernelOrderPolynomial(blocks) * QuotientOrderPolynomial(blocks);
end function;

OrbitSizePolynomial := function(n, blocks)
    return (K!GLOrderPolynomial(n)) / (K!StabilizerOrderPolynomial(blocks));
end function;

for n in [4..7] do
    total := K!0;
    for T in FiniteOrbitTypesNLe7(n) do
        total +:= OrbitSizePolynomial(n, T`blocks);
    end for;

    m := n * (n - 1) div 2;
    target := GaussianTwoSubspacesPolynomial(m);
    assert total eq K!target;

    printf "[ok] n=%o: orbit-stabilizer sum matches the Grassmann polynomial for 2-spaces in F_q^%o\n", n, m;
    printf "     target polynomial = %o\n", Factorization(target);
end for;

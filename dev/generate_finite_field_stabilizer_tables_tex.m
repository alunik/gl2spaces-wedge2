SetColumns(0);
load "src/paper_data.m";

ReplaceAllStrings := function(s, old, new)
    if old eq "" then
        return s;
    end if;

    out := "";
    rest := s;
    pos := Position(rest, old);
    while pos ne 0 do
        if pos gt 1 then
            out cat:= rest[1..pos - 1];
        end if;
        out cat:= new;
        start := pos + #old;
        if start le #rest then
            rest := rest[start..#rest];
            pos := Position(rest, old);
        else
            rest := "";
            pos := 0;
        end if;
    end while;

    return out cat rest;
end function;

LatexTypeString := function(s)
    out := s;
    out := ReplaceAllStrings(out, "lambda", "\\lambda");
    out := ReplaceAllStrings(out, "mu", "\\mu");
    out := ReplaceAllStrings(out, "nu", "\\nu");
    out := ReplaceAllStrings(out, "!=", "\\neq");
    return "$" cat out cat "$";
end function;

SeqEqual := function(a, b)
    return #a eq #b and &and[a[i] eq b[i] : i in [1..#a]];
end function;

MatchesTwoSupports := function(supports, s1, s2)
    return (#supports eq 2) and
        ((SeqEqual(supports[1], s1) and SeqEqual(supports[2], s2)) or
         (SeqEqual(supports[1], s2) and SeqEqual(supports[2], s1)));
end function;

BlockDimension := function(block)
    if block`kind eq "S" then
        return 2 * block`d - 1;
    elif block`kind eq "J" then
        return 2 * block`d;
    elif block`kind eq "P" then
        return 2 * block`degree;
    end if;

    error "Unknown block kind";
end function;

TotalDimension := function(blocks)
    return &+[BlockDimension(block) : block in blocks];
end function;

SplitOffS1Blocks := function(blocks)
    r := 0;
    core := [];
    for block in blocks do
        if block`kind eq "S" and block`d eq 1 then
            r +:= 1;
        else
            Append(~core, block);
        end if;
    end for;
    return r, core;
end function;

RegularSupportData := function(blocks)
    a_support := [];
    b_support := [];
    c_support := [];
    poly_degrees := [];

    for block in blocks do
        if block`kind eq "J" then
            if block`label eq "a" then
                Append(~a_support, block`d);
            elif block`label eq "b" then
                Append(~b_support, block`d);
            elif block`label eq "c" then
                Append(~c_support, block`d);
            end if;
        elif block`kind eq "P" then
            Append(~poly_degrees, block`degree);
        end if;
    end for;

    Sort(~a_support);
    Reverse(~a_support);
    Sort(~b_support);
    Reverse(~b_support);
    Sort(~c_support);
    Reverse(~c_support);

    supports := [seq : seq in [a_support, b_support, c_support] | #seq gt 0];
    return supports, poly_degrees;
end function;

QuotientStringFromBlocks := function(blocks)
    supports, poly_degrees := RegularSupportData([block : block in blocks | block`kind ne "S"]);

    if #supports eq 0 and #poly_degrees eq 0 then
        return "$\\GL_2(q)$";
    end if;

    if #poly_degrees eq 0 then
        if #supports eq 1 then
            return "$B(q)$";
        elif #supports eq 2 then
            if MatchesTwoSupports(supports, [1], [1]) then
                return "$N_s(q)$";
            end if;
            return "$T_s(q)$";
        elif #supports eq 3 then
            return "$\\widetilde{S}_3(q)$";
        end if;
    end if;

    if #poly_degrees eq 1 and poly_degrees[1] eq 2 then
        if #supports eq 0 then
            return "$N_{ns}(q)$";
        elif #supports eq 1 then
            return "$\\widetilde{C}_2(q)$";
        end if;
    end if;

    if #poly_degrees eq 1 and poly_degrees[1] eq 3 and #supports eq 0 then
        return "$\\widetilde{C}_3(q)$";
    end if;

    error "Unsupported quotient pattern";
end function;

FqPowerString := function(m)
    return Sprintf("\\Fq^{\\,%o}", m);
end function;

RegularKernelString := function(blocks)
    supports, poly_degrees := RegularSupportData(blocks);

    if #poly_degrees eq 0 then
        if #supports eq 1 then
            ds := supports[1];
            if #ds eq 1 then
                if ds[1] eq 1 then
                    return "\\SL_2(q)";
                end if;
                return Sprintf("\\SL_2(\\Fq[t]/(t^%o))", ds[1]);
            elif SeqEqual(ds, [1, 1]) then
                return "\\Sp_4(q)";
            elif SeqEqual(ds, [2, 1]) then
                return "U_{2,1}(q)\\rtimes (\\SL_2(q)\\times \\SL_2(q))";
            end if;
        elif MatchesTwoSupports(supports, [1], [1]) then
            return "\\SL_2(q)\\times \\SL_2(q)";
        elif MatchesTwoSupports(supports, [2], [1]) then
            return "\\SL_2(\\Fq[t]/(t^2))\\times \\SL_2(q)";
        elif MatchesTwoSupports(supports, [1, 1], [1]) then
            return "\\Sp_4(q)\\times \\SL_2(q)";
        elif #supports eq 3 and &and[SeqEqual(ds, [1]) : ds in supports] then
            return "\\SL_2(q)\\times \\SL_2(q)\\times \\SL_2(q)";
        end if;
    elif #poly_degrees eq 1 and poly_degrees[1] eq 2 then
        if #supports eq 0 then
            return "\\SL_2(q^2)";
        elif #supports eq 1 and SeqEqual(supports[1], [1]) then
            return "\\SL_2(q)\\times \\SL_2(q^2)";
        end if;
    elif #poly_degrees eq 1 and poly_degrees[1] eq 3 and #supports eq 0 then
        return "\\SL_2(q^3)";
    end if;

    error "Unsupported regular core in stabilizer-table generator";
end function;

KernelStringNoRadicals := function(blocks)
    if &and[block`kind eq "S" : block in blocks] then
        if #blocks eq 1 then
            d := blocks[1]`d;
            return FqPowerString(2 * d - 2) cat "\\rtimes \\Fq^\\times";
        elif #blocks eq 2 and blocks[1]`d eq 2 and blocks[2]`d eq 2 then
            return "(\\Sym_2(\\Fq)\\oplus \\Sym_2(\\Fq))\\rtimes \\GL_2(q)";
        end if;
        error "Unsupported singular core in stabilizer-table generator";
    end if;

    nontrivial_s := [block : block in blocks | block`kind eq "S" and block`d gt 1];
    regular_core := [block : block in blocks | block`kind ne "S"];

    if #nontrivial_s eq 0 then
        return RegularKernelString(blocks);
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 2 and #regular_core gt 0 then
        return "U_L(q)\\rtimes (\\Fq^\\times\\times (" cat RegularKernelString(regular_core) cat "))";
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 3 and #regular_core eq 1 then
        return "U_L(q)\\rtimes (\\Fq^\\times\\times \\SL_2(q))";
    end if;

    error "Unsupported non-radical finite-field family in stabilizer-table generator";
end function;

KernelStringFromBlocks := function(blocks)
    r, core := SplitOffS1Blocks(blocks);
    if r gt 0 then
        m := TotalDimension(core);
        return "$" cat FqPowerString(r * m) cat "\\rtimes (\\GL_" cat Sprint(r)
            cat "(q)\\times (" cat KernelStringNoRadicals(core) cat "))$";
    end if;

    return "$" cat KernelStringNoRadicals(blocks) cat "$";
end function;

PrintFiniteStabilizerTable := procedure(n)
    types := FiniteOrbitTypesNLe7(n);

    print "";
    print "\\paragraph{$n=" cat Sprint(n) cat "$.}";
    print "\\begin{longtable}{@{}c p{0.28\\textwidth} p{0.47\\textwidth} p{0.12\\textwidth}@{}}";
    print "\\toprule";
    print "No. & Finite-field type & $K_L$ in $\\Stab_{\\GL_n(q)}(L)=K_L\\rtimes Q_L$ & $Q_L$ \\\\";
    print "\\midrule";
    print "\\endfirsthead";
    print "\\toprule";
    print "No. & Finite-field type & $K_L$ in $\\Stab_{\\GL_n(q)}(L)=K_L\\rtimes Q_L$ & $Q_L$ \\\\";
    print "\\midrule";
    print "\\endhead";
    print "\\midrule";
    print "\\multicolumn{4}{r}{Continued on next page}\\\\";
    print "\\midrule";
    print "\\endfoot";
    print "\\bottomrule";
    print "\\endlastfoot";

    for i in [1..#types] do
        kernel := KernelStringFromBlocks(types[i]`blocks);
        quotient := QuotientStringFromBlocks(types[i]`blocks);
        print Sprint(i) cat " & " cat LatexTypeString(types[i]`name)
            cat " & " cat kernel cat " & " cat quotient cat " \\\\";
    end for;

    print "\\end{longtable}";
end procedure;

print "% This file is generated by magma/examples/generate_geometric_finite_field_stabilizer_tables_tex.m.";
for n in [4..7] do
    PrintFiniteStabilizerTable(n);
end for;

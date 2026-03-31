SetColumns(0);
load "src/paper_data.m";

LatexTypeString := function(s)
    return "$" cat s cat "$";
end function;

QuotientLatex := function(s)
    if s eq "GL_2(k)" then
        return "$\\GL_2(k)$";
    elif s eq "B" then
        return "$B$";
    elif s eq "N_T" then
        return "$N_T$";
    elif s eq "T" then
        return "$T$";
    elif s eq "preimage of S_3" then
        return "$S_3$ preimage";
    end if;

    error "Unsupported geometric quotient string";
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
    end if;

    error "Unknown geometric block kind";
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

    for block in blocks do
        if block`kind eq "J" then
            if block`label eq "a" then
                Append(~a_support, block`d);
            elif block`label eq "b" then
                Append(~b_support, block`d);
            elif block`label eq "c" then
                Append(~c_support, block`d);
            end if;
        end if;
    end for;

    Sort(~a_support);
    Reverse(~a_support);
    Sort(~b_support);
    Reverse(~b_support);
    Sort(~c_support);
    Reverse(~c_support);

    return [seq : seq in [a_support, b_support, c_support] | #seq gt 0];
end function;

GaPowerString := function(m)
    return Sprintf("\\Ga^{%o}", m);
end function;

RegularKernelString := function(blocks)
    supports := RegularSupportData(blocks);

    if #supports eq 1 then
        ds := supports[1];
        if #ds eq 1 then
            if ds[1] eq 1 then
                return "\\SL_2(k)";
            end if;
            return Sprintf("\\SL_2(k[t]/(t^%o))", ds[1]);
        elif SeqEqual(ds, [1, 1]) then
            return "\\Sp_4(k)";
        elif SeqEqual(ds, [2, 1]) then
            return "U_{2,1}(k)\\rtimes (\\SL_2(k)\\times \\SL_2(k))";
        end if;
    elif MatchesTwoSupports(supports, [1], [1]) then
        return "\\SL_2(k)\\times \\SL_2(k)";
    elif MatchesTwoSupports(supports, [2], [1]) then
        return "\\SL_2(k[t]/(t^2))\\times \\SL_2(k)";
    elif MatchesTwoSupports(supports, [1, 1], [1]) then
        return "\\Sp_4(k)\\times \\SL_2(k)";
    elif #supports eq 3 and &and[SeqEqual(ds, [1]) : ds in supports] then
        return "\\SL_2(k)\\times \\SL_2(k)\\times \\SL_2(k)";
    end if;

    error "Unsupported geometric regular core";
end function;

KernelStringNoRadicals := function(blocks)
    if &and[block`kind eq "S" : block in blocks] then
        if #blocks eq 1 then
            d := blocks[1]`d;
            return GaPowerString(2 * d - 2) cat "\\rtimes \\Gm";
        elif #blocks eq 2 and blocks[1]`d eq 2 and blocks[2]`d eq 2 then
            return "(\\Sym_2(k)\\oplus \\Sym_2(k))\\rtimes \\GL_2(k)";
        end if;

        error "Unsupported geometric singular core";
    end if;

    nontrivial_s := [block : block in blocks | block`kind eq "S" and block`d gt 1];
    regular_core := [block : block in blocks | block`kind ne "S"];

    if #nontrivial_s eq 0 then
        return RegularKernelString(blocks);
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 2 and #regular_core gt 0 then
        return "U_L\\rtimes (\\Gm\\times (" cat RegularKernelString(regular_core) cat "))";
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 3 and #regular_core eq 1 then
        return "U_L\\rtimes (\\Gm\\times \\SL_2(k))";
    end if;

    error "Unsupported geometric non-radical family";
end function;

KernelLatexFromBlocks := function(blocks)
    r, core := SplitOffS1Blocks(blocks);
    if r gt 0 then
        m := TotalDimension(core);
        return "$" cat GaPowerString(r * m) cat "\\rtimes (\\GL_" cat Sprint(r)
            cat "(k)\\times (" cat KernelStringNoRadicals(core) cat "))$";
    end if;
    return "$" cat KernelStringNoRadicals(blocks) cat "$";
end function;

PrintGeometricStabilizerTable := procedure(n)
    types := GeometricOrbitTypesNLe7(n);
    print "";
    print "\\paragraph{$n=" cat Sprint(n) cat "$.}";
    print "\\begin{longtable}{@{}c p{0.28\\textwidth} p{0.47\\textwidth} p{0.12\\textwidth}@{}}";
    print "\\toprule";
    print "No. & Geometric type & $K_L$ in $\\Stab_{\\GL(V)}(L)=K_L\\rtimes Q_L$ & $Q_L$ \\\\";
    print "\\midrule";
    print "\\endfirsthead";
    print "\\toprule";
    print "No. & Geometric type & $K_L$ in $\\Stab_{\\GL(V)}(L)=K_L\\rtimes Q_L$ & $Q_L$ \\\\";
    print "\\midrule";
    print "\\endhead";
    print "\\midrule";
    print "\\multicolumn{4}{r}{Continued on next page}\\\\";
    print "\\midrule";
    print "\\endfoot";
    print "\\bottomrule";
    print "\\endlastfoot";

    for i in [1..#types] do
        kernel := KernelLatexFromBlocks(types[i]`blocks);
        quotient := QuotientLatex(GeometricQuotientStringFromBlocks(types[i]`blocks));
        print Sprint(i) cat " & " cat LatexTypeString(types[i]`name)
            cat " & " cat kernel cat " & " cat quotient cat " \\\\";
    end for;

    print "\\end{longtable}";
end procedure;

print "% This file is generated by magma/examples/generate_geometric_stabilizer_tables_tex.m.";
for n in [4..7] do
    PrintGeometricStabilizerTable(n);
end for;

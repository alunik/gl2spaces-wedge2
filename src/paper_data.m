BlockSpecRF := recformat< kind, d, label, degree >;
PaperOrbitRF := recformat< name, blocks >;
PaperOrbitDataRF := recformat< name, blocks, pair, wedge1, wedge2, kernel, quotient >;

JoinStrings := function(parts, sep)
    if #parts eq 0 then
        return "";
    end if;

    out := parts[1];
    for i in [2..#parts] do
        out cat:= sep cat parts[i];
    end for;
    return out;
end function;

SBlock := function(d)
    return rec<BlockSpecRF | kind := "S", d := d, label := "", degree := 0>;
end function;

JBlock := function(label, d)
    return rec<BlockSpecRF | kind := "J", d := d, label := label, degree := 0>;
end function;

QBlock := function()
    return rec<BlockSpecRF | kind := "P", d := 0, label := "", degree := 2>;
end function;

CBlock := function()
    return rec<BlockSpecRF | kind := "P", d := 0, label := "", degree := 3>;
end function;

CrossBlockPair := function(X, Y)
    K := BaseRing(X);
    r := Nrows(X);
    c := Ncols(X);

    if not (Nrows(Y) eq r and Ncols(Y) eq c) then
        error "Cross blocks must have the same shape";
    end if;

    n := r + c;
    A := ZeroMatrix(K, n, n);
    B := ZeroMatrix(K, n, n);

    for i in [1..r] do
        for j in [1..c] do
            if X[i, j] ne 0 then
                A[i, r + j] := X[i, j];
                A[r + j, i] := -X[i, j];
            end if;

            if Y[i, j] ne 0 then
                B[i, r + j] := Y[i, j];
                B[r + j, i] := -Y[i, j];
            end if;
        end for;
    end for;

    return <A, B>;
end function;

SingularBlockPair := function(F, d)
    if d eq 1 then
        return <ZeroMatrix(F, 1, 1), ZeroMatrix(F, 1, 1)>;
    end if;

    X := ZeroMatrix(F, d, d - 1);
    Y := ZeroMatrix(F, d, d - 1);
    for i in [1..d - 1] do
        X[i, i] := F!1;
        Y[i + 1, i] := F!1;
    end for;

    return CrossBlockPair(X, Y);
end function;

JordanBlockMatrix := function(F, d, lambda)
    J := ZeroMatrix(F, d, d);
    for i in [1..d] do
        J[i, i] := lambda;
    end for;
    for i in [1..d - 1] do
        J[i, i + 1] := F!1;
    end for;
    return J;
end function;

RegularRationalBlockPair := function(F, label, d)
    if label eq "a" then
        X := IdentityMatrix(F, d);
        Y := JordanBlockMatrix(F, d, F!0);
    elif label eq "b" then
        X := IdentityMatrix(F, d);
        Y := JordanBlockMatrix(F, d, F!1);
    elif label eq "c" then
        X := JordanBlockMatrix(F, d, F!0);
        Y := IdentityMatrix(F, d);
    else
        error "Rational labels must be a, b, c";
    end if;

    return CrossBlockPair(X, Y);
end function;

RegularPolynomialBlockPair := function(F, degree)
    g := IrreduciblePolynomial(F, degree);
    C := CompanionMatrix(g);
    return CrossBlockPair(IdentityMatrix(F, degree), C);
end function;

BlockPairFromSpec := function(F, spec)
    if spec`kind eq "S" then
        return SingularBlockPair(F, spec`d);
    elif spec`kind eq "J" then
        return RegularRationalBlockPair(F, spec`label, spec`d);
    elif spec`kind eq "P" then
        return RegularPolynomialBlockPair(F, spec`degree);
    end if;

    error Sprintf("Unknown block kind %o", spec`kind);
end function;

OrthogonalDirectSumPairs := function(pairs)
    F := BaseRing(pairs[1][1]);
    dims := [Nrows(pair[1]) : pair in pairs];
    n := &+dims;
    A := ZeroMatrix(F, n, n);
    B := ZeroMatrix(F, n, n);

    offset := 0;
    for pair in pairs do
        d := Nrows(pair[1]);
        InsertBlock(~A, pair[1], offset + 1, offset + 1);
        InsertBlock(~B, pair[2], offset + 1, offset + 1);
        offset +:= d;
    end for;

    return <A, B>;
end function;

RepresentativePairFromBlocks := function(F, blocks)
    pairs := [* *];
    for block in blocks do
        Append(~pairs, BlockPairFromSpec(F, block));
    end for;
    return OrthogonalDirectSumPairs(pairs);
end function;

WedgeTermString := function(i, j, coeff)
    term := Sprintf("e%o^e%o", i, j);
    if coeff eq 1 then
        return term;
    elif coeff eq -1 then
        return "-" cat term;
    end if;
    return "(" cat Sprint(coeff) cat ")*" cat term;
end function;

WedgeMatrixToString := function(M)
    terms := [];
    for i in [1..Nrows(M)] do
        for j in [i + 1..Ncols(M)] do
            if M[i, j] ne 0 then
                Append(~terms, WedgeTermString(i, j, M[i, j]));
            end if;
        end for;
    end for;

    if #terms eq 0 then
        return "0";
    end if;
    return JoinStrings(terms, " + ");
end function;

FiniteOrbitTypesNLe7 := function(n)
    if n eq 4 then
        return [
            rec<PaperOrbitRF | name := "J_{lambda,2}", blocks := [JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "J_{lambda,1} + J_{mu,1} (lambda != mu)", blocks := [JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "Q", blocks := [QBlock()]>,
            rec<PaperOrbitRF | name := "S_1 + S_2", blocks := [SBlock(1), SBlock(2)]>
        ];
    elif n eq 5 then
        return [
            rec<PaperOrbitRF | name := "S_3", blocks := [SBlock(3)]>,
            rec<PaperOrbitRF | name := "S_2 + S_1^2", blocks := [SBlock(2), SBlock(1), SBlock(1)]>,
            rec<PaperOrbitRF | name := "S_2 + J_{lambda,1}", blocks := [SBlock(2), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{lambda,2}", blocks := [SBlock(1), JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{lambda,1} + J_{mu,1} (lambda != mu)", blocks := [SBlock(1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + Q", blocks := [SBlock(1), QBlock()]>
        ];
    elif n eq 6 then
        return [
            rec<PaperOrbitRF | name := "J_{lambda,3}", blocks := [JBlock("a", 3)]>,
            rec<PaperOrbitRF | name := "J_{lambda,2} + J_{lambda,1}", blocks := [JBlock("a", 2), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "J_{lambda,2} + J_{mu,1} (lambda != mu)", blocks := [JBlock("a", 2), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "2J_{lambda,1} + J_{mu,1} (lambda != mu)", blocks := [JBlock("a", 1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "J_{lambda,1} + J_{mu,1} + J_{nu,1} (distinct)", blocks := [JBlock("a", 1), JBlock("b", 1), JBlock("c", 1)]>,
            rec<PaperOrbitRF | name := "J_{lambda,1} + Q", blocks := [JBlock("a", 1), QBlock()]>,
            rec<PaperOrbitRF | name := "C", blocks := [CBlock()]>,
            rec<PaperOrbitRF | name := "S_1^2 + J_{lambda,2}", blocks := [SBlock(1), SBlock(1), JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "S_1^2 + J_{lambda,1} + J_{mu,1} (lambda != mu)", blocks := [SBlock(1), SBlock(1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1^2 + Q", blocks := [SBlock(1), SBlock(1), QBlock()]>,
            rec<PaperOrbitRF | name := "S_1 + S_2 + J_{lambda,1}", blocks := [SBlock(1), SBlock(2), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_1^3 + S_2", blocks := [SBlock(1), SBlock(1), SBlock(1), SBlock(2)]>,
            rec<PaperOrbitRF | name := "S_1 + S_3", blocks := [SBlock(1), SBlock(3)]>,
            rec<PaperOrbitRF | name := "S_2^2", blocks := [SBlock(2), SBlock(2)]>
        ];
    elif n eq 7 then
        return [
            rec<PaperOrbitRF | name := "S_4", blocks := [SBlock(4)]>,
            rec<PaperOrbitRF | name := "S_3 + S_1^2", blocks := [SBlock(3), SBlock(1), SBlock(1)]>,
            rec<PaperOrbitRF | name := "S_2^2 + S_1", blocks := [SBlock(2), SBlock(2), SBlock(1)]>,
            rec<PaperOrbitRF | name := "S_2 + S_1^4", blocks := [SBlock(2), SBlock(1), SBlock(1), SBlock(1), SBlock(1)]>,
            rec<PaperOrbitRF | name := "S_3 + J_{lambda,1}", blocks := [SBlock(3), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_2 + S_1^2 + J_{lambda,1}", blocks := [SBlock(2), SBlock(1), SBlock(1), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_2 + J_{lambda,2}", blocks := [SBlock(2), JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "S_2 + 2J_{lambda,1}", blocks := [SBlock(2), JBlock("a", 1), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_2 + J_{lambda,1} + J_{mu,1} (lambda != mu)", blocks := [SBlock(2), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_2 + Q", blocks := [SBlock(2), QBlock()]>,
            rec<PaperOrbitRF | name := "S_1^3 + J_{lambda,2}", blocks := [SBlock(1), SBlock(1), SBlock(1), JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "S_1^3 + J_{lambda,1} + J_{mu,1} (lambda != mu)", blocks := [SBlock(1), SBlock(1), SBlock(1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1^3 + Q", blocks := [SBlock(1), SBlock(1), SBlock(1), QBlock()]>,
            rec<PaperOrbitRF | name := "S_1 + J_{lambda,3}", blocks := [SBlock(1), JBlock("a", 3)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{lambda,2} + J_{lambda,1}", blocks := [SBlock(1), JBlock("a", 2), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{lambda,2} + J_{mu,1} (lambda != mu)", blocks := [SBlock(1), JBlock("a", 2), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + 2J_{lambda,1} + J_{mu,1} (lambda != mu)", blocks := [SBlock(1), JBlock("a", 1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{lambda,1} + J_{mu,1} + J_{nu,1} (distinct)", blocks := [SBlock(1), JBlock("a", 1), JBlock("b", 1), JBlock("c", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{lambda,1} + Q", blocks := [SBlock(1), JBlock("a", 1), QBlock()]>,
            rec<PaperOrbitRF | name := "S_1 + C", blocks := [SBlock(1), CBlock()]>
        ];
    end if;

    error "This package only covers n = 4,5,6,7";
end function;

GeometricOrbitTypesNLe7 := function(n)
    if n eq 4 then
        return [
            rec<PaperOrbitRF | name := "J_{a,2}", blocks := [JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "J_{a,1} + J_{b,1}", blocks := [JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + S_2", blocks := [SBlock(1), SBlock(2)]>
        ];
    elif n eq 5 then
        return [
            rec<PaperOrbitRF | name := "S_3", blocks := [SBlock(3)]>,
            rec<PaperOrbitRF | name := "S_2 + S_1^2", blocks := [SBlock(2), SBlock(1), SBlock(1)]>,
            rec<PaperOrbitRF | name := "S_2 + J_{a,1}", blocks := [SBlock(2), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{a,2}", blocks := [SBlock(1), JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{a,1} + J_{b,1}", blocks := [SBlock(1), JBlock("a", 1), JBlock("b", 1)]>
        ];
    elif n eq 6 then
        return [
            rec<PaperOrbitRF | name := "J_{a,3}", blocks := [JBlock("a", 3)]>,
            rec<PaperOrbitRF | name := "J_{a,2} + J_{a,1}", blocks := [JBlock("a", 2), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "J_{a,2} + J_{b,1}", blocks := [JBlock("a", 2), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "2J_{a,1} + J_{b,1}", blocks := [JBlock("a", 1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "J_{a,1} + J_{b,1} + J_{c,1}", blocks := [JBlock("a", 1), JBlock("b", 1), JBlock("c", 1)]>,
            rec<PaperOrbitRF | name := "S_1^2 + J_{a,2}", blocks := [SBlock(1), SBlock(1), JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "S_1^2 + J_{a,1} + J_{b,1}", blocks := [SBlock(1), SBlock(1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + S_2 + J_{a,1}", blocks := [SBlock(1), SBlock(2), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_1^3 + S_2", blocks := [SBlock(1), SBlock(1), SBlock(1), SBlock(2)]>,
            rec<PaperOrbitRF | name := "S_1 + S_3", blocks := [SBlock(1), SBlock(3)]>,
            rec<PaperOrbitRF | name := "S_2^2", blocks := [SBlock(2), SBlock(2)]>
        ];
    elif n eq 7 then
        return [
            rec<PaperOrbitRF | name := "S_4", blocks := [SBlock(4)]>,
            rec<PaperOrbitRF | name := "S_3 + S_1^2", blocks := [SBlock(3), SBlock(1), SBlock(1)]>,
            rec<PaperOrbitRF | name := "S_2^2 + S_1", blocks := [SBlock(2), SBlock(2), SBlock(1)]>,
            rec<PaperOrbitRF | name := "S_2 + S_1^4", blocks := [SBlock(2), SBlock(1), SBlock(1), SBlock(1), SBlock(1)]>,
            rec<PaperOrbitRF | name := "S_3 + J_{a,1}", blocks := [SBlock(3), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_2 + S_1^2 + J_{a,1}", blocks := [SBlock(2), SBlock(1), SBlock(1), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_2 + J_{a,2}", blocks := [SBlock(2), JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "S_2 + 2J_{a,1}", blocks := [SBlock(2), JBlock("a", 1), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_2 + J_{a,1} + J_{b,1}", blocks := [SBlock(2), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1^3 + J_{a,2}", blocks := [SBlock(1), SBlock(1), SBlock(1), JBlock("a", 2)]>,
            rec<PaperOrbitRF | name := "S_1^3 + J_{a,1} + J_{b,1}", blocks := [SBlock(1), SBlock(1), SBlock(1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{a,3}", blocks := [SBlock(1), JBlock("a", 3)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{a,2} + J_{a,1}", blocks := [SBlock(1), JBlock("a", 2), JBlock("a", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{a,2} + J_{b,1}", blocks := [SBlock(1), JBlock("a", 2), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + 2J_{a,1} + J_{b,1}", blocks := [SBlock(1), JBlock("a", 1), JBlock("a", 1), JBlock("b", 1)]>,
            rec<PaperOrbitRF | name := "S_1 + J_{a,1} + J_{b,1} + J_{c,1}", blocks := [SBlock(1), JBlock("a", 1), JBlock("b", 1), JBlock("c", 1)]>
        ];
    end if;

    error "This package only covers n = 4,5,6,7";
end function;

SeqEqual := function(a, b)
    return #a eq #b and &and[a[i] eq b[i] : i in [1..#a]];
end function;

MatchesTwoSupports := function(supports, s1, s2)
    return (#supports eq 2) and
        ((SeqEqual(supports[1], s1) and SeqEqual(supports[2], s2)) or
         (SeqEqual(supports[1], s2) and SeqEqual(supports[2], s1)));
end function;

RegularSupportDataFinite := function(blocks)
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

    Sort(~a_support); Reverse(~a_support);
    Sort(~b_support); Reverse(~b_support);
    Sort(~c_support); Reverse(~c_support);
    return [seq : seq in [a_support, b_support, c_support] | #seq gt 0], poly_degrees;
end function;

RegularSupportDataGeometric := function(blocks)
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

    Sort(~a_support); Reverse(~a_support);
    Sort(~b_support); Reverse(~b_support);
    Sort(~c_support); Reverse(~c_support);
    return [seq : seq in [a_support, b_support, c_support] | #seq gt 0];
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

FiniteQuotientStringFromBlocks := function(blocks)
    supports, poly_degrees := RegularSupportDataFinite([block : block in blocks | block`kind ne "S"]);

    if #supports eq 0 and #poly_degrees eq 0 then
        return "GL_2(q)";
    end if;
    if #poly_degrees eq 0 then
        if #supports eq 1 then
            return "B(q)";
        elif #supports eq 2 then
            if MatchesTwoSupports(supports, [1], [1]) then
                return "N_s(q)";
            end if;
            return "T_s(q)";
        elif #supports eq 3 then
            return "~S_3(q)";
        end if;
    end if;
    if #poly_degrees eq 1 and poly_degrees[1] eq 2 then
        if #supports eq 0 then
            return "N_ns(q)";
        elif #supports eq 1 then
            return "~C_2(q)";
        end if;
    end if;
    if #poly_degrees eq 1 and poly_degrees[1] eq 3 and #supports eq 0 then
        return "~C_3(q)";
    end if;

    error "Unsupported finite-field quotient pattern";
end function;

GeometricQuotientStringFromBlocks := function(blocks)
    supports := RegularSupportDataGeometric([block : block in blocks | block`kind ne "S"]);
    if #supports eq 0 then
        return "GL_2(k)";
    elif #supports eq 1 then
        return "B";
    elif #supports eq 2 then
        if MatchesTwoSupports(supports, [1], [1]) then
            return "N_T";
        end if;
        return "T";
    elif #supports eq 3 then
        return "preimage of S_3";
    end if;

    error "Unsupported geometric quotient pattern";
end function;

FqPowerString := function(m)
    return Sprintf("F_q^%o", m);
end function;

GaPowerString := function(m)
    return Sprintf("G_a^%o", m);
end function;

FiniteRegularKernelString := function(blocks)
    supports, poly_degrees := RegularSupportDataFinite(blocks);

    if #poly_degrees eq 0 then
        if #supports eq 1 then
            ds := supports[1];
            if #ds eq 1 then
                if ds[1] eq 1 then
                    return "SL_2(q)";
                end if;
                return Sprintf("SL_2(F_q[t]/(t^%o))", ds[1]);
            elif SeqEqual(ds, [1, 1]) then
                return "Sp_4(q)";
            elif SeqEqual(ds, [2, 1]) then
                return "U_{2,1}(q) rtimes (SL_2(q) x SL_2(q))";
            end if;
        elif MatchesTwoSupports(supports, [1], [1]) then
            return "SL_2(q) x SL_2(q)";
        elif MatchesTwoSupports(supports, [2], [1]) then
            return "SL_2(F_q[t]/(t^2)) x SL_2(q)";
        elif MatchesTwoSupports(supports, [1, 1], [1]) then
            return "Sp_4(q) x SL_2(q)";
        elif #supports eq 3 and &and[SeqEqual(ds, [1]) : ds in supports] then
            return "SL_2(q) x SL_2(q) x SL_2(q)";
        end if;
    elif #poly_degrees eq 1 and poly_degrees[1] eq 2 then
        if #supports eq 0 then
            return "SL_2(q^2)";
        elif #supports eq 1 and SeqEqual(supports[1], [1]) then
            return "SL_2(q) x SL_2(q^2)";
        end if;
    elif #poly_degrees eq 1 and poly_degrees[1] eq 3 and #supports eq 0 then
        return "SL_2(q^3)";
    end if;

    error "Unsupported finite-field regular core";
end function;

GeometricRegularKernelString := function(blocks)
    supports := RegularSupportDataGeometric(blocks);

    if #supports eq 1 then
        ds := supports[1];
        if #ds eq 1 then
            if ds[1] eq 1 then
                return "SL_2(k)";
            end if;
            return Sprintf("SL_2(k[t]/(t^%o))", ds[1]);
        elif SeqEqual(ds, [1, 1]) then
            return "Sp_4(k)";
        elif SeqEqual(ds, [2, 1]) then
            return "U_{2,1}(k) rtimes (SL_2(k) x SL_2(k))";
        end if;
    elif MatchesTwoSupports(supports, [1], [1]) then
        return "SL_2(k) x SL_2(k)";
    elif MatchesTwoSupports(supports, [2], [1]) then
        return "SL_2(k[t]/(t^2)) x SL_2(k)";
    elif MatchesTwoSupports(supports, [1, 1], [1]) then
        return "Sp_4(k) x SL_2(k)";
    elif #supports eq 3 and &and[SeqEqual(ds, [1]) : ds in supports] then
        return "SL_2(k) x SL_2(k) x SL_2(k)";
    end if;

    error "Unsupported geometric regular core";
end function;

FiniteKernelStringNoRadicals := function(blocks)
    if &and[block`kind eq "S" : block in blocks] then
        if #blocks eq 1 then
            return FqPowerString(2 * blocks[1]`d - 2) cat " rtimes F_q^x";
        elif #blocks eq 2 and blocks[1]`d eq 2 and blocks[2]`d eq 2 then
            return "(Sym_2(F_q) + Sym_2(F_q)) rtimes GL_2(q)";
        end if;
        error "Unsupported finite singular core";
    end if;

    nontrivial_s := [block : block in blocks | block`kind eq "S" and block`d gt 1];
    regular_core := [block : block in blocks | block`kind ne "S"];

    if #nontrivial_s eq 0 then
        return FiniteRegularKernelString(blocks);
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 2 and #regular_core gt 0 then
        return "U_L(q) rtimes (F_q^x x (" cat FiniteRegularKernelString(regular_core) cat "))";
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 3 and #regular_core eq 1 then
        return "U_L(q) rtimes (F_q^x x SL_2(q))";
    end if;

    error "Unsupported finite non-radical family";
end function;

GeometricKernelStringNoRadicals := function(blocks)
    if &and[block`kind eq "S" : block in blocks] then
        if #blocks eq 1 then
            return GaPowerString(2 * blocks[1]`d - 2) cat " rtimes G_m";
        elif #blocks eq 2 and blocks[1]`d eq 2 and blocks[2]`d eq 2 then
            return "(Sym_2(k) + Sym_2(k)) rtimes GL_2(k)";
        end if;
        error "Unsupported geometric singular core";
    end if;

    nontrivial_s := [block : block in blocks | block`kind eq "S" and block`d gt 1];
    regular_core := [block : block in blocks | block`kind ne "S"];

    if #nontrivial_s eq 0 then
        return GeometricRegularKernelString(blocks);
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 2 and #regular_core gt 0 then
        return "U_L rtimes (G_m x (" cat GeometricRegularKernelString(regular_core) cat "))";
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 3 and #regular_core eq 1 then
        return "U_L rtimes (G_m x SL_2(k))";
    end if;

    error "Unsupported geometric non-radical family";
end function;

FiniteKernelStringFromBlocks := function(blocks)
    r, core := SplitOffS1Blocks(blocks);
    if r gt 0 then
        return FqPowerString(r * TotalDimension(core)) cat " rtimes (GL_" cat Sprint(r) cat "(q) x (" cat FiniteKernelStringNoRadicals(core) cat "))";
    end if;
    return FiniteKernelStringNoRadicals(blocks);
end function;

GeometricKernelStringFromBlocks := function(blocks)
    r, core := SplitOffS1Blocks(blocks);
    if r gt 0 then
        return GaPowerString(r * TotalDimension(core)) cat " rtimes (GL_" cat Sprint(r) cat "(k) x (" cat GeometricKernelStringNoRadicals(core) cat "))";
    end if;
    return GeometricKernelStringNoRadicals(blocks);
end function;

FiniteOrbitDataNLe7 := function(F, n)
    types := FiniteOrbitTypesNLe7(n);
    out := [* *];
    for T in types do
        pair := RepresentativePairFromBlocks(F, T`blocks);
        Append(~out, rec<PaperOrbitDataRF |
            name := T`name,
            blocks := T`blocks,
            pair := pair,
            wedge1 := WedgeMatrixToString(pair[1]),
            wedge2 := WedgeMatrixToString(pair[2]),
            kernel := FiniteKernelStringFromBlocks(T`blocks),
            quotient := FiniteQuotientStringFromBlocks(T`blocks)
        >);
    end for;
    return out;
end function;

GeometricOrbitDataNLe7 := function(F, n)
    types := GeometricOrbitTypesNLe7(n);
    out := [* *];
    for T in types do
        pair := RepresentativePairFromBlocks(F, T`blocks);
        Append(~out, rec<PaperOrbitDataRF |
            name := T`name,
            blocks := T`blocks,
            pair := pair,
            wedge1 := WedgeMatrixToString(pair[1]),
            wedge2 := WedgeMatrixToString(pair[2]),
            kernel := GeometricKernelStringFromBlocks(T`blocks),
            quotient := GeometricQuotientStringFromBlocks(T`blocks)
        >);
    end for;
    return out;
end function;

PrintPaperOrbitData := procedure(data : Title := "")
    if Title ne "" then
        print Title;
    end if;
    for i in [1..#data] do
        reci := data[i];
        printf "[%o] %o\n", i, reci`name;
        printf "  omega_1 = %o\n", reci`wedge1;
        printf "  omega_2 = %o\n", reci`wedge2;
        printf "  K_L = %o\n", reci`kernel;
        printf "  Q_L = %o\n", reci`quotient;
    end for;
end procedure;

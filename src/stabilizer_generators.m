load "src/paper_data.m";
load "src/local_matrix_models.m";

FiniteStabilizerSpecRF := recformat<
    name,
    blocks,
    pair,
    wedge1,
    wedge2,
    kernel,
    quotient,
    kernel_generators,
    quotient_generators,
    quotient_lifts,
    full_generators
>;

forward J21LocalKernelGenerators;
forward RegularClusterKernelGenerators;
forward RegularClustersWithOffsets;
forward DirectSumKernelGenerators;
forward S2ExtensionKernelGeneratorsFromBlocks;
forward RadicalExtensionKernelGenerators;
forward KernelGeneratorsNoRadicals;
forward KernelGeneratorsFromBlocks;
forward KernelOrderNoRadicals;
forward KernelOrderFromBlocks;
forward FiniteStabilizerSpecsNLe7;
forward PrintFiniteStabilizerSpecs;
forward RadicalConjugationData;
forward ActionMatrixOnPair;
forward SymPower2x2Matrix;
forward RegularLocalQuotientLift;
forward SingularLocalQuotientLift;
forward LocalBlockQuotientLift;
forward DirectSumQuotientLiftFromMatrix;
forward NSplitSwapLift;
forward S3Swap01Lift;
forward S3Swap0InfLift;
forward MixedQuadraticInvolutionMatrix;
forward CubicCycleMatrix;
forward MobiusMatrixFromImages;
forward FullStabilizerGeneratorsFromBlocks;

BlockMat := function(P, Q, R, S, T, U, V, W, X)
    row1 := HorizontalJoin(P, HorizontalJoin(Q, R));
    row2 := HorizontalJoin(S, HorizontalJoin(T, U));
    row3 := HorizontalJoin(V, HorizontalJoin(W, X));
    return VerticalJoin(row1, VerticalJoin(row2, row3));
end function;

SupportKey := function(block)
    if block`kind eq "J" then
        return "R:" cat block`label;
    elif block`kind eq "P" then
        return Sprintf("P:%o", block`degree);
    end if;

    error "SupportKey only applies to regular blocks";
end function;

RadicalConjugationData := function(F, blocks)
    dims := [BlockDimension(block) : block in blocks];
    starts := [];
    pos := 1;
    for d in dims do
        Append(~starts, pos);
        pos +:= d;
    end for;

    s1_ids := [i : i in [1..#blocks] | blocks[i]`kind eq "S" and blocks[i]`d eq 1];
    core_ids := [i : i in [1..#blocks] | not (blocks[i]`kind eq "S" and blocks[i]`d eq 1)];
    new_ids := s1_ids cat core_ids;

    old_positions := [];
    reordered_blocks := [];
    for i in new_ids do
        for j in [starts[i]..starts[i] + dims[i] - 1] do
            Append(~old_positions, j);
        end for;
        Append(~reordered_blocks, blocks[i]);
    end for;

    n := &+dims;
    P := ZeroMatrix(F, n, n);
    for i in [1..n] do
        P[i, old_positions[i]] := F!1;
    end for;

    return P, reordered_blocks;
end function;

EmbedInIdentity := function(F, total_dim, offset, g)
    H := IdentityMatrix(F, total_dim);
    InsertBlock(~H, g, offset + 1, offset + 1);
    return H;
end function;

SimpleSL2Generators := function(F)
    gens := [
        Matrix(F, 2, 2, [1, 1, 0, 1]),
        Matrix(F, 2, 2, [1, 0, 1, 1])
    ];

    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, Matrix(F, 2, 2, [w, 0, 0, w^-1]));
    end if;

    return gens;
end function;

UpperHankelBasis := function(F, d)
    mats := [];
    for k in [1..d] do
        H := ZeroMatrix(F, d, d);
        for i in [1..d] do
            j := k + 1 - i;
            if j ge 1 and j le d then
                H[i, j] := F!1;
            end if;
        end for;
        Append(~mats, H);
    end for;
    return mats;
end function;

LowerHankelBasis := function(F, d)
    mats := [];
    for k in [1..d] do
        H := ZeroMatrix(F, d, d);
        for i in [1..d] do
            j := d + k - i;
            if j ge 1 and j le d then
                H[i, j] := F!1;
            end if;
        end for;
        Append(~mats, H);
    end for;
    return mats;
end function;

SolvePolynomialUpperBasis := function(F, C)
    m := Nrows(C);
    V := VectorSpace(F, m * m);
    eqs := [];

    for i in [1..m] do
        for j in [1..m] do
            coeffs := [F | 0 : t in [1..m * m]];
            coeffs[(i - 1) * m + j] +:= F!1;
            coeffs[(j - 1) * m + i] -:= F!1;
            Append(~eqs, V!coeffs);
        end for;
    end for;

    for i in [1..m] do
        for j in [1..m] do
            coeffs := [F | 0 : t in [1..m * m]];
            for k in [1..m] do
                coeffs[(k - 1) * m + j] +:= C[i, k];
                coeffs[(i - 1) * m + k] -:= C[j, k];
            end for;
            Append(~eqs, V!coeffs);
        end for;
    end for;

    U := Nullspace(Transpose(Matrix(F, #eqs, m * m, &cat[Eltseq(e) : e in eqs])));
    return [Matrix(F, m, m, Eltseq(b)) : b in Basis(U)];
end function;

SolvePolynomialLowerBasis := function(F, C)
    m := Nrows(C);
    V := VectorSpace(F, m * m);
    eqs := [];

    for i in [1..m] do
        for j in [1..m] do
            coeffs := [F | 0 : t in [1..m * m]];
            coeffs[(i - 1) * m + j] +:= F!1;
            coeffs[(j - 1) * m + i] -:= F!1;
            Append(~eqs, V!coeffs);
        end for;
    end for;

    for i in [1..m] do
        for j in [1..m] do
            coeffs := [F | 0 : t in [1..m * m]];
            for k in [1..m] do
                coeffs[(i - 1) * m + k] +:= C[k, j];
                coeffs[(k - 1) * m + j] -:= C[k, i];
            end for;
            Append(~eqs, V!coeffs);
        end for;
    end for;

    U := Nullspace(Transpose(Matrix(F, #eqs, m * m, &cat[Eltseq(e) : e in eqs])));
    return [Matrix(F, m, m, Eltseq(b)) : b in Basis(U)];
end function;

JdLocalKernelGenerators := function(F, d)
    I := IdentityMatrix(F, d);
    Z := ZeroMatrix(F, d, d);
    gens := [];

    for H in UpperHankelBasis(F, d) do
        Append(~gens, VerticalJoin(HorizontalJoin(I, H), HorizontalJoin(Z, I)));
    end for;
    for H in LowerHankelBasis(F, d) do
        Append(~gens, VerticalJoin(HorizontalJoin(I, Z), HorizontalJoin(H, I)));
    end for;
    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, VerticalJoin(HorizontalJoin(w * I, Z), HorizontalJoin(Z, w^-1 * I)));
    end if;

    return gens;
end function;

Sp4LocalKernelGenerators := function(F)
    I := IdentityMatrix(F, 4);
    gens := [];

    g := I; g[1, 2] := F!1; Append(~gens, g);
    g := I; g[2, 1] := F!1; Append(~gens, g);
    g := I; g[3, 4] := F!1; Append(~gens, g);
    g := I; g[4, 3] := F!1; Append(~gens, g);
    g := I; g[1, 3] := F!1; g[4, 2] := -F!1; Append(~gens, g);
    g := I; g[1, 4] := F!1; g[3, 2] := F!1; Append(~gens, g);

    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, DiagonalMatrix(F, [w, w^-1, F!1, F!1]));
        Append(~gens, DiagonalMatrix(F, [F!1, F!1, w, w^-1]));
    end if;

    return gens;
end function;

PolynomialLocalKernelGenerators := function(F, m)
    C := CompanionMatrix(IrreduciblePolynomial(F, m));
    I := IdentityMatrix(F, m);
    Z := ZeroMatrix(F, m, m);
    gens := [];

    for X in SolvePolynomialUpperBasis(F, C) do
        Append(~gens, VerticalJoin(HorizontalJoin(I, X), HorizontalJoin(Z, I)));
    end for;
    for X in SolvePolynomialLowerBasis(F, C) do
        Append(~gens, VerticalJoin(HorizontalJoin(I, Z), HorizontalJoin(X, I)));
    end for;
    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, VerticalJoin(HorizontalJoin(w * I, Z), HorizontalJoin(Z, w^-1 * I)));
    end if;

    return gens;
end function;

J21LocalKernelGenerators := function(F)
    I2 := IdentityMatrix(F, 2);
    Z2 := ZeroMatrix(F, 2, 2);
    gens := [];

    for H in UpperHankelBasis(F, 2) do
        Append(~gens, BlockMat(I2, H, Z2, Z2, I2, Z2, Z2, Z2, I2));
    end for;
    for H in LowerHankelBasis(F, 2) do
        Append(~gens, BlockMat(I2, Z2, Z2, H, I2, Z2, Z2, Z2, I2));
    end for;
    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, BlockMat(w * I2, Z2, Z2, Z2, w^-1 * I2, Z2, Z2, Z2, I2));
    end if;

    for h in SimpleSL2Generators(F) do
        Append(~gens, BlockMat(I2, Z2, Z2, Z2, I2, Z2, Z2, Z2, h));
    end for;

    Append(~gens, BlockMat(I2, Z2, Z2, Matrix(F, 2, 2, [0, 1, 1, 0]), I2,
        ZeroMatrix(F, 2, 2), ZeroMatrix(F, 2, 2), Z2, I2));
    Append(~gens, BlockMat(I2, Z2, Z2, Matrix(F, 2, 2, [0, 0, 0, 1]), I2,
        ZeroMatrix(F, 2, 2), ZeroMatrix(F, 2, 2), Z2, I2));
    Append(~gens, BlockMat(I2, Z2, Z2, ZeroMatrix(F, 2, 2), I2,
        Matrix(F, 2, 2, [0, 0, 1, 0]), Matrix(F, 2, 2, [0, 0, 0, 1]), Z2, I2));
    Append(~gens, BlockMat(I2, Z2, Z2, ZeroMatrix(F, 2, 2), I2,
        Matrix(F, 2, 2, [0, 0, 0, 1]), Matrix(F, 2, 2, [0, -1, 0, 0]), Z2, I2));

    return gens;
end function;

SingularBlockKernelGenerators := function(F, d)
    gens := [];
    for i in [1..2 * d - 2] do
        coeffs := [F | 0 : j in [1..2 * d - 2]];
        coeffs[i] := F!1;
        Append(~gens, SingularKernelCandidateMatrix(F, d, F!1, coeffs));
    end for;
    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, SingularKernelCandidateMatrix(F, d, w,
            [F | 0 : i in [1..2 * d - 2]]));
    end if;
    return gens;
end function;

S22LocalKernelGenerators := function(F)
    I2 := IdentityMatrix(F, 2);
    Z := ZeroMatrix(F, 2, 2);
    sym_basis := [
        Matrix(F, 2, 2, [1, 0, 0, 0]),
        Matrix(F, 2, 2, [0, 1, 1, 0]),
        Matrix(F, 2, 2, [0, 0, 0, 1])
    ];

    old_gens := [];
    for C in sym_basis do
        Append(~old_gens, S22CandidateMatrix(F, I2, C, Z));
        Append(~old_gens, S22CandidateMatrix(F, I2, Z, C));
    end for;
    Append(~old_gens, S22CandidateMatrix(F, Matrix(F, 2, 2, [1, 1, 0, 1]), Z, Z));
    Append(~old_gens, S22CandidateMatrix(F, Matrix(F, 2, 2, [0, -1, 1, 0]), Z, Z));
    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~old_gens, S22CandidateMatrix(F, Matrix(F, 2, 2, [w, 0, 0, 1]), Z, Z));
    end if;

    P := PermutationMatrix(F, Sym(6)!(2, 4, 5, 3));
    return [Transpose(P) * g * P : g in old_gens];
end function;

S3PlusJ1KernelGenerators := function(F)
    gens := [
        S3PlusJ1UnipotentMatrix(F, F!1, F!0, F!0, F!0, F!0, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!1, F!0, F!0, F!0, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!0, F!1, F!0, F!0, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!0, F!0, F!1, F!0, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!0, F!0, F!0, F!1, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!0, F!0, F!0, F!0, F!1),
        S3PlusJ1LeviMatrix(F, F!1, Matrix(F, 2, 2, [1, 1, 0, 1])),
        S3PlusJ1LeviMatrix(F, F!1, Matrix(F, 2, 2, [1, 0, -1, 1]))
    ];

    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, S3PlusJ1LeviMatrix(F, w, IdentityMatrix(F, 2)));
        Append(~gens, S3PlusJ1LeviMatrix(F, F!1, Matrix(F, 2, 2, [w, 0, 0, w^-1])));
    end if;

    return gens;
end function;

RegularClustersWithOffsets := function(blocks)
    if #blocks eq 0 then
        return [];
    end if;

    clusters := [];
    start := 1;
    current_blocks := [blocks[1]];
    current_key := SupportKey(blocks[1]);
    offset := 0;
    current_offset := 0;

    for i in [2..#blocks] do
        offset +:= BlockDimension(blocks[i - 1]);
        key := SupportKey(blocks[i]);
        if key eq current_key then
            Append(~current_blocks, blocks[i]);
        else
            Append(~clusters, <current_offset, current_blocks>);
            current_blocks := [blocks[i]];
            current_key := key;
            current_offset := offset;
        end if;
    end for;

    Append(~clusters, <current_offset, current_blocks>);
    return clusters;
end function;

RegularClusterKernelGenerators := function(F, cluster_blocks)
    if #cluster_blocks eq 0 then
        return [];
    end if;

    if cluster_blocks[1]`kind eq "P" then
        return PolynomialLocalKernelGenerators(F, cluster_blocks[1]`degree);
    end if;

    ds := [block`d : block in cluster_blocks];
    if #ds eq 1 then
        return JdLocalKernelGenerators(F, ds[1]);
    elif ds eq [1, 1] then
        return Sp4LocalKernelGenerators(F);
    elif ds eq [2, 1] then
        return J21LocalKernelGenerators(F);
    end if;

    error "Unsupported regular cluster";
end function;

DirectSumKernelGenerators := function(F, block_clusters)
    total_dim := &+[&+[BlockDimension(block) : block in entry[2]] : entry in block_clusters];
    gens := [];
    for entry in block_clusters do
        offset := entry[1];
        local_blocks := entry[2];
        local_gens := RegularClusterKernelGenerators(F, local_blocks);
        local_dim := TotalDimension(local_blocks);
        for g in local_gens do
            Append(~gens, EmbedInIdentity(F, total_dim, offset, g));
        end for;
    end for;
    return gens;
end function;

S2ExtensionElementRow := function(a, x, y, w0, h, A0, B0)
    F := BaseRing(A0);
    m := Nrows(A0);
    col1 := -a * h * A0 * Transpose(w0);
    col2 := -a * h * B0 * Transpose(w0);

    row1 := HorizontalJoin(Matrix(F, 1, 1, [a]),
        HorizontalJoin(Matrix(F, 1, 1, [0]),
        HorizontalJoin(Matrix(F, 1, 1, [0]), ZeroMatrix(F, 1, m))));
    row2 := HorizontalJoin(Matrix(F, 1, 1, [0]),
        HorizontalJoin(Matrix(F, 1, 1, [a]),
        HorizontalJoin(Matrix(F, 1, 1, [0]), ZeroMatrix(F, 1, m))));
    row3 := HorizontalJoin(Matrix(F, 1, 1, [x]),
        HorizontalJoin(Matrix(F, 1, 1, [y]),
        HorizontalJoin(Matrix(F, 1, 1, [a^-1]), w0)));

    bottom_left := HorizontalJoin(col1,
        HorizontalJoin(col2, ZeroMatrix(F, m, 1)));
    bottom := HorizontalJoin(bottom_left, h);

    return VerticalJoin(row1, VerticalJoin(row2, VerticalJoin(row3, bottom)));
end function;

S2ExtensionKernelGeneratorsFromBlocks := function(F, core_blocks)
    core_pair := RepresentativePairFromBlocks(F, core_blocks);
    A0 := core_pair[1];
    B0 := core_pair[2];
    m := Nrows(A0);
    core_gens := KernelGeneratorsNoRadicals(F, core_blocks);
    gens := [];

    Append(~gens, S2ExtensionElementRow(F!1, F!1, F!0,
        ZeroMatrix(F, 1, m), IdentityMatrix(F, m), A0, B0));
    Append(~gens, S2ExtensionElementRow(F!1, F!0, F!1,
        ZeroMatrix(F, 1, m), IdentityMatrix(F, m), A0, B0));

    for i in [1..m] do
        w0 := ZeroMatrix(F, 1, m);
        w0[1, i] := F!1;
        Append(~gens, S2ExtensionElementRow(F!1, F!0, F!0,
            w0, IdentityMatrix(F, m), A0, B0));
    end for;

    for h in core_gens do
        Append(~gens, S2ExtensionElementRow(F!1, F!0, F!0,
            ZeroMatrix(F, 1, m), h, A0, B0));
    end for;

    if #F gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, S2ExtensionElementRow(w, F!0, F!0,
            ZeroMatrix(F, 1, m), IdentityMatrix(F, m), A0, B0));
    end if;

    return gens;
end function;

RadicalExtensionKernelGenerators := function(F, r, core_blocks)
    core_pair := RepresentativePairFromBlocks(F, core_blocks);
    core_dim := Nrows(core_pair[1]);
    total_dim := r + core_dim;
    gens := [];

    for g in Generators(GL(r, F)) do
        H := IdentityMatrix(F, total_dim);
        InsertBlock(~H, g, 1, 1);
        Append(~gens, H);
    end for;

    for h in KernelGeneratorsNoRadicals(F, core_blocks) do
        H := IdentityMatrix(F, total_dim);
        InsertBlock(~H, h, r + 1, r + 1);
        Append(~gens, H);
    end for;

    for i in [1..core_dim] do
        for j in [1..r] do
            H := IdentityMatrix(F, total_dim);
            H[r + i, j] := F!1;
            Append(~gens, H);
        end for;
    end for;

    return gens;
end function;

KernelGeneratorsNoRadicals := function(F, blocks)
    if &and[block`kind eq "S" : block in blocks] then
        if #blocks eq 1 then
            return SingularBlockKernelGenerators(F, blocks[1]`d);
        elif #blocks eq 2 and blocks[1]`d eq 2 and blocks[2]`d eq 2 then
            return S22LocalKernelGenerators(F);
        end if;
        error "Unsupported pure singular core";
    end if;

    nontrivial_s := [block : block in blocks | block`kind eq "S" and block`d gt 1];
    regular_core := [block : block in blocks | block`kind ne "S"];

    if #nontrivial_s eq 0 then
        return DirectSumKernelGenerators(F, RegularClustersWithOffsets(blocks));
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 2 and #regular_core gt 0 then
        return S2ExtensionKernelGeneratorsFromBlocks(F, regular_core);
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 3 and #regular_core eq 1 and
        regular_core[1]`kind eq "J" and regular_core[1]`d eq 1 then
        return S3PlusJ1KernelGenerators(F);
    end if;

    error "Unsupported non-radical family";
end function;

KernelGeneratorsFromBlocks := function(F, blocks)
    r, core := SplitOffS1Blocks(blocks);
    if r gt 0 then
        P, reordered_blocks := RadicalConjugationData(F, blocks);
        _, reordered_core := SplitOffS1Blocks(reordered_blocks);
        gens := RadicalExtensionKernelGenerators(F, r, reordered_core);
        return [Transpose(P) * g * P : g in gens];
    end if;
    return KernelGeneratorsNoRadicals(F, blocks);
end function;

ExtensionMultiplicationMatrix := function(F, m, alpha)
    basis := [F!1];
    if m gt 1 then
        z := alpha;
        basis := [F!1] cat [z^i : i in [1..m - 1]];
    end if;

    cols := [];
    for b in basis do
        Append(~cols, Coordinates(F, alpha * b));
    end for;
    return Matrix(F, m, m, &cat cols);
end function;

PrimitiveExtensionTorusGenerator := function(F, m)
    E<z> := ext<F | m>;
    alpha := PrimitiveElement(E);
    basis := [E!1] cat [z^i : i in [1..m - 1]];
    cols := [];
    for b in basis do
        Append(~cols, Eltseq(E ! (alpha * b)));
    end for;
    return Matrix(F, m, m, &cat [Eltseq(v) : v in cols]);
end function;

ExtensionFrobeniusMatrix := function(F, m)
    E<z> := ext<F | m>;
    basis := [E!1] cat [z^i : i in [1..m - 1]];
    cols := [];
    for b in basis do
        Append(~cols, Eltseq(b^(#F)));
    end for;
    return Matrix(F, m, m, &cat [Eltseq(v) : v in cols]);
end function;

QuotientGeneratorsFromBlocks := function(F, blocks)
    qtype := FiniteQuotientStringFromBlocks(blocks);
    q := #F;
    w := F!1;
    if q gt 2 then
        w := PrimitiveElement(F);
    end if;

    U := Matrix(F, 2, 2, [1, 1, 0, 1]);
    D := Matrix(F, 2, 2, [w, 0, 0, 1]);
    C := Matrix(F, 2, 2, [w, 0, 0, w]);
    S := Matrix(F, 2, 2, [0, 1, 1, 0]);
    R3 := Matrix(F, 2, 2, [0, -1, 1, -1]);

    if qtype eq "GL_2(q)" then
        gens := [U, Matrix(F, 2, 2, [0, -1, 1, 0])];
        if q gt 2 then
            Append(~gens, D);
        end if;
        return gens;
    elif qtype eq "B(q)" then
        gens := [U];
        if q gt 2 then
            Append(~gens, D);
            Append(~gens, C);
        end if;
        return gens;
    elif qtype eq "T_s(q)" then
        if q eq 2 then
            return [IdentityMatrix(F, 2)];
        end if;
        return [D, C];
    elif qtype eq "N_s(q)" then
        gens := [S];
        if q gt 2 then
            Append(~gens, D);
            Append(~gens, C);
        end if;
        return gens;
    elif qtype eq "N_ns(q)" then
        T := PrimitiveExtensionTorusGenerator(F, 2);
        gens := [T, ExtensionFrobeniusMatrix(F, 2)];
        if q gt 2 then
            Append(~gens, C);
        end if;
        return gens;
    elif qtype eq "~C_2(q)" then
        gens := [S];
        if q gt 2 then
            Append(~gens, C);
        end if;
        return gens;
    elif qtype eq "~C_3(q)" then
        gens := [R3];
        if q gt 2 then
            Append(~gens, C);
        end if;
        return gens;
    elif qtype eq "~S_3(q)" then
        gens := [S, R3];
        if q gt 2 then
            Append(~gens, C);
        end if;
        return gens;
    end if;

    error "Unsupported quotient type";
end function;

ActionMatrixOnPair := function(pair, g)
    F := BaseRing(pair[1]);
    A := pair[1];
    B := pair[2];
    Ag := g * A * Transpose(g);
    Bg := g * B * Transpose(g);
    M := ZeroMatrix(F, 2, 2);

    for j in [1..2] do
        target := [Ag, Bg][j];
        found := false;
        for x in F do
            for y in F do
                if x * A + y * B eq target then
                    M[1, j] := x;
                    M[2, j] := y;
                    found := true;
                    break x;
                end if;
            end for;
        end for;
        assert found;
    end for;

    return M;
end function;

SymPower2x2Matrix := function(F, M, n)
    if n eq 0 then
        return Matrix(F, 1, 1, [F!1]);
    end if;

    mons := [<n - i, i> : i in [0..n]];
    P := ZeroMatrix(F, n + 1, n + 1);
    a := M[1, 1];
    b := M[1, 2];
    c := M[2, 1];
    d := M[2, 2];

    for j in [1..n + 1] do
        e1 := mons[j][1];
        e2 := mons[j][2];
        coeffs := [F | 0 : i in [1..n + 1]];
        for i1 in [0..e1] do
            coeff1 := Binomial(e1, i1) * a^(e1 - i1) * b^i1;
            for i2 in [0..e2] do
                coeff2 := Binomial(e2, i2) * c^(e2 - i2) * d^i2;
                xpow := (e1 - i1) + (e2 - i2);
                idx := n - xpow + 1;
                coeffs[idx] +:= coeff1 * coeff2;
            end for;
        end for;
        for i in [1..n + 1] do
            P[i, j] := coeffs[i];
        end for;
    end for;

    return P;
end function;

SingularLocalQuotientLift := function(F, d, M)
    if d eq 1 then
        return IdentityMatrix(F, 1);
    end if;

    P := SymPower2x2Matrix(F, M, d - 1);
    Q := Transpose(SymPower2x2Matrix(F, M, d - 2)^-1);
    return DiagonalJoin(P, Q);
end function;

RegularLocalQuotientLift := function(F, Y, M)
    m := Nrows(Y);
    I := IdentityMatrix(F, m);
    a := M[1, 1];
    b := M[1, 2];
    c := M[2, 1];
    d := M[2, 2];
    S := a * I + b * Y;
    if Determinant(S) eq 0 then
        error "Chosen codomain matrix does not preserve this regular support block";
    end if;

    Yp := (c * I + d * Y) * S^-1;
    flag, P := IsSimilar(Y, Yp);
    if not flag then
        error "Failed to construct a local quotient lift for a regular block";
    end if;

    Q := Transpose(P^-1 * S);
    return DiagonalJoin(P, Q);
end function;

CoordsOverBase := function(F, E, x)
    if Type(x) eq FldFinElt and Parent(x) cmpeq F then
        return [F | x];
    end if;

    coeffs := Eltseq(E ! x);
    return [F | c : c in coeffs];
end function;

MobiusMatrixFromImages := function(F, E, src, dst)
    RowsFromCoeffs := function(coeffs)
        seqs := [CoordsOverBase(F, E, coeff) : coeff in coeffs];
        m := #seqs[1];
        return [[F | seqs[j][i] : j in [1..#coeffs]] : i in [1..m]];
    end function;

    rows := [];

    for k in [1..#src] do
        s := src[k];
        t := dst[k];

        if Type(s) eq MonStgElt and s eq "inf" then
            if Type(t) eq MonStgElt and t eq "inf" then
                Append(~rows, [F | 0, 1, 0, 0]);
            else
                rows cat:= RowsFromCoeffs([E!0, -(E!t), E!0, E!1]);
            end if;
        else
            if Type(t) eq MonStgElt and t eq "inf" then
                rows cat:= RowsFromCoeffs([E!1, E!s, E!0, E!0]);
            else
                rows cat:= RowsFromCoeffs([-(E!t), -(E!t) * (E!s), E!1, E!s]);
            end if;
        end if;
    end for;

    Msys := Matrix(F, #rows, 4, &cat rows);
    N := Nullspace(Transpose(Msys));
    for v in Basis(N) do
        a0 := v[1];
        b0 := v[2];
        c0 := v[3];
        d0 := v[4];
        if a0 * d0 - b0 * c0 ne 0 then
            return Matrix(F, 2, 2, [a0, b0, c0, d0]);
        end if;
    end for;

    error "No invertible Möbius matrix found for the requested support action";
end function;

MixedQuadraticInvolutionMatrix := function(F)
    g := IrreduciblePolynomial(F, 2);
    E<a> := ext<F | g>;
    return MobiusMatrixFromImages(F, E, [E!0, a, a^#F], [E!0, a^#F, a]);
end function;

CubicCycleMatrix := function(F)
    g := IrreduciblePolynomial(F, 3);
    E<a> := ext<F | g>;
    q := #F;
    return MobiusMatrixFromImages(F, E, [a, a^q, a^(q^2)], [a^q, a^(q^2), a]);
end function;

LocalBlockQuotientLift := function(F, block, M)
    if block`kind eq "S" then
        return SingularLocalQuotientLift(F, block`d, M);
    elif block`kind eq "J" then
        if block`label eq "a" then
            Y := JordanBlockMatrix(F, block`d, F!0);
            return RegularLocalQuotientLift(F, Y, M);
        elif block`label eq "b" then
            if block`d ne 1 then
                error "Only the simple b-support block occurs in range";
            end if;
            Y := Matrix(F, 1, 1, [F!1]);
            return RegularLocalQuotientLift(F, Y, M);
        elif block`label eq "c" then
            if block`d ne 1 then
                error "Only the simple c-support block occurs in range";
            end if;
            if M[1, 2] ne 0 then
                error "The c-support block requires a permutation lift here";
            end if;
            return DiagonalMatrix(F, [F!1, M[2, 2]]);
        end if;
    elif block`kind eq "P" then
        Y := CompanionMatrix(IrreduciblePolynomial(F, block`degree));
        return RegularLocalQuotientLift(F, Y, M);
    end if;

    error "Unsupported local quotient lift";
end function;

DirectSumQuotientLiftFromMatrix := function(F, blocks, M)
    total_dim := TotalDimension(blocks);
    H := IdentityMatrix(F, total_dim);
    offset := 0;
    for block in blocks do
        d := BlockDimension(block);
        G := IdentityMatrix(F, d);
        if not (block`kind eq "S" and block`d eq 1) then
            G := LocalBlockQuotientLift(F, block, M);
        end if;
        InsertBlock(~H, G, offset + 1, offset + 1);
        offset +:= d;
    end for;
    return H;
end function;

BlockStartPositions := function(blocks)
    starts := [];
    pos := 1;
    for block in blocks do
        Append(~starts, pos);
        pos +:= BlockDimension(block);
    end for;
    return starts;
end function;

SwapTwoBlocksMatrix := function(F, blocks, i, j)
    starts := BlockStartPositions(blocks);
    dims := [BlockDimension(block) : block in blocks];
    if dims[i] ne dims[j] then
        error "Can only swap blocks of equal dimension";
    end if;

    n := TotalDimension(blocks);
    perm := [k : k in [1..n]];
    di := dims[i];
    dj := dims[j];
    _ := dj;
    for t in [0..di - 1] do
        perm[starts[i] + t] := starts[j] + t;
        perm[starts[j] + t] := starts[i] + t;
    end for;

    P := ZeroMatrix(F, n, n);
    for r in [1..n] do
        P[r, perm[r]] := F!1;
    end for;
    return P;
end function;

NSplitSwapLift := function(F, blocks)
    M := Matrix(F, 2, 2, [1, 0, 1, -1]);
    ia := 0;
    ib := 0;
    starts := BlockStartPositions(blocks);
    n := TotalDimension(blocks);
    H := IdentityMatrix(F, n);
    for i in [1..#blocks] do
        if blocks[i]`kind eq "J" and blocks[i]`label eq "a" and blocks[i]`d eq 1 then
            ia := i;
        elif blocks[i]`kind eq "J" and blocks[i]`label eq "b" and blocks[i]`d eq 1 then
            ib := i;
        elif not (blocks[i]`kind eq "S" and blocks[i]`d eq 1) then
            InsertBlock(~H, LocalBlockQuotientLift(F, blocks[i], M), starts[i], starts[i]);
        end if;
    end for;
    if ia eq 0 or ib eq 0 then
        error "Expected one a-simple block and one b-simple block";
    end if;
    return SwapTwoBlocksMatrix(F, blocks, ia, ib) * H;
end function;

S3Swap01Lift := function(F, blocks)
    M := Matrix(F, 2, 2, [1, 0, 1, -1]);
    ia := 0;
    ib := 0;
    ic := 0;
    starts := BlockStartPositions(blocks);
    n := TotalDimension(blocks);
    H := IdentityMatrix(F, n);
    for i in [1..#blocks] do
        if blocks[i]`kind eq "J" and blocks[i]`label eq "a" and blocks[i]`d eq 1 then
            ia := i;
        elif blocks[i]`kind eq "J" and blocks[i]`label eq "b" and blocks[i]`d eq 1 then
            ib := i;
        elif blocks[i]`kind eq "J" and blocks[i]`label eq "c" and blocks[i]`d eq 1 then
            ic := i;
            InsertBlock(~H, LocalBlockQuotientLift(F, blocks[i], M), starts[i], starts[i]);
        elif not (blocks[i]`kind eq "S" and blocks[i]`d eq 1) then
            InsertBlock(~H, LocalBlockQuotientLift(F, blocks[i], M), starts[i], starts[i]);
        end if;
    end for;
    _ := ic;
    if ia eq 0 or ib eq 0 then
        error "Expected one a-simple block and one b-simple block";
    end if;
    return SwapTwoBlocksMatrix(F, blocks, ia, ib) * H;
end function;

S3Swap0InfLift := function(F, blocks)
    ia := 0;
    ib := 0;
    ic := 0;
    starts := BlockStartPositions(blocks);
    n := TotalDimension(blocks);
    H := IdentityMatrix(F, n);

    for i in [1..#blocks] do
        if blocks[i]`kind eq "J" and blocks[i]`label eq "a" and blocks[i]`d eq 1 then
            ia := i;
        elif blocks[i]`kind eq "J" and blocks[i]`label eq "b" and blocks[i]`d eq 1 then
            ib := i;
        elif blocks[i]`kind eq "J" and blocks[i]`label eq "c" and blocks[i]`d eq 1 then
            ic := i;
        elif not (blocks[i]`kind eq "S" and blocks[i]`d eq 1) then
            InsertBlock(~H, IdentityMatrix(F, BlockDimension(blocks[i])), starts[i], starts[i]);
        end if;
    end for;

    if ib ne 0 then
        InsertBlock(~H, IdentityMatrix(F, 2), starts[ib], starts[ib]);
    end if;

    if ia eq 0 or ic eq 0 then
        error "Expected one a-simple block and one c-simple block";
    end if;

    return SwapTwoBlocksMatrix(F, blocks, ia, ic) * H;
end function;

CodomainLiftMatricesFromBlocks := function(F, blocks)
    qtype := FiniteQuotientStringFromBlocks(blocks);
    lifts := [];

    if #F gt 2 then
        w := PrimitiveElement(F);
        C := Matrix(F, 2, 2, [w, 0, 0, w]);
    else
        C := IdentityMatrix(F, 2);
    end if;

    if qtype eq "GL_2(q)" then
        U := Matrix(F, 2, 2, [1, 1, 0, 1]);
        W := Matrix(F, 2, 2, [0, -1, 1, 0]);
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, U));
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, W));
        if #F gt 2 then
            D := Matrix(F, 2, 2, [PrimitiveElement(F), 0, 0, 1]);
            Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, D));
        end if;
        return lifts;
    elif qtype eq "B(q)" then
        U := Matrix(F, 2, 2, [1, 1, 0, 1]);
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, U));
        if #F gt 2 then
            D := Matrix(F, 2, 2, [PrimitiveElement(F), 0, 0, 1]);
            Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, D));
            Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, C));
        end if;
        return lifts;
    elif qtype eq "T_s(q)" then
        if #F eq 2 then
            return [];
        end if;
        w := PrimitiveElement(F);
        T := Matrix(F, 2, 2, [1, w - 1, 0, w]);
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, T));
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, C));
        return lifts;
    elif qtype eq "N_s(q)" then
        if #F gt 2 then
            w := PrimitiveElement(F);
            T := Matrix(F, 2, 2, [1, w - 1, 0, w]);
            Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, T));
            Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, C));
        end if;
        Append(~lifts, NSplitSwapLift(F, blocks));
        return lifts;
    elif qtype eq "N_ns(q)" then
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, PrimitiveExtensionTorusGenerator(F, 2)));
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, ExtensionFrobeniusMatrix(F, 2)));
        return lifts;
    elif qtype eq "~C_2(q)" then
        M := MixedQuadraticInvolutionMatrix(F);
        if #F gt 2 then
            Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, C));
        end if;
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, M));
        return lifts;
    elif qtype eq "~C_3(q)" then
        M := CubicCycleMatrix(F);
        if #F gt 2 then
            Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, C));
        end if;
        Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, M));
        return lifts;
    elif qtype eq "~S_3(q)" then
        if #F gt 2 then
            Append(~lifts, DirectSumQuotientLiftFromMatrix(F, blocks, C));
        end if;
        Append(~lifts, S3Swap01Lift(F, blocks));
        Append(~lifts, S3Swap0InfLift(F, blocks));
        return lifts;
    end if;

    error "Unsupported codomain-lift family";
end function;

FullStabilizerGeneratorsFromBlocks := function(F, blocks)
    return [g : g in KernelGeneratorsFromBlocks(F, blocks)]
        cat [h : h in CodomainLiftMatricesFromBlocks(F, blocks)];
end function;

QuotientOrderFromBlocks := function(F, blocks)
    qtype := FiniteQuotientStringFromBlocks(blocks);
    q := #F;
    if qtype eq "GL_2(q)" then
        return #GL(2, F);
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

KernelOrderNoRadicals := function(F, blocks)
    q := #F;

    if &and[block`kind eq "S" : block in blocks] then
        if #blocks eq 1 then
            d := blocks[1]`d;
            return q^(2 * d - 2) * (q - 1);
        elif #blocks eq 2 and blocks[1]`d eq 2 and blocks[2]`d eq 2 then
            return q^6 * #GL(2, F);
        end if;
        error "Unsupported pure singular core";
    end if;

    nontrivial_s := [block : block in blocks | block`kind eq "S" and block`d gt 1];
    regular_core := [block : block in blocks | block`kind ne "S"];

    if #nontrivial_s eq 0 then
        clusters := RegularClustersWithOffsets(blocks);
        order := 1;
        for entry in clusters do
            cluster := entry[2];
            if cluster[1]`kind eq "P" then
                m := cluster[1]`degree;
                order *:= q^m * (q^(2 * m) - 1);
            else
                ds := [block`d : block in cluster];
                if #ds eq 1 then
                    d := ds[1];
                    order *:= q^(3 * d - 3) * q * (q^2 - 1);
                elif ds eq [1, 1] then
                    order *:= #SymplecticGroup(4, F);
                elif ds eq [2, 1] then
                    order *:= q^7 * (q * (q^2 - 1))^2;
                else
                    error "Unsupported regular cluster";
                end if;
            end if;
        end for;
        return order;
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 2 and #regular_core gt 0 then
        return q^(TotalDimension(regular_core) + 2) * (q - 1) * KernelOrderNoRadicals(F, regular_core);
    elif #nontrivial_s eq 1 and nontrivial_s[1]`d eq 3 and #regular_core eq 1 and
        regular_core[1]`kind eq "J" and regular_core[1]`d eq 1 then
        return q^6 * (q - 1) * q * (q^2 - 1);
    end if;

    error "Unsupported non-radical family";
end function;

KernelOrderFromBlocks := function(F, blocks)
    r, core := SplitOffS1Blocks(blocks);
    if r gt 0 then
        return (#F)^(r * TotalDimension(core)) * #GL(r, F) * KernelOrderNoRadicals(F, core);
    end if;
    return KernelOrderNoRadicals(F, blocks);
end function;

FiniteStabilizerSpecsNLe7 := function(F, n)
    data := FiniteOrbitDataNLe7(F, n);
    out := [* *];

    for reci in data do
        Append(~out, rec<FiniteStabilizerSpecRF |
            name := reci`name,
            blocks := reci`blocks,
            pair := reci`pair,
            wedge1 := reci`wedge1,
            wedge2 := reci`wedge2,
            kernel := reci`kernel,
            quotient := reci`quotient,
            kernel_generators := KernelGeneratorsFromBlocks(F, reci`blocks),
            quotient_generators := QuotientGeneratorsFromBlocks(F, reci`blocks),
            quotient_lifts := CodomainLiftMatricesFromBlocks(F, reci`blocks),
            full_generators := FullStabilizerGeneratorsFromBlocks(F, reci`blocks)
        >);
    end for;

    return out;
end function;

PrintFiniteStabilizerSpecs := procedure(F, n)
    specs := FiniteStabilizerSpecsNLe7(F, n);
    printf "n = %o, q = %o, finite-field types = %o\n", n, #F, #specs;
    for i in [1..#specs] do
        S := specs[i];
        printf "[%o] %o\n", i, S`name;
        printf "  omega_1 = %o\n", S`wedge1;
        printf "  omega_2 = %o\n", S`wedge2;
        printf "  K_L = %o\n", S`kernel;
        printf "  Q_L = %o\n", S`quotient;
        printf "  #kernel generators = %o\n", #S`kernel_generators;
        printf "  #quotient generators = %o\n", #S`quotient_generators;
        printf "  #quotient lifts = %o\n", #S`quotient_lifts;
        printf "  #full stabilizer generators = %o\n", #S`full_generators;
    end for;
end procedure;

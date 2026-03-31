load "src/paper_data.m";

HankelMatrix := function(F, d, coeffs)
    H := ZeroMatrix(F, d - 1, d);
    for r in [1..d - 1] do
        for s in [1..d] do
            H[r, s] := coeffs[r + s - 1];
        end for;
    end for;
    return H;
end function;

SingularKernelCandidateMatrix := function(F, d, a, coeffs)
    H := HankelMatrix(F, d, coeffs);
    top := HorizontalJoin(a * IdentityMatrix(F, d), ZeroMatrix(F, d, d - 1));
    bottom := HorizontalJoin(H, a^-1 * IdentityMatrix(F, d - 1));
    return VerticalJoin(top, bottom);
end function;

S22Pair := function(F)
    A := ZeroMatrix(F, 6, 6);
    B := ZeroMatrix(F, 6, 6);

    A[1, 5] := F!1; A[5, 1] := -F!1;
    A[2, 6] := F!1; A[6, 2] := -F!1;
    B[3, 5] := F!1; B[5, 3] := -F!1;
    B[4, 6] := F!1; B[6, 4] := -F!1;

    return <A, B>;
end function;

S22CandidateMatrix := function(F, T, C1, C2)
    Z := ZeroMatrix(F, 2, 2);
    Tinvt := Transpose(T^-1);
    row1 := HorizontalJoin(T, HorizontalJoin(Z, Z));
    row2 := HorizontalJoin(Z, HorizontalJoin(T, Z));
    row3 := HorizontalJoin(Tinvt * C1, HorizontalJoin(Tinvt * C2, Tinvt));
    return VerticalJoin(row1, VerticalJoin(row2, row3));
end function;

S3PlusJ1Pair := function(F)
    A := ZeroMatrix(F, 7, 7);
    B := ZeroMatrix(F, 7, 7);

    A[1, 4] := F!1; A[4, 1] := -F!1;
    A[2, 5] := F!1; A[5, 2] := -F!1;
    A[6, 7] := F!1; A[7, 6] := -F!1;
    B[2, 4] := F!1; B[4, 2] := -F!1;
    B[3, 5] := F!1; B[5, 3] := -F!1;

    return <A, B>;
end function;

S3PlusJ1UnipotentMatrix := function(F, x, y, z, t, r, s)
    return Matrix(F, 7, 7, [
        1, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0,
        x, y, z, 1, 0, r, s,
        y, z, t, 0, 1, 0, 0,
        -s, 0, 0, 0, 0, 1, 0,
        r, 0, 0, 0, 0, 0, 1
    ]);
end function;

S3PlusJ1LeviMatrix := function(F, a, h)
    G := IdentityMatrix(F, 7);
    for i in [1..3] do
        G[i, i] := a;
    end for;
    G[4, 4] := a^-1;
    G[5, 5] := a^-1;
    InsertBlock(~G, h, 6, 6);
    return G;
end function;

CheckSingularBlockKernelModel := procedure(q, d)
    F := GF(q);
    pair := SingularBlockPair(F, d);
    A := pair[1];
    B := pair[2];

    gens := [];
    for i in [1..2 * d - 2] do
        coeffs := [F | 0 : j in [1..2 * d - 2]];
        coeffs[i] := F!1;
        Append(~gens, SingularKernelCandidateMatrix(F, d, F!1, coeffs));
    end for;
    if q gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, SingularKernelCandidateMatrix(F, d, w,
            [F | 0 : i in [1..2 * d - 2]]));
    end if;

    G := sub<GL(2 * d - 1, F) | gens>;
    expected := q^(2 * d - 2) * (q - 1);
    assert &and[g * A * Transpose(g) eq A and g * B * Transpose(g) eq B : g in gens];
    assert #G eq expected;
    printf "S_%o over GF(%o): order %o verified\n", d, q, #G;
end procedure;

CheckS22KernelModel := procedure(q)
    F := GF(q);
    pair := S22Pair(F);
    A := pair[1];
    B := pair[2];
    I2 := IdentityMatrix(F, 2);
    Z := ZeroMatrix(F, 2, 2);

    sym_basis := [
        Matrix(F, 2, 2, [1, 0, 0, 0]),
        Matrix(F, 2, 2, [0, 1, 1, 0]),
        Matrix(F, 2, 2, [0, 0, 0, 1])
    ];

    gens := [];
    for C in sym_basis do
        Append(~gens, S22CandidateMatrix(F, I2, C, Z));
        Append(~gens, S22CandidateMatrix(F, I2, Z, C));
    end for;
    Append(~gens, S22CandidateMatrix(F, Matrix(F, 2, 2, [1, 1, 0, 1]), Z, Z));
    Append(~gens, S22CandidateMatrix(F, Matrix(F, 2, 2, [0, -1, 1, 0]), Z, Z));
    if q gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, S22CandidateMatrix(F, Matrix(F, 2, 2, [w, 0, 0, 1]), Z, Z));
    end if;

    G := sub<GL(6, F) | gens>;
    expected := q^6 * #GL(2, F);
    assert &and[g * A * Transpose(g) eq A and g * B * Transpose(g) eq B : g in gens];
    assert #G eq expected;
    printf "S_2^2 over GF(%o): order %o verified\n", q, #G;
end procedure;

CheckS3PlusJ1LocalModel := procedure(q)
    F := GF(q);
    pair := S3PlusJ1Pair(F);
    A := pair[1];
    B := pair[2];
    I2 := IdentityMatrix(F, 2);

    gens := [
        S3PlusJ1UnipotentMatrix(F, F!1, F!0, F!0, F!0, F!0, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!1, F!0, F!0, F!0, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!0, F!1, F!0, F!0, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!0, F!0, F!1, F!0, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!0, F!0, F!0, F!1, F!0),
        S3PlusJ1UnipotentMatrix(F, F!0, F!0, F!0, F!0, F!0, F!1),
        S3PlusJ1LeviMatrix(F, F!1, Matrix(F, 2, 2, [1, 1, 0, 1])),
        S3PlusJ1LeviMatrix(F, F!1, Matrix(F, 2, 2, [0, -1, 1, 0]))
    ];
    if q gt 2 then
        w := PrimitiveElement(F);
        Append(~gens, S3PlusJ1LeviMatrix(F, w, I2));
        Append(~gens, S3PlusJ1LeviMatrix(F, F!1, Matrix(F, 2, 2, [w, 0, 0, w^-1])));
    end if;

    G := sub<GL(7, F) | gens>;
    expected := q^6 * (q - 1) * q * (q^2 - 1);
    assert &and[g * A * Transpose(g) in sub<MatrixAlgebra(F, 7) | A, B> and
                g * B * Transpose(g) in sub<MatrixAlgebra(F, 7) | A, B> : g in gens];
    assert #G eq expected;
    printf "S_3 + J_{a,1} over GF(%o): order %o verified\n", q, #G;
end procedure;

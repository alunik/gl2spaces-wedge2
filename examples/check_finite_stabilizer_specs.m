load "src/stabilizer_generators.m";

VerifyFullField := procedure(q)
    F := GF(q);
    for n in [4..7] do
        specs := FiniteStabilizerSpecsNLe7(F, n);
        for S in specs do
            A := S`pair[1];
            B := S`pair[2];

            for g in S`kernel_generators do
                assert g * A * Transpose(g) eq A;
                assert g * B * Transpose(g) eq B;
            end for;

            K := sub<GL(Nrows(A), F) | S`kernel_generators>;
            assert #K eq KernelOrderFromBlocks(F, S`blocks);

            image_gens := [];
            for g in S`full_generators do
                M := ActionMatrixOnPair(S`pair, g);
                Append(~image_gens, M);
                assert g * A * Transpose(g) eq M[1, 1] * A + M[2, 1] * B;
                assert g * B * Transpose(g) eq M[1, 2] * A + M[2, 2] * B;
            end for;

            Qim := sub<GL(2, F) | image_gens>;
            assert #Qim eq QuotientOrderFromBlocks(F, S`blocks);

            G := sub<GL(Nrows(A), F) | S`full_generators>;
            assert #G eq #K * #Qim;
        end for;
        printf "q=%o n=%o verified %o full stabilizer specs\n", q, n, #specs;
    end for;
end procedure;

VerifyLiftImages := procedure(q)
    F := GF(q);
    for n in [4..7] do
        specs := FiniteStabilizerSpecsNLe7(F, n);
        for S in specs do
            mats := [ActionMatrixOnPair(S`pair, g) : g in S`full_generators];
            Qim := sub<GL(2, F) | mats>;
            assert #Qim eq QuotientOrderFromBlocks(F, S`blocks);
        end for;
        printf "q=%o n=%o verified %o quotient lift images\n", q, n, #specs;
    end for;
end procedure;

for q in [2, 3] do
    VerifyFullField(q);
end for;

for q in [5, 7] do
    VerifyLiftImages(q);
end for;

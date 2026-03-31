load "src/paper_data.m";
load "src/local_matrix_models.m";

expected_geometric := AssociativeArray();
expected_geometric[4] := 3; expected_geometric[5] := 5;
expected_geometric[6] := 11; expected_geometric[7] := 16;

expected_finite := AssociativeArray();
expected_finite[4] := 4; expected_finite[5] := 6;
expected_finite[6] := 14; expected_finite[7] := 20;

F := GF(5);
for n in [4..7] do
    assert #GeometricOrbitDataNLe7(F, n) eq expected_geometric[n];
    assert #FiniteOrbitDataNLe7(F, n) eq expected_finite[n];
    printf "Counts verified for n = %o\n", n;
end for;

for q in [2, 3, 5] do
    for d in [2, 3, 4] do
        CheckSingularBlockKernelModel(q, d);
    end for;
    CheckS22KernelModel(q);
end for;

for q in [2, 3, 4, 5] do
    CheckS3PlusJ1LocalModel(q);
end for;

load "src/paper_data.m";

F := GF(5);
for n in [4..7] do
    data := GeometricOrbitDataNLe7(F, n);
    PrintPaperOrbitData(data : Title := Sprintf("Geometric data for n = %o over GF(%o)", n, #F));
    print "";
end for;

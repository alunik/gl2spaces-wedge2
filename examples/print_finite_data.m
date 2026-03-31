load "src/paper_data.m";

F := GF(5);
for n in [4..7] do
    data := FiniteOrbitDataNLe7(F, n);
    PrintPaperOrbitData(data : Title := Sprintf("Finite-field data for n = %o over GF(%o)", n, #F));
    print "";
end for;

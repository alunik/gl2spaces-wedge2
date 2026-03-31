load "src/stabilizer_generators.m";

q := 3;
n := 6;

F := GF(q);
specs := FiniteStabilizerSpecsNLe7(F, n);

printf "q = %o, n = %o\n", q, n;
for i in [1..#specs] do
    S := specs[i];
    printf "[%o] %o\n", i, S`name;
    printf "  omega_1 = %o\n", S`wedge1;
    printf "  omega_2 = %o\n", S`wedge2;
    printf "  K_L = %o\n", S`kernel;
    printf "  Q_L = %o\n", S`quotient;
    printf "  full stabilizer generators:\n";
    for g in S`full_generators do
        print g;
    end for;
end for;

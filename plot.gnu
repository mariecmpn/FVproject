# plot GNU pour le transport scalaire

set yrange [0:4]
set title "Solution à T_{fin} du modèle scalaire"
plot 'solution.dat' with lines, 'solution_ex.dat' with lines

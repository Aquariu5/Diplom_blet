filename(n) = sprintf("fracture%i.txt", n)
set terminal gif animate delay 2
set output "anim_fracture.gif"
set title "Fracture distrib" font "Times New Roman, 12"
do for[ii=0:5000] {
	set zrange[0:1]
	set ticslevel 0
	set view 0, 90, 1, 1
	splot filename(ii) matrix w pm3d
}
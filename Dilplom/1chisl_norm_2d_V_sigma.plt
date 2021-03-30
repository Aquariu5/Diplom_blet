filename(n) = sprintf("ch_norm_2d_sigma_V%i.txt", n)
set terminal gif animate delay 8
set output "anim_2d_norm_sigma_V.gif"
set title "V numerical" font "Times New Roman, 12"
do for[ii=0:500] {
	unset key
	set view 0, 90, 1, 1
	set zrange[-1000:1000]
	set ticslevel 0
	splot filename(ii) matrix w pm3d
}
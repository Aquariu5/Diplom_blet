filename(n) = sprintf("ch_norm_2d_sigma_P%i.txt", n)
set terminal gif animate delay 10
set output "anim_2d_norm_sigma_P.gif"
set title "P numerical" font "Times New Roman, 12"
do for[ii=0:150] {
	unset key
	#set zrange[0:500]
	set zrange[-50000:50000]
	set view 0, 90, 1, 1
	set ticslevel 0
	splot filename(ii) matrix w pm3d
}
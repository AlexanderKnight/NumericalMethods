set term png size 1080, 1080
set output "3.a.png"
set title "dx vs log(f^~-f), f(x)=(x^2-1)/(x+2)"
a = 1.1
f(x) = a*(x-2.3)+6.2
plot '3.a.gen.dat' u 1:2, f(x) title 'x^{1.1}';

clear

set term png size 1080, 1080
set output "3.b.png"
set title "dx vs log(f^~-f), f(x)=H(x-0.15)"
plot '3.b.gen.dat' u 1:2

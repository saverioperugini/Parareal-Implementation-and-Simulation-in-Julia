PARAREAL ALGORITHM IMPLEMENTATION AND SIMULATION IN JULIA

Tyler M. Masthay and Saverio Perugini 
Department of Computer Science 
University of Dayton 
300 College Park 
Dayton, Ohio  45469–2160 
(937) 229–4079 
{tmasthay1,saverio}@udayton.edu

ABSTRACT

We present a full implementation of the parareal algorithm—an integration
technique to solve diﬀerential equations in parallel—in the Julia programming
language for a fully general, ﬁrst-order, initial-value problem. We provide a
brief overview of Julia—a concurrent programming language for scientiﬁc
computing. Our implementation of the parareal algorithm accepts both coarse and
ﬁne integrators as functional arguments. We use Euler’s method and another
Runge-Kutta integration technique as the integrators in our experiments. We
also present a simulation of the algorithm for purposes of pedagogy and as a
tool for investigating the performance of the algorithm.

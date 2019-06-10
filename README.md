PARAREAL ALGORITHM IMPLEMENTATION AND SIMULATION IN JULIA

Tyler M. Masthay and Saverio Perugini

Institute for Computational and Engineering Sciences, University of Texas Austin, tyler@ices.utexas.edu

Department of Computer Science, University of Dayton, 300 College Park Dayton, Ohio  45469–2160, saverio@udayton.edu

# Summary

We present an implementation of the parareal algorithm—an integration technique
to solve diﬀerential equations in parallel—ﬁrst proposed in 2001 by Lions,
Maday, and Turinici [@Lions:2001]—in the Julia programming language
[@Julia:2014] for a fully general, ﬁrst-order, initial-value problem.  We also
provide a graphical simulation of the parareal algorithm intended to be used in
a numerical analysis course to both introduce the parareal algorithm to
students and aid them in an investigation of the types of curves for which the
parareal algorithm might be practical.  Our implementation of the parareal
algorithm accepts both coarse and ﬁne integrators as functional arguments.  We
provide implementations of Euler’s method and another Runge-Kutta integration
technique as the integrators.  The final two functions in the source code file
\texttt{parareal.jl} are functions implementing Euler’s method and another
Runge-Kutta integration technique that can be used as examples to be passed as
ﬁrst-class functions as coarse or ﬁne integration techniques to the
\texttt{parareal} or \texttt{simulate} functions.  A Git repository of both the
implementation and graphical simulation is available in GitHub at
<https://github.com/sperugin/Parareal-Implementation-and-Simulation-in-Julia.git>.
All of the graphical plots are generated with the Julia Plots package available
at <https://juliaplots.github.io/>.  A video describing this application of
Julia is available on YouTube at <https://www.youtube.com/watch?v=MtgbeLO6ZM4>.
The purpose of this software is pedagogical: as a simulation to introduce
students to the parareal algorithm and the concept of concurrency, and as a
tool for (graphically) investigating the performance of the algorithm.

Language: Julia
Version: 0.5.1
Compiler: Julia Built-in
Interpreter: Julia Built-in
Test Platform: Windows (64-bit)
Download page: http://julia-lang.org/downloads/

parareal.jl: accepts any fine and coarse integrator (e.g., it can accept any
Runge-Kutta integration technique).

# SuperCoV
Simulation code for superspreading of SARS-CoV-2, written in FORTRAN.

Instructions:

The program runs with parameters hardcoded in the source code,
the parameters can be found under the comment which begins with
c.....Normal parameters

Parameter combination for non-superspreader simulation:
Nsupermax=1
superfreq=0
beta0=1.1
helpfac=1.

Parameter combination for superspreader simulation:
beta0=1
Nsuperlim=5000
Superpow=0.9    !(dispersion=0.1)
Rep0=1.1
helpfac=10

Mitigation is controlled later in program by fudges0
fudges0=1  Normal unmitigated epidemic
fudges0=0  Close down random sector when 1% attack rate is reached. 

Can be compiled with the GNU Fortan compiler (gfortran) by executing the following command:
gfortran -O cs.f -o c
The program may then be executed by running:
./c>c.d
with the number of infected people during the epidemic being recorded in column 3 of the data file c.d, and the 1st column containing the time variable.


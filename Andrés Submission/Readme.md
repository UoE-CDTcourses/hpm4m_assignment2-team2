# Results

As we don't have an HPC right now, I used another machine with 20 cores. For the serial code, I made the following experiments while changing the value of `dt` so that we could study how the Mean Average Error (MAE) evolves.

| $\frac{\mathrm{d}t}{\mathrm{d}x \times \mathrm{d}x}$ | `MAE` |
| -----:|-----:|
|0.000001 | 4.42154 |
|0.00001 | 4.12776 |
|0.0001 | 2.07546 |
|0.001 | 0.00214094 |
|0.01 | 3.28387e-07 |
|0.05 | 1.975e-09 |
|0.1 | 4.98021e-11 |
|0.5 | 5.91209e-22 |
|0.505 | 2.46342e+141 |
|0.51 | `inf` |
|0.75 | `nan` |

As we can see, a value between $(10^{-3}, 0.5)$ yields good numerical results. 

The parallel code is in `heat.cpp`. The code computes the values of $M$ and $N$ to achieve a ratio  $\frac{\mathrm{d}t}{\mathrm{d}x \times \mathrm{d}x} = \frac{1}{3}$, where the user can just change $J$. The quantities are computed as follows:
\begin{align*}
    M &= \mathtt{size} \times (J - 2) + 2;
    \\
    N &= \lceil 3 \times T \times M \times M \rceil 
\end{align*}
These come from the fact that $\Delta t$ is equal to $T/N$ and $\Delta x$ is equal to $1/M$.

The problem with the above scheme, which I realise now after spending so much time on it, is that as we select $T = 0.1$, the value of $N$ is too big in comparison to the value of $M$. In that sense, we are not taking advantage of the parallelism as looping $N$ steps in parallel takes too much time. Anyways, I believe I could flip the scheme in the future and see how it goes.

As you might guess, the code takes a while to excecute. The best `MAE` at the time of this report is  , which resulted after calling 
```
    mpirun -np 8 ./Heat 0.1
```
The output reads as follows:
```
    dx = 0.000125219, dt = 5.22661e-09, dt/dx² = 0.333333

    True and numerical values at M=7986 space points, N=19132859 time points, at time T=0.1:

    True values           Numerical solutions

    Printing first 5 values
    0                      0
    1.51818e-05            2.61662e-30
    3.03637e-05            5.23321e-30
    4.55455e-05            7.84975e-30
    6.07273e-05            1.04662e-29
    Printing last 5 values
    -6.07273e-05            -0.0818092
    -4.55455e-05            -0.0613619
    -3.03637e-05            -0.0409104
    -1.51818e-05            -0.0204559
    -4.72623e-18            0
    
    MAE: 2.0834

    
```

Something good to note here is that the code also allows the possibility of printing out the solution at an early stage. Then I checked whether the solution was actually evolving in the way that it should. If we iterate for the first `1000` time-steps, `MAE` is less than `0.01`.  After a while, the resulting scheme is usually poor in terms of $\Delta x$, as the number of spatial terms is small.

Another idea that might come up here is to repeat the process starting at other times:
```
    mpirun -np 10 ./Heat .0001

    dx = 0.00010018, dt = 3.34526e-09, dt/dx² = 0.333323

    True and numerical values at M=9982 space points, N=29893 time points, at time T=0.0001:

    True values           Numerical solutions

Printing first 5 values
    0            0
    0.0164217            0.0164216
    0.0328428            0.0328428
    0.0492629            0.0492629
    0.0656815            0.0656814
    Printing last 5 values
    -0.0411172            -0.0654551
    -0.0308397            -0.0490939
    -0.0205606            -0.0327305
    -0.0102806            -0.0163656
    6.89697e-15            0
    MAE: 0.363139
```

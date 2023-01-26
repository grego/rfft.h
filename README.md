# rfft.h
Public domain single header fast Fourier transform for arbitrary array sizes,
in about 100 lines of C code, which should be straightforward to understand.

## Algorithms
The classic [Cooley-Turkey algorithm](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)
is used in place (without additional allocations) for arrays of size `2^k`.

For more more general ones, [Bluestein algorithm](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm)
is used. It utilizes the binomial identity `2nk = n^2 + k^2 - (k - n)^2` to
express the Fourier transform as a convolution of two sequences,
which can be computed using the algorith for the power of 2 sizes.
It needs to allocate two auxillary array of size at most `4n + 3`.

The runnig time is always `O(n log(n))`. However, if the speed is crucial,
more optimized libraries like [FFTW](http://fftw.org/) are recommended.

## Inspiration
Inspired by [Project Nayuki](https://www.nayuki.io/page/free-small-fft-in-multiple-languages),
but written in a simpler and arguably more straightforward way.

## License
Public domain.

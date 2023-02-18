#ifndef RFFT_HPP
#define RFFT_HPP
// Reasonably Fast Fourier Transform (in public domain)

// Perform the FFT in place for an array of size 2^k.
// No normalization is done.
void fft_transform_radix2(std::vector<std::complex<double>>& vec, bool inverse);

// Perform the FFT for an arbitrary array using the Bluestein's algorithm.
// The code needs to allocate two supplementary buffers (of size < 4 * n - 3).
// No normalization is done.
void fft_transform_bluestein(std::vector<std::complex<double>>& vec, bool inverse);

// Perform the FFT, choosing the suitable algorithm from the two above.
void fft_transform(std::vector<std::complex<double>>& vec, bool inverse);

#endif // RFFT_HPP

#ifdef RFFT_IMPLEMENTATION
#include <complex>

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif 

void fft_transform_radix2(std::vector<std::complex<double>>& vec, bool inverse) {
	size_t n = vec.size();
	int levels = 0;	 // Compute levels = floor(log2(n))
	for (size_t k =	1; (k &	n) == 0; k <<= 1)
		levels++;
	
	// Permute vec by reversing the bits of addresses
	for (size_t i =	0; i < n; i++) {
		// Reverse the bits of i
		size_t j = 0, ii = i;
		for (int k = 0;	k < levels; k++) {
			j = (j << 1) | (ii & 1);
			ii >>= 1;
		}
	
		if (j >	i) {
			std::complex<double> tmp = vec[i];
			vec[i] = vec[j];
			vec[j] = tmp;
		}
	}
	
	// Cooley-Tukey	in place
	for (size_t half = 1; half < n;	half *=	2) {
		size_t size = 2	* half;
		size_t step = n	/ size;
		for (size_t j =	0, k = 0; j < half; j++, k += step) {
			double angle = (inverse	? 2 : -2) * M_PI * k / n;
			std::complex<double> omega = std::polar(1.0, angle);
        		for (size_t i =	0; i < n; i += size) {
				std::complex<double> tmp = vec[i + j + half] * omega;
				vec[i + j + half] = vec[i + j] - tmp;
				vec[i + j] += tmp;
			}
		}
	}
}

void fft_transform_bluestein(std::vector<std::complex<double>>& vec, bool inverse) {
	size_t n = vec.size();
	// Find m = 2^k such that m >= 2 * n + 1
	size_t m = 1;
	while (m <= 2 * n) {
		m *= 2;
	}
	
	// Extended vectors of size m
	std::vector<std::complex<double>> avec(m), bvec(m);
	avec[0]	= vec[0];
	bvec[0]	= 1.0;
	for (size_t i =	1; i < n; i++) {
		size_t k = (i * i) % (2 * n);
		double angle = (inverse	? M_PI : -M_PI)	* k / n;
		std::complex<double> omega = std::polar(1.0, angle);
		avec[i]	= vec[i] * omega;
		bvec[i]	= bvec[m - i] =	std::conj(omega);
	}

	// Convolution
	fft_transform_radix2(avec, false);
	fft_transform_radix2(bvec, false);
	for (size_t i =	0; i < m; i++) {
		avec[i]	*= bvec[i];
	}
	fft_transform_radix2(avec, true);

	for (size_t i =	0; i < n; i++) {
		size_t k = (i * i) % (2 * n);
		double angle = (inverse	? M_PI : -M_PI)	* k / n;
		std::complex<double> omega = std::polar(1.0, angle);
		vec[i] = avec[i] * omega / double(m);
	}
}

void fft_transform(std::vector<std::complex<double>>& vec, bool inverse) {
	size_t n = vec.size();
	if (n <= 1)
		return;
	else if	((n & (n - 1)) == 0)  // Power of 2
		fft_transform_radix2(vec, inverse);
	else
		fft_transform_bluestein(vec, inverse);
}

#endif // RFFT_IMPLEMENTATION

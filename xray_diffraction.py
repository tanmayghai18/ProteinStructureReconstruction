# xray pseudocode

# we did not have time to do the full x-ray code
# instead, we made pseudocode for a related problem:
# Fraunhofer diffraction

# suppose we have a 1D slit of width a, and we send through it
# a plane wave with amplitude E_o traveling perpendicular
# what does the diffraction pattern look like?

# think of a plane wave as a wave whose wavefronts are represented by
# parallel planes

# http://kmdouglass.github.io/posts/approximating-diffraction-patterns-of-rectangular-apertures-with-the-fft.html

# we want to find the "irradiance," or
# the power received by a surface

# the irradiance is proportional to U(x, z)
# which is the amplitude of the field
# because we are working in 1D, x is the horizontal dimension
# z is how far it is

# the amplitude U(x, z) can be found as
# a convolution
# by convolution theorem, we can instead
# find two Fourier transforms
# in our case, we find the FT of
# our aperture function

width = a

# our "aperture function"
# this represents our plane wave through samples
Gin = [array of values such that -width/2 < x < width/2 = 1, all others 0 ]
# it looks like a piecewise function with all 0s except in the center where it's 1
# this is like rect

dx = sampling period of field
diffractedField = dx * fft(fftshift(field)) # an FT, it should resemble sinc

# we have to shift back
Gout(x, z) = U(x, z) = fftshift(diffractedField)

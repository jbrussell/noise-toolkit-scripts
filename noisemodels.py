import numpy as np

def accNLNM(f):
    """The Peterson New Low-Noise Model.

    Returns an acceleration ASD.

    """
    Pl = np.array([1.00e-02, 1.00e-01, 1.70e-01, 4.00e-01, 8.00e-01, 1.24e+00,
       2.40e+00, 4.30e+00, 5.00e+00, 6.00e+00, 1.00e+01, 1.20e+01,
       1.56e+01, 2.19e+01, 3.16e+01, 4.50e+01, 7.00e+01, 1.01e+02,
       1.54e+02, 3.28e+02, 6.00e+02, 1.00e+04])
    Al = np.array([-156.72, -162.36, -166.7 , -170.  , -166.4 , -168.6 , -159.98,
       -141.1 ,  -71.36,  -97.26, -132.18, -205.27,  -37.65, -114.37,
       -160.58, -187.5 , -216.47, -185.  , -168.34, -217.43, -258.28,
       -346.88])
    Bl = np.array([   5.64,    5.64,    0.  ,   -8.3 ,   28.9 ,   52.48,   29.81,
          0.  ,  -99.77,  -66.49,  -31.57,   36.16, -104.33,  -47.1 ,
        -16.28,    0.  ,   15.7 ,    0.  ,   -7.61,   11.9 ,   26.6 ,
         48.75])
    nlnm_dB = np.interp(1/f, Pl, Al+Bl*np.log10(Pl))
    nlnm = 10**(nlnm_dB/10) # convert dB rel 1 m**2/s**4/Hz --> m**2/s**4/Hz (not log)
    # nlnm = 10**(nlnm/20) # convert dB rel 1 m**2/s**4/Hz --> m/s**2/Hz**0.5 (not log)
        
    return nlnm



def accNHNM(f):
    """The Peterson New High-Noise Model.

    Returns an acceleration ASD.

    """
    Pl = np.array([0.1, 0.22, 0.32, 0.8, 3.8, 4.6, 6.3, 7.9, 15.4, 20.0, 354.8])
    Al = np.array([-108.73, -150.34, -122.31, -116.85, -108.48, -74.66, 0.66, -93.37, 73.54, -151.52, -206.66])
    Bl = np.array([-17.23, -80.5, -23.87, 32.51, 18.08, -32.95, -127.18, -22.42, -162.98, 10.01, 31.63])
    nhnm_dB = np.interp(1/f, Pl, Al+Bl*np.log10(Pl))
    nhnm = 10**(nhnm_dB/10) # convert dB rel 1 m**2/s**4/Hz --> m**2/s**4/Hz (not log)
    # nhnm = 10**(nhnm/20) # convert dB rel 1 m**2/s**4/Hz --> m/s**2/Hz**0.5 (not log)
    
    return nhnm
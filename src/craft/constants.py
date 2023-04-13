DM_CONSTANT = 4.15
MHZ_TO_GHZ = 1e-3
GHZ_TO_MHZ = 1e3
S_TO_MS = 1e3
MS_TO_S = 1e-3
S_TO_US = 1e6
US_TO_S = 1e-6
S_IN_A_DAY = 86400
C = 299792458
PI = 3.14159265358979323846
BOLTZMAN_CONSTANT = 1.380649 * 1e-23

def convert_hms_to_hours(hms):
    '''
    Converts time in hours:minutes:seconds.ff format to hours.ff format

    Params
    ------
    hms : str
        A string containing the hours:minutes:seconds.ff

    Returns
    -------
    hours : float
        A float containing the hours in hours.ff format

    '''

    hours, minutes, seconds = hms.strip().split(":")

    hours = float(hours) + float(minutes)/60 + float(seconds)/3600

    return hours



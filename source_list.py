"""
Calculate angular distances between sources in a WSClean sources file
and a target RA/DEC position using astropy.
"""

from pathlib import Path
from astropy.coordinates import SkyCoord, EarthLocation, get_body
from astropy.time import Time
import astropy.units as u
from casatools import table

src_coord_lst = {
    'CygA':  ['19:59:28.3', '+40.44.02'],
    'CasA':  ['-00:36:20.6', '+58.48.54'],
    'VirA':  ['12:30:49.4', '+12.23.28'],
    'TauA':  ['05:34:31.9', '+22.00.52'],
    'HerA':  ['16:51:08.1', '+04.59.34'],
    'HydA':  ['09:18:05.7', '-12.05.44']}

def parse_wsclean_coordinates(ra_str, dec_str):
    """Parse WSClean coordinates to SkyCoord object."""
    # Convert DEC format DD.MM.SS to DD:MM:SS for astropy
    dec_str_astropy = dec_str.replace('.', ':', 2)
    return SkyCoord(ra_str, dec_str_astropy, unit=(u.hourangle, u.deg))


def load_wsclean_sources(filename):
    """Load sources from WSClean sources file."""
    sources = []
    with open(filename, 'r') as f:
        for line in f.readlines()[1:]:  # Skip header
            parts = line.strip().split(',')
            if len(parts) < 4:
                continue
            
            name, source_type, ra_str, dec_str = parts[:4]
            flux = float(parts[4]) if len(parts) > 4 else 0.0
            
            try:
                coord = parse_wsclean_coordinates(ra_str, dec_str)
                sources.append({
                    'name': name, 'coord': coord, 'flux': flux,
                    'ra_deg': coord.ra.deg, 'dec_deg': coord.dec.deg
                })
            except (ValueError, IndexError):
                continue
    return sources


def distance_to_src_list(sourcelist_fname, ra_deg, dec_deg):
    """Calculate distances from sources to target position.
    
    Args:
        sourcelist_fname: Path to WSClean sources file
        ra_deg, dec_deg: Target coordinates in degrees
    
    Returns:
        List of dicts with source info and distances
    """
    sourcelist_file = Path(sourcelist_fname)
    if not sourcelist_file.exists():
        raise FileNotFoundError(f"Sources file {sourcelist_file} not found")
    
    target_coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
    sources = load_wsclean_sources(sourcelist_file)
    
    results = []
    for source in sources:
        sep = source['coord'].separation(target_coord)
        results.append({
            **source,  # Include all source data
            'distance_deg': sep.deg
        })
    
    return results


def get_time_mjd(msname):
    """Get the time in Modified Julian Days from a MS file."""
    tb = table()
    tb.open(msname+'/OBSERVATION')
    start_mjd = tb.getcol('TIME_RANGE')[0][0] / 86400.0  
    tb.close()
    return start_mjd



def get_Sun_RA_DEC(time_mjd, observatory='OVRO'):
    """Get the RA and DEC of the Sun at a given time.
    
    Args:
        time_mjd (float): Time in Modified Julian Days
        observatory (str): Observatory name for astropy EarthLocation
    
    Returns:
        tuple: (ra_deg, dec_deg) Sun coordinates in degrees
    
    Example:
        >>> ra, dec = get_Sun_RA_DEC(59000.5)  # MJD
        >>> print(f"Sun at RA={ra:.4f}°, DEC={dec:.4f}°")
    """
    # Convert MJD to astropy Time object
    time = Time(time_mjd, format='mjd')
    
    # Get observatory location
    location = EarthLocation.of_site(observatory)
    
    # Get Sun position as seen from the observatory
    sun_coord = get_body('sun', time, location)
    
    ra_hourangle = sun_coord.ra.degree
    dec_hourangle = sun_coord.dec.degree
    
    return str(ra_hourangle)+'deg', str(dec_hourangle)+'deg'

def mask_far_Sun_sources(sourcelis_fname, fname_out, ra_deg, dec_deg, distance_deg=8.0):
    """
    read in a sourcelist file, generate a new sourcelist txt file with only the
    source NOT within distance_deg from the Sun.
    """

    dist_to_sun = distance_to_src_list(sourcelis_fname, ra_deg, dec_deg)
    sources_to_remove = [source['name'] for source in dist_to_sun if source['distance_deg'] <= distance_deg]
    
    with open(sourcelis_fname, 'r') as f:
        lines = f.readlines()

    # remove the sources from the output file
    with open(fname_out, 'w') as f:
        for i, line in enumerate(lines):
            # Always include the header line (first line)
            if i == 0:
                f.write(line)
                continue
                
            name = line.split(',')[0]
            if name not in sources_to_remove:
                f.write(line)
                
    return fname_out
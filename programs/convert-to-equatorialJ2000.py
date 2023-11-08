import astropy.coordinates as coord
import astropy.units as u
import argparse

parser = argparse.ArgumentParser(
    description="""Plot side-by-side RGB images of sources""")


parser.add_argument("--ra", type=float, default=None,
                    help="""RA in degrees""")
parser.add_argument("--dec", type=float, default=None,
                    help="""DEC in degrees""")

cmd_args = parser.parse_args()

# Example coordinates (replace with your actual values)
ra_degrees =  cmd_args.ra # Replace with your RA in degrees
dec_degrees = cmd_args.dec  # Replace with your DEC in degrees

# Calculate RA and format it
ra = coord.Angle(ra_degrees, unit=u.deg).to_string(u.hour, sep=':', precision=2, pad=True)

# Calculate DEC and format it
dec = coord.Angle(dec_degrees, unit=u.deg).to_string(sep=':', precision=2, pad=True)

# Print the formatted coordinates
print(f'RA: {ra}')
print(f'DEC: {dec}')
